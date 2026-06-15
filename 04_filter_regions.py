#!/usr/bin/env python3
# =============================================================================
# Step 4 — Filter Lethal Candidates by Genomic Region
# =============================================================================
# Takes the lethal candidate sites from Step 3 and filters them to only those
# overlapping user-defined genomic regions (e.g. exons, CDS, promoters).
#
# Supports BED, GFF, GFF3, and GTF input formats — detected automatically
# from the file extension. Both compressed (.gz) and uncompressed files
# are supported.
#
# For GFF/GFF3/GTF files, only features matching --feature are used
# (default: exon). This prevents duplicate matches when a position overlaps
# multiple annotation types (gene, exon, CDS) for the same region.
#
# Parallelized — the candidates file is split into N chunks and each chunk
# is processed simultaneously. The regions dictionary is built once and
# shared across all workers (read-only, no memory multiplication).
#
# This script is called by the LethAl Finder master script and should not
# be run directly unless you know what you are doing.
#
# Usage:
#   python3 04_filter_regions.py <threshold> <output_dir> <regions_file>
#                                [--feature <type>] [--cores <n>]
#   e.g. python3 04_filter_regions.py 0.25 /results genome.gff --feature exon
#        python3 04_filter_regions.py 0.25 /results regions.bed --cores 4
# =============================================================================

import sys
import os
import gzip
import datetime
import argparse
import math
from multiprocessing import Pool


# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Filter lethal candidates by genomic region."
)
parser.add_argument("threshold",     type=str, help="Frequency threshold (e.g. 0.25)")
parser.add_argument("output_dir",    type=str, help="Directory containing lethal candidates file")
parser.add_argument("regions_file",  type=str, help="BED, GFF, GFF3, or GTF regions file")
parser.add_argument("--feature",     type=str, default="exon",
                    help="Feature type to extract from GFF/GTF files (default: exon). "
                         "Ignored for BED files.")
parser.add_argument("--cores",       type=int, default=1,
                    help="Number of CPU cores for parallelization (default: 1)")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

THRESHOLD      = str(float(args.threshold))
OUTPUT_DIR     = args.output_dir
REGIONS_FILE   = args.regions_file
FEATURE        = args.feature
CORES          = args.cores

LOG_DIR        = os.path.join(OUTPUT_DIR, "logs")
CANDIDATES_IN  = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}.tsv")
REGIONAL_OUT   = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}_regional.tsv")
TEMP_DIR       = os.path.join(OUTPUT_DIR, "temp_region_chunks")
CHUNK_SIZE     = 100_000  # rows per chunk per worker


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

def log(msg: str):
    """Print message with timestamp."""
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")
    sys.stdout.flush()


def validate_paths():
    """Check all required input paths exist before processing begins."""
    if not os.path.isfile(CANDIDATES_IN):
        print(f"[ERROR] Lethal candidates file not found: {CANDIDATES_IN}")
        print(f"        Make sure Step 3 has completed successfully.")
        sys.exit(1)
    if not os.path.isfile(REGIONS_FILE):
        print(f"[ERROR] Regions file not found: {REGIONS_FILE}")
        sys.exit(1)
    if not os.path.isdir(OUTPUT_DIR):
        print(f"[ERROR] Output directory not found: {OUTPUT_DIR}")
        sys.exit(1)


def detect_format(file_path: str) -> str:
    """
    Detect the regions file format from its extension.
    Handles both compressed (.gz) and uncompressed files.
    Returns one of: 'bed', 'gff', 'gtf'
    Exits with error if the extension is not recognized.
    """
    # Strip .gz before checking format extension
    ext = file_path.lower()
    if ext.endswith(".gz"):
        ext = ext[:-3]

    if ext.endswith(".bed"):
        return "bed"
    elif ext.endswith(".gff") or ext.endswith(".gff3"):
        return "gff"
    elif ext.endswith(".gtf"):
        return "gtf"
    else:
        print(f"[ERROR] Unrecognized regions file format: {file_path}")
        print(f"        Supported formats: .bed, .gff, .gff3, .gtf, and .gz compressed versions")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Regions File Parsing
# ---------------------------------------------------------------------------

def parse_bed(file_path: str) -> dict:
    """
    Parse a BED file and return a dictionary of intervals per chromosome.

    BED format (0-based, half-open intervals):
        col 0: chromosome
        col 1: start (0-based)
        col 2: end (exclusive)

    Converts to 1-based closed intervals internally for consistent comparison
    with genomic positions from HEARTY files.

    Returns:
        {chr: [(start, end), ...]}
    """
    regions = {}
    skipped = 0

    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            # Skip comment and track lines
            if not line or line.startswith("#") or line.startswith("track") \
                    or line.startswith("browser"):
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                skipped += 1
                continue

            chrom = parts[0]
            try:
                # Convert BED 0-based to 1-based
                start = int(parts[1]) + 1
                end   = int(parts[2])
            except ValueError:
                skipped += 1
                continue

            if chrom not in regions:
                regions[chrom] = []
            regions[chrom].append((start, end))

    if skipped > 0:
        log(f"  [WARNING] Skipped {skipped} malformed BED lines.")

    return regions


def parse_gff_gtf(file_path: str, feature: str) -> dict:
    """
    Parse a GFF/GFF3/GTF file and return a dictionary of intervals per
    chromosome, filtered to only the specified feature type.

    GFF/GTF format (1-based, closed intervals):
        col 0: chromosome
        col 2: feature type (e.g. gene, exon, CDS, mRNA)
        col 3: start (1-based)
        col 4: end (inclusive)

    Lines starting with # are comments and are skipped.

    Returns:
        {chr: [(start, end), ...]}
    """
    regions = {}
    skipped = 0
    matched = 0

    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            # Skip comment lines
            if not line or line.startswith("#"):
                continue

            parts = line.split('\t')
            if len(parts) < 5:
                skipped += 1
                continue

            chrom        = parts[0]
            feature_type = parts[2]

            # Only keep lines matching the requested feature type
            if feature_type != feature:
                continue

            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                skipped += 1
                continue

            if chrom not in regions:
                regions[chrom] = []
            regions[chrom].append((start, end))
            matched += 1

    if matched == 0:
        print(f"[ERROR] No features of type '{feature}' found in: {file_path}")
        if file_path.endswith('.gz'):
            print(f"        Check available feature types with: zcat {file_path} | cut -f3 | sort -u")
        else:
            print(f"        Check available feature types with: cut -f3 {file_path} | sort -u")
        sys.exit(1)

    if skipped > 0:
        log(f"  [WARNING] Skipped {skipped} malformed GFF/GTF lines.")

    return regions


def load_regions() -> dict:
    """
    Load and parse the regions file based on detected format.
    Returns a dictionary of intervals per chromosome.
    """
    fmt = detect_format(REGIONS_FILE)
    log(f"Detected format: {fmt.upper()}")

    if fmt == "bed":
        regions = parse_bed(REGIONS_FILE)
    else:
        log(f"Extracting feature type: '{FEATURE}'")
        regions = parse_gff_gtf(REGIONS_FILE, FEATURE)

    total_intervals = sum(len(v) for v in regions.values())
    log(f"  -> {total_intervals:,} intervals loaded across "
        f"{len(regions):,} chromosomes.")
    return regions


# ---------------------------------------------------------------------------
# Interval Overlap Check
# ---------------------------------------------------------------------------

def overlaps_any_region(chrom: str, pos: int, regions: dict) -> bool:
    """
    Check if a genomic position (chrom, pos) falls within any interval
    on that chromosome.

    Regions are pre-grouped by chromosome so only intervals on the same
    chromosome are checked — avoids comparing across chromosomes entirely.

    Returns True if the position overlaps any region, False otherwise.
    """
    if chrom not in regions:
        return False
    for start, end in regions[chrom]:
        if start <= pos <= end:
            return True
    return False


# ---------------------------------------------------------------------------
# Parallelized Filtering
# ---------------------------------------------------------------------------

def filter_chunk(args_tuple) -> tuple:
    """
    Worker function — filters one chunk of candidate rows against the regions.

    Each worker:
        - Reads its assigned chunk file
        - Checks each site against the regions dictionary
        - Writes matching rows to a temp output file

    Returns (chunk_id, match_count, total_count).
    """
    chunk_file, chunk_id, regions = args_tuple
    out_file = os.path.join(TEMP_DIR, f"chunk_{chunk_id}_filtered.tsv")

    matches = 0
    total   = 0
    buffer  = []

    with open(chunk_file, 'r') as f, open(out_file, 'w') as out:
        for line in f:
            line = line.strip()
            if not line:
                continue

            total += 1
            parts = line.split('\t')
            site  = parts[0]

            # Parse site ID: chr_pos format e.g. NC_056679.1_11381
            # Split on last underscore to separate chr from pos
            last_underscore = site.rfind('_')
            if last_underscore == -1:
                continue

            chrom = site[:last_underscore]
            try:
                pos = int(site[last_underscore + 1:])
            except ValueError:
                continue

            if overlaps_any_region(chrom, pos, regions):
                buffer.append(line)
                matches += 1

                if len(buffer) >= CHUNK_SIZE:
                    out.write("\n".join(buffer) + "\n")
                    buffer = []

        if buffer:
            out.write("\n".join(buffer) + "\n")

    return chunk_id, matches, total


def split_candidates_into_chunks(header: str) -> list:
    """
    Split the candidates file into N chunks for parallel processing.
    Each chunk is written to a temp file (without the header).
    Returns list of chunk file paths.
    """
    os.makedirs(TEMP_DIR, exist_ok=True)

    # Count total rows first to size chunks evenly
    with open(CANDIDATES_IN, 'r') as f:
        next(f)  # skip header
        total_rows = sum(1 for _ in f)

    rows_per_chunk = max(1, math.ceil(total_rows / CORES))
    chunk_files    = []
    chunk_id       = 0
    buffer         = []

    with open(CANDIDATES_IN, 'r') as f:
        next(f)  # skip header

        for line in f:
            buffer.append(line)

            if len(buffer) >= rows_per_chunk:
                chunk_file = os.path.join(TEMP_DIR, f"chunk_{chunk_id}.tsv")
                with open(chunk_file, 'w') as out:
                    out.writelines(buffer)
                chunk_files.append(chunk_file)
                chunk_id += 1
                buffer = []

        # Write remaining rows
        if buffer:
            chunk_file = os.path.join(TEMP_DIR, f"chunk_{chunk_id}.tsv")
            with open(chunk_file, 'w') as out:
                out.writelines(buffer)
            chunk_files.append(chunk_file)

    return chunk_files


def merge_chunk_outputs(chunk_files: list, header: str):
    """
    Merge all filtered chunk output files into the final regional output file.
    Writes the header first, then appends each chunk's filtered rows.
    Cleans up chunk files after merging.
    """
    with open(REGIONAL_OUT, 'w') as out:
        out.write(header + "\n")

        for i, chunk_file in enumerate(chunk_files):
            filtered_file = os.path.join(TEMP_DIR, f"chunk_{i}_filtered.tsv")
            if os.path.exists(filtered_file):
                with open(filtered_file, 'r') as f:
                    for line in f:
                        if line.strip():
                            out.write(line)

    # Cleanup
    for i, chunk_file in enumerate(chunk_files):
        for f in [chunk_file,
                  os.path.join(TEMP_DIR, f"chunk_{i}_filtered.tsv")]:
            if os.path.exists(f):
                os.remove(f)

    if os.path.exists(TEMP_DIR) and not os.listdir(TEMP_DIR):
        os.rmdir(TEMP_DIR)


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def main():
    os.makedirs(LOG_DIR, exist_ok=True)

    print("=============================================")
    print(" Step 4 — Filter Lethal Candidates by Region")
    print("=============================================")
    print(f"Threshold        : {THRESHOLD}")
    print(f"Candidates input : {CANDIDATES_IN}")
    print(f"Regions file     : {REGIONS_FILE}")
    print(f"Feature type     : {FEATURE} (ignored for BED files)")
    print(f"Cores            : {CORES}")
    print(f"Output file      : {REGIONAL_OUT}")
    print(f"Timestamp        : {datetime.datetime.now()}")
    print("")

    validate_paths()

    # Load regions into memory — built once, shared across all workers
    log("Loading regions file...")
    regions = load_regions()
    print("")

    # Read header from candidates file
    with open(CANDIDATES_IN, 'r') as f:
        header = f.readline().strip()

    # Split candidates into chunks
    log(f"Splitting candidates into {CORES} chunk(s)...")
    chunk_files = split_candidates_into_chunks(header)
    log(f"  -> {len(chunk_files)} chunk(s) created.")
    print("")

    # Parallel filtering
    log(f"Filtering chunks using {CORES} core(s)...")
    worker_args = [
        (chunk_file, i, regions)
        for i, chunk_file in enumerate(chunk_files)
    ]

    with Pool(processes=CORES) as pool:
        results = pool.map(filter_chunk, worker_args)

    total_evaluated = sum(r[2] for r in results)
    total_matched   = sum(r[1] for r in results)

    for chunk_id, matches, total in sorted(results):
        log(f"  Chunk {chunk_id}: {matches:,} / {total:,} candidates matched")
    print("")

    # Merge outputs
    log("Merging filtered chunks into final output...")
    merge_chunk_outputs(chunk_files, header)

    print("")
    print("=============================================")
    print(" Step 4 Complete")
    print("=============================================")
    print(f"Candidates evaluated : {total_evaluated:,}")
    print(f"Regional candidates  : {total_matched:,} "
          f"({100*total_matched/total_evaluated:.2f}% of lethal candidates)")
    print(f"Output saved to      : {REGIONAL_OUT}")
    print(f"Timestamp            : {datetime.datetime.now()}")
    print("=============================================")


if __name__ == "__main__":
    main()

