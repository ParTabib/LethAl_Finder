#!/usr/bin/env python3
# =============================================================================
# Step 2 — Build Genotype Matrix
# =============================================================================
# Reads the HET site list produced by Step 1 into a Python set for fast O(1)
# lookups, then streams through all HEARTY sample files extracting base calls
# only at HET sites. Writes one temporary file per sample, sorts them, then
# merges all into the final genotype matrix.
#
# Parallelization:
#   Phase 2 — N sample files processed simultaneously (each worker gets its
#              own copy of the HET set)
#   Phase 3 — N temp files sorted simultaneously
#   Phase 4 — Sequential merge (inherently single-threaded)
#
# Memory usage scales with cores × HET set size:
#   0.25 threshold: ~1GB per core  (e.g. 4 cores = ~4GB)
#   0.05 threshold: ~30GB per core (e.g. 4 cores = ~120GB)
#
# Filtering rules:
#   - Missing data is written as NA
#   - UNKNOWN base calls are normalized to NA
#   - Sites where fewer than --min-coverage percent of samples have data
#     are filtered out (default: 15%)
#
# This script is called by the LethAl Finder master script and should not
# be run directly unless you know what you are doing.
#
# Usage:
#   python3 02_build_matrix.py <threshold> <input_dir> <output_dir>
#                              [--min-coverage <percent>]
#                              [--cores <n>]
#                              [--sample-sheet <file>]
# =============================================================================

import sys
import os
import gzip
import datetime
import subprocess
import argparse
import math
from multiprocessing import Pool


# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Build genotype matrix from HEARTY sample files."
)
parser.add_argument("threshold",   type=str, help="Frequency threshold (e.g. 0.25)")
parser.add_argument("input_dir",   type=str, help="Directory with .basecall.txt.gz files")
parser.add_argument("output_dir",  type=str, help="Directory to write output files into")
parser.add_argument("--min-coverage", type=float, default=15.0,
                    help="Min %% of samples with non-NA data per site (default: 15)")
parser.add_argument("--cores",     type=int, default=1,
                    help="Number of CPU cores for parallelization (default: 1)")
parser.add_argument("--sample-sheet", type=str, default=None,
                    help="Optional TSV: filename TAB sample_id")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

THRESHOLD    = str(float(args.threshold))
INPUT_DIR    = args.input_dir
OUTPUT_DIR   = args.output_dir
MIN_COVERAGE = args.min_coverage
CORES        = args.cores
SAMPLE_SHEET = args.sample_sheet

HET_SITES_FILE = os.path.join(OUTPUT_DIR, f"het_sites_{THRESHOLD}.txt")
MATRIX_OUT     = os.path.join(OUTPUT_DIR, f"master_genotype_matrix_{THRESHOLD}.tsv")
TEMP_DIR       = os.path.join(OUTPUT_DIR, "temp_sample_calls")
LOG_DIR        = os.path.join(OUTPUT_DIR, "logs")
BASE_COL       = f"Base_{THRESHOLD}"
STATUS_COL     = f"Status_{THRESHOLD}"
CHUNK_SIZE     = 1_000_000


# ---------------------------------------------------------------------------
# Validate arguments
# ---------------------------------------------------------------------------

if not (0 < MIN_COVERAGE <= 100):
    print(f"[ERROR] --min-coverage must be between 0 and 100 (got {MIN_COVERAGE})")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

def log(msg: str):
    """Print message with timestamp."""
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")
    sys.stdout.flush()


def validate_paths():
    """Check all required input paths exist before processing begins."""
    if not os.path.isfile(HET_SITES_FILE):
        print(f"[ERROR] HET sites file not found: {HET_SITES_FILE}")
        print(f"        Make sure Step 1 has completed successfully.")
        sys.exit(1)
    if not os.path.isdir(INPUT_DIR):
        print(f"[ERROR] Input directory not found: {INPUT_DIR}")
        sys.exit(1)
    if not os.path.isdir(OUTPUT_DIR):
        print(f"[ERROR] Output directory not found: {OUTPUT_DIR}")
        sys.exit(1)
    if SAMPLE_SHEET and not os.path.isfile(SAMPLE_SHEET):
        print(f"[ERROR] Sample sheet not found: {SAMPLE_SHEET}")
        sys.exit(1)


def get_sample_files() -> list:
    """Get sorted list of all gzipped HEARTY files in the input directory."""
    files = sorted([
        os.path.join(INPUT_DIR, f)
        for f in os.listdir(INPUT_DIR)
        if f.endswith(".basecall.txt.gz")
    ])
    if not files:
        print(f"[ERROR] No .basecall.txt.gz files found in: {INPUT_DIR}")
        sys.exit(1)
    return files


def strip_common_suffix(filenames: list) -> dict:
    """
    Automatically derive sample IDs by stripping the longest common suffix.
    e.g. ['lion01.basecall.txt.gz', 'lion02.basecall.txt.gz']
         -> {'lion01.basecall.txt.gz': 'lion01', ...}
    """
    reversed_names = [name[::-1] for name in filenames]
    common = os.path.commonprefix(reversed_names)[::-1]
    mapping = {}
    for name in filenames:
        sample_id = name[:-len(common)] if common else name
        mapping[name] = sample_id if sample_id else name
    return mapping


def load_sample_sheet(sample_files: list) -> dict:
    """Load sample ID mapping from user-provided TSV file."""
    mapping = {}
    with open(SAMPLE_SHEET, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]

    missing = [os.path.basename(f) for f in sample_files
               if os.path.basename(f) not in mapping]
    if missing:
        print(f"[ERROR] Files missing from sample sheet:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)
    return mapping


def build_id_mapping(sample_files: list) -> dict:
    """Build filename -> sample_id mapping from sheet or auto-detection."""
    filenames = [os.path.basename(f) for f in sample_files]
    if SAMPLE_SHEET:
        return load_sample_sheet(sample_files)
    else:
        return strip_common_suffix(filenames)


def get_sample_id(file_path: str, id_mapping: dict) -> str:
    """Return the sample ID for a given file path."""
    return id_mapping[os.path.basename(file_path)]


def validate_threshold_in_file(file_path: str):
    """Check threshold columns exist in the file header."""
    with gzip.open(file_path, 'rt') as f:
        header = f.readline().strip().split('\t')
    if BASE_COL not in header or STATUS_COL not in header:
        print(f"[ERROR] Columns '{BASE_COL}' or '{STATUS_COL}' not found.")
        print(f"        Available threshold columns:")
        for col in header:
            if col.startswith("Base_") or col.startswith("Status_"):
                print(f"          {col}")
        sys.exit(1)


def get_temp_file(sample_id: str) -> str:
    return os.path.join(TEMP_DIR, f"{sample_id}_calls.txt")


def get_sorted_temp_file(sample_id: str) -> str:
    return os.path.join(TEMP_DIR, f"{sample_id}_calls_sorted.txt")


def compute_min_samples(n_samples: int) -> int:
    """Compute minimum samples required per site from coverage percentage."""
    return max(1, math.ceil(n_samples * MIN_COVERAGE / 100))


def cleanup(sample_ids: list):
    """Remove all temporary files and temp directory."""
    log("Cleaning up temporary files...")
    for sid in sample_ids:
        for f in [get_temp_file(sid), get_sorted_temp_file(sid)]:
            if os.path.exists(f):
                os.remove(f)
    if os.path.exists(TEMP_DIR) and not os.listdir(TEMP_DIR):
        os.rmdir(TEMP_DIR)


# ---------------------------------------------------------------------------
# Phase 1 — Load HET Sites into Python Set
# ---------------------------------------------------------------------------

def load_het_sites() -> set:
    """
    Load all HET site IDs from Step 1 output into a Python set.
    Provides O(1) lookup speed during sample file scanning.

    Memory usage scales with the number of HET sites.
    When parallelizing Phase 2, each worker gets its own copy of this set.
    Total RAM = cores × set size.
    """
    log("Phase 1 — Loading HET sites into memory...")
    print(f"  Input: {HET_SITES_FILE}")
    print("")

    het_sites = set()
    total     = 0

    with open(HET_SITES_FILE, 'r') as f:
        for line in f:
            site = line.strip()
            if site:
                het_sites.add(site)
                total += 1
                if total % 10_000_000 == 0:
                    log(f"  Loaded {total:,} sites so far...")

    log(f"  -> {total:,} unique HET sites loaded into memory.")
    print("")
    return het_sites


# ---------------------------------------------------------------------------
# Phase 2 — Write Temporary Files Per Sample (parallelized)
# ---------------------------------------------------------------------------

def process_single_sample(args_tuple) -> tuple:
    """
    Worker function for Phase 2 parallelization.
    Each worker receives its own copy of the HET sites set.
    Streams through one sample file and writes matching sites to a temp file.

    Memory per worker: size of HET set (~1GB at 0.25, ~30GB at 0.05)
    Returns (sample_id, match_count).
    """
    file_path, sample_id, het_sites, base_col, temp_dir, chunk_size = args_tuple

    tmp_file = os.path.join(temp_dir, f"{sample_id}_calls.txt")

    with gzip.open(file_path, 'rt') as f:
        header = f.readline().strip().split('\t')
    chr_idx  = header.index('chr')
    pos_idx  = header.index('pos')
    base_idx = header.index(base_col)

    matches = 0
    buffer  = []

    with gzip.open(file_path, 'rt') as f:
        next(f)  # skip header
        with open(tmp_file, 'w') as out:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) <= base_idx:
                    continue
                site      = parts[chr_idx] + "_" + parts[pos_idx]
                base_call = parts[base_idx]

                if site in het_sites:
                    # Normalize UNKNOWN base calls to NA
                    base_call = "NA" if base_call.strip().upper() == "UNKNOWN" \
                                else base_call
                    buffer.append(f"{site}\t{base_call}\n")
                    matches += 1

                    if len(buffer) >= chunk_size:
                        out.writelines(buffer)
                        buffer = []

            if buffer:
                out.writelines(buffer)

    return sample_id, matches


def write_sample_temp_files(sample_files: list, het_sites: set,
                             id_mapping: dict):
    """
    Process all sample files in parallel using --cores workers.
    Each worker gets its own copy of the HET set for O(1) lookups.
    """
    log(f"Phase 2 — Processing {len(sample_files)} sample files "
        f"using {CORES} core(s)...")
    print("")

    os.makedirs(TEMP_DIR, exist_ok=True)

    worker_args = [
        (f, get_sample_id(f, id_mapping), het_sites,
         BASE_COL, TEMP_DIR, CHUNK_SIZE)
        for f in sample_files
    ]

    with Pool(processes=CORES) as pool:
        results = pool.map(process_single_sample, worker_args)

    for sample_id, matches in sorted(results):
        log(f"  {sample_id}: {matches:,} HET sites matched")

    print("")


# ---------------------------------------------------------------------------
# Phase 3 — Sort Temporary Files (parallelized)
# ---------------------------------------------------------------------------

def sort_single_file(args_tuple):
    """Worker function — sorts one temp file and removes the unsorted version."""
    sample_id, temp_dir = args_tuple
    tmp_file    = os.path.join(temp_dir, f"{sample_id}_calls.txt")
    sorted_file = os.path.join(temp_dir, f"{sample_id}_calls_sorted.txt")
    subprocess.run(["sort", "-k1,1", tmp_file, "-o", sorted_file], check=True)
    os.remove(tmp_file)
    return sample_id


def sort_temp_files(sample_ids: list):
    """
    Sort all sample temp files simultaneously using --cores workers.
    Reduces Phase 3 from ~N × sort_time to ~sort_time of one file.
    """
    log(f"Phase 3 — Sorting {len(sample_ids)} temp files "
        f"using {CORES} core(s)...")
    print("")

    sort_args = [(sid, TEMP_DIR) for sid in sample_ids]

    with Pool(processes=CORES) as pool:
        sorted_ids = pool.map(sort_single_file, sort_args)

    for sid in sorted(sorted_ids):
        log(f"  Sorted: {sid}")

    log("  -> All temporary files sorted.")
    print("")


# ---------------------------------------------------------------------------
# Phase 4 — Merge Into Final Matrix (sequential)
# ---------------------------------------------------------------------------

def sorted_file_generator(file_path: str):
    """Generator that streams a sorted temp file line by line."""
    if not os.path.exists(file_path):
        return
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                yield parts[0], parts[1]


def merge_into_matrix(sample_ids: list, min_samples: int):
    """
    Merge all sorted temporary sample files into the final genotype matrix.

    Opens all sorted temp files as generators simultaneously.
    For each site in the HET sites file, advances each sample's pointer
    to find a match — memory usage is ~N lines total at any point.

    Writes the final matrix to disk row by row in chunks.
    This phase is inherently sequential — the merge walks all files
    in lockstep and cannot be parallelized.
    """
    log("Phase 4 — Merging into final genotype matrix...")
    print(f"  Output                            : {MATRIX_OUT}")
    print(f"  Min coverage                      : {MIN_COVERAGE}%")
    print(f"  Min samples required per site     : {min_samples} / {len(sample_ids)}")
    print("")

    header_line = "Site\t" + "\t".join(sample_ids) + "\n"

    generators = {
        sid: sorted_file_generator(get_sorted_temp_file(sid))
        for sid in sample_ids
    }

    current = {}
    for sid, gen in generators.items():
        try:
            current[sid] = next(gen)
        except StopIteration:
            current[sid] = None

    written  = 0
    filtered = 0
    buffer   = []

    with open(MATRIX_OUT, 'w') as out:
        out.write(header_line)

        with open(HET_SITES_FILE, 'r') as sites_f:
            for line in sites_f:
                site = line.strip()
                if not site:
                    continue

                calls  = []
                non_na = 0

                for sid in sample_ids:
                    while (current[sid] is not None and
                           current[sid][0] < site):
                        try:
                            current[sid] = next(generators[sid])
                        except StopIteration:
                            current[sid] = None

                    if current[sid] is not None and current[sid][0] == site:
                        calls.append(current[sid][1])
                        non_na += 1
                    else:
                        calls.append("NA")

                if non_na < min_samples:
                    filtered += 1
                    continue

                buffer.append(site + "\t" + "\t".join(calls))
                written += 1

                if len(buffer) >= CHUNK_SIZE:
                    out.write("\n".join(buffer) + "\n")
                    buffer = []

                if written % 10_000_000 == 0:
                    log(f"  Written {written:,} rows so far "
                        f"(filtered {filtered:,} low-coverage sites)")

        if buffer:
            out.write("\n".join(buffer) + "\n")

    log(f"  -> Matrix saved to                : {MATRIX_OUT}")
    log(f"  -> Rows written                   : {written:,}")
    log(f"  -> Rows filtered (< {MIN_COVERAGE}% coverage) : {filtered:,}")
    print("")


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def main():
    os.makedirs(LOG_DIR, exist_ok=True)

    sample_files = get_sample_files()
    id_mapping   = build_id_mapping(sample_files)
    sample_ids   = [get_sample_id(f, id_mapping) for f in sample_files]
    min_samples  = compute_min_samples(len(sample_files))

    # Memory estimation
    set_size_gb  = round(sum(1 for _ in open(HET_SITES_FILE)) * 75 / 1e9, 2) \
                   if os.path.isfile(HET_SITES_FILE) else "?"
    total_ram_gb = round(set_size_gb * CORES, 2) \
                   if isinstance(set_size_gb, float) else "?"

    print("=============================================")
    print(" Step 2 — Build Genotype Matrix")
    print("=============================================")
    print(f"Threshold        : {THRESHOLD}")
    print(f"Min coverage     : {MIN_COVERAGE}% "
          f"(= {min_samples} samples out of {len(sample_files)})")
    print(f"Cores            : {CORES}")
    print(f"Est. RAM (Phase 2): ~{total_ram_gb} GB "
          f"({set_size_gb} GB × {CORES} core(s))")
    print(f"HET sites file   : {HET_SITES_FILE}")
    print(f"Input directory  : {INPUT_DIR}")
    print(f"Output file      : {MATRIX_OUT}")
    print(f"Timestamp        : {datetime.datetime.now()}")
    print("")

    validate_paths()
    print(f"Sample files detected: {len(sample_files)}")
    print("")

    log("Sample ID mapping:")
    for file_path in sample_files:
        fname     = os.path.basename(file_path)
        sample_id = get_sample_id(file_path, id_mapping)
        print(f"  {fname} -> {sample_id}")
    print("")

    validate_threshold_in_file(sample_files[0])

    # Phase 1 — load HET sites
    het_sites = load_het_sites()

    # Phase 2 — parallel sample processing
    write_sample_temp_files(sample_files, het_sites, id_mapping)

    # Phase 3 — parallel sorting
    sort_temp_files(sample_ids)

    # Phase 4 — sequential merge
    merge_into_matrix(sample_ids, min_samples)

    # Cleanup
    cleanup(sample_ids)

    print("=============================================")
    print(" Step 2 Complete")
    print("=============================================")
    print(f"Output saved to  : {MATRIX_OUT}")
    print(f"Timestamp        : {datetime.datetime.now()}")
    print("=============================================")


if __name__ == "__main__":
    main()
