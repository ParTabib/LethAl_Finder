#!/usr/bin/env python3
# =============================================================================
# Step 5 — Annotate Lethal Candidates with Functional Consequences
# =============================================================================
# Takes the regional lethal candidates from Step 4 and annotates each site
# with its predicted functional consequence — implemented entirely in pure
# Python without any external tools or dependencies.
#
# For each candidate site the script:
#   1. Determines the reference allele by reading directly from the FASTA
#      file using the .fai index for fast byte-offset seeking
#   2. Identifies which CDS feature(s) from the GFF contain the position
#   3. Extracts the codon containing the variant from the reference FASTA
#   4. Substitutes the alternate allele into the codon
#   5. Translates both reference and alternate codons using the standard
#      genetic code
#   6. Classifies the consequence: synonymous, missense, stop_gained,
#      stop_lost, or start_lost
#
# Parallelized — the regional candidates file is split into N chunks and
# each chunk is annotated simultaneously. The CDS dictionary and FASTA
# index are built once before the pool starts and shared read-only.
#
# Only GFF/GFF3 format is supported — CDS features are required to
# determine reading frame (phase). BED and GTF files are not supported
# and the step exits with a clear message if a non-GFF file is provided.
#
# No external tools required — pure Python standard library only.
#
# This script is called by the LethAl Finder master script and should not
# be run directly unless you know what you are doing.
#
# Usage:
#   python3 05_annotate_variants.py <threshold> <output_dir> <regions_file>
#                                   --fasta <ref.fna>
#                                   [--cores <n>]
#   e.g. python3 05_annotate_variants.py 0.25 /results genome.gff
#            --fasta ref.fna --cores 4
# =============================================================================

import sys
import os
import gzip
import math
import datetime
import argparse
from multiprocessing import Pool


# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Annotate lethal candidates with functional consequences."
)
parser.add_argument("threshold",    type=str, help="Frequency threshold (e.g. 0.25)")
parser.add_argument("output_dir",   type=str, help="Directory containing regional candidates file")
parser.add_argument("regions_file", type=str, help="GFF3 annotation file used in Step 4")
parser.add_argument("--fasta",      type=str, required=True,
                    help="Reference FASTA file (.fna) with .fai index")
parser.add_argument("--cores",      type=int, default=4,
                    help="Number of CPU cores for parallelization (default: 4)")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

THRESHOLD     = str(float(args.threshold))
OUTPUT_DIR    = args.output_dir
REGIONS_FILE  = args.regions_file
FASTA         = args.fasta
CORES         = args.cores

LOG_DIR       = os.path.join(OUTPUT_DIR, "logs")
REGIONAL_IN   = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}_regional.tsv")
ANNOTATED_OUT = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}_annotated.tsv")
TEMP_DIR      = os.path.join(OUTPUT_DIR, "temp_annotation_chunks")
CHUNK_SIZE    = 100_000  # rows per write chunk


# ---------------------------------------------------------------------------
# Standard Genetic Code
# ---------------------------------------------------------------------------
# Maps every 64 three-base codons to a single-letter amino acid code.
# * denotes a stop codon.
# Used by translate() to convert reference and alternate codons.

CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Complement mapping for reverse strand handling
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C",
              "a": "t", "t": "a", "c": "g", "g": "c"}


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

def log(msg: str):
    """Print message with timestamp."""
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")
    sys.stdout.flush()


def detect_format(file_path: str) -> str:
    """
    Detect file format from extension. Strips .gz before checking.
    Returns one of: 'bed', 'gff', 'gtf', 'unknown'.
    """
    ext = file_path.lower()
    if ext.endswith(".gz"):
        ext = ext[:-3]
    if ext.endswith(".bed"):
        return "bed"
    elif ext.endswith(".gff") or ext.endswith(".gff3"):
        return "gff"
    elif ext.endswith(".gtf"):
        return "gtf"
    return "unknown"


def validate_inputs():
    """Check all required input paths and files exist before processing."""
    errors = []
    if not os.path.isfile(REGIONAL_IN):
        errors.append(f"Regional candidates file not found: {REGIONAL_IN}\n"
                      f"        Make sure Step 4 has completed successfully.")
    if not os.path.isfile(FASTA):
        errors.append(f"Reference FASTA not found: {FASTA}")
    if not os.path.isfile(FASTA + ".fai"):
        errors.append(f"FASTA index not found: {FASTA}.fai\n"
                      f"        Run: samtools faidx {FASTA}")
    if not os.path.isfile(REGIONS_FILE):
        errors.append(f"Regions file not found: {REGIONS_FILE}")
    if errors:
        print("\n[ERROR] Validation failed:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    Used for genes on the minus strand — the coding sequence runs
    in the opposite direction so the codon must be reverse complemented
    before translation.
    """
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))


def translate(codon: str) -> str:
    """
    Translate a three-base codon to a single-letter amino acid.
    Returns '?' for unknown or degenerate codons.
    """
    return CODON_TABLE.get(codon.upper(), "?")


# ---------------------------------------------------------------------------
# FASTA Index Reader
# ---------------------------------------------------------------------------

def load_fasta_index(fai_path: str) -> dict:
    """
    Read the .fai FASTA index file and return a dictionary mapping
    chromosome names to their index information.

    The .fai format has five tab-separated columns:
        col 0: chromosome name
        col 1: number of bases in the sequence
        col 2: byte offset of the first base in the FASTA file
        col 3: number of bases per line
        col 4: number of bytes per line (includes the newline character)

    This information allows us to calculate the exact byte position of
    any base in the FASTA file without reading the whole file into memory.

    Returns:
        {chrom: (length, offset, bases_per_line, bytes_per_line)}
    """
    index = {}
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom           = parts[0]
            length          = int(parts[1])
            offset          = int(parts[2])
            bases_per_line  = int(parts[3])
            bytes_per_line  = int(parts[4])
            index[chrom] = (length, offset, bases_per_line, bytes_per_line)
    return index


def get_sequence(fasta_path: str, fasta_index: dict,
                 chrom: str, start: int, end: int) -> str:
    """
    Extract a subsequence from the reference FASTA using the .fai index
    for direct byte-offset seeking — no external tools required.

    Coordinates are 1-based closed intervals (same as GFF).

    The .fai index tells us:
        - Where the chromosome starts in the file (byte offset)
        - How many bases per line and bytes per line

    From these we calculate the exact byte position of any base:
        line_number  = (pos - 1) // bases_per_line
        col_in_line  = (pos - 1) %  bases_per_line
        byte_offset  = chrom_offset + line_number * bytes_per_line + col_in_line

    We then seek to that position and read the required number of bytes,
    stripping newline characters.

    Returns the sequence as an uppercase string, or empty string on error.
    """
    if chrom not in fasta_index:
        return ""

    length, chrom_offset, bases_per_line, bytes_per_line = fasta_index[chrom]

    # Clamp to chromosome boundaries
    start = max(1, start)
    end   = min(length, end)
    if start > end:
        return ""

    seq = []
    try:
        with open(fasta_path, 'rb') as f:
            pos = start
            while pos <= end:
                # Calculate byte offset for this position
                line_num    = (pos - 1) // bases_per_line
                col_in_line = (pos - 1) %  bases_per_line
                byte_pos    = chrom_offset + line_num * bytes_per_line + col_in_line

                # How many bases remain on this line
                bases_on_line = bases_per_line - col_in_line
                # How many bases we still need
                bases_needed  = end - pos + 1
                # Read the minimum of the two
                to_read = min(bases_on_line, bases_needed)

                f.seek(byte_pos)
                raw = f.read(to_read).decode('ascii').replace('\n', '').replace('\r', '')
                seq.append(raw)
                pos += len(raw)
    except Exception:
        return ""

    return "".join(seq).upper()


# ---------------------------------------------------------------------------
# GFF CDS Parser
# ---------------------------------------------------------------------------

def parse_cds_from_gff(file_path: str) -> dict:
    """
    Parse a GFF3 file and extract all CDS features.

    CDS features are required (not exon) because only CDS lines carry
    the phase column (col 7) — the reading frame offset needed to
    determine which codon position a variant falls in.

    GFF3 CDS format:
        col 0: chromosome
        col 2: feature type (must be 'CDS')
        col 3: start (1-based)
        col 4: end (inclusive)
        col 6: strand ('+' or '-')
        col 7: phase (0, 1, or 2 — number of bases to skip at the start
                of this CDS to reach the first complete codon)
        col 8: attributes (contains gene= or Parent= for gene name)

    Returns:
        {chrom: [(start, end, strand, phase, gene_name), ...]}
    """
    cds_dict = {}
    skipped  = 0
    matched  = 0

    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                skipped += 1
                continue

            feature_type = parts[2]
            if feature_type != "CDS":
                continue

            chrom  = parts[0]
            strand = parts[6]
            try:
                start = int(parts[3])
                end   = int(parts[4])
                phase = int(parts[7]) if parts[7] != '.' else 0
            except ValueError:
                skipped += 1
                continue

            # Extract gene name from attributes column
            # Attributes look like: ID=cds-XM_001;Parent=rna-XM_001;gene=LOC122228025
            gene_name = "."
            attrs = parts[8]
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith("gene="):
                    gene_name = attr[5:]
                    break
                if attr.startswith("Parent=") and gene_name == ".":
                    gene_name = attr[7:]

            if chrom not in cds_dict:
                cds_dict[chrom] = []
            cds_dict[chrom].append((start, end, strand, phase, gene_name))
            matched += 1

    if matched == 0:
        print(f"[ERROR] No CDS features found in: {file_path}")
        sys.exit(1)

    if skipped > 0:
        log(f"  [WARNING] Skipped {skipped} malformed GFF lines.")

    return cds_dict


# ---------------------------------------------------------------------------
# Variant Consequence Prediction
# ---------------------------------------------------------------------------

def find_cds_for_position(chrom: str, pos: int, cds_dict: dict) -> list:
    """
    Find all CDS features that contain the given position.
    Returns a list of (start, end, strand, phase, gene_name) tuples.
    Multiple CDS features can overlap a position when transcripts overlap.
    """
    if chrom not in cds_dict:
        return []
    return [(s, e, strand, phase, gene)
            for s, e, strand, phase, gene in cds_dict[chrom]
            if s <= pos <= e]


def get_codon_for_position(fasta_path: str, fasta_index: dict,
                           chrom: str, pos: int,
                           cds_start: int, cds_end: int,
                           strand: str, phase: int) -> tuple:
    """
    Extract the codon containing the given position within a CDS feature.

    For a plus-strand gene:
        - The first complete codon base is cds_start + phase
        - codon_index = (pos - first_codon_base) // 3
        - codon starts at: first_codon_base + codon_index * 3

    For a minus-strand gene:
        - The CDS runs right to left from cds_end
        - The first complete codon base is cds_end - phase
        - The codon must be reverse complemented before translation

    Returns (ref_codon, codon_pos_in_codon, codon_number) or
            (None, None, None) if position is outside coding region.
    """
    if strand == "+":
        first_codon_base = cds_start + phase
        if pos < first_codon_base:
            return None, None, None

        offset             = pos - first_codon_base
        codon_index        = offset // 3
        codon_pos_in_codon = offset %  3
        codon_number       = codon_index + 1
        codon_start        = first_codon_base + codon_index * 3
        codon_end          = codon_start + 2

        ref_codon = get_sequence(fasta_path, fasta_index, chrom, codon_start, codon_end)
        return ref_codon, codon_pos_in_codon, codon_number

    else:  # minus strand
        first_codon_base = cds_end - phase
        if pos > first_codon_base:
            return None, None, None

        offset              = first_codon_base - pos
        codon_index         = offset // 3
        codon_pos_in_codon  = offset %  3
        codon_number        = codon_index + 1
        codon_end_genomic   = first_codon_base - codon_index * 3
        codon_start_genomic = codon_end_genomic - 2

        ref_seq   = get_sequence(fasta_path, fasta_index, chrom,
                                 codon_start_genomic, codon_end_genomic)
        ref_codon = reverse_complement(ref_seq)
        # Mirror the codon position for minus strand
        codon_pos_in_codon = 2 - codon_pos_in_codon
        return ref_codon, codon_pos_in_codon, codon_number


def predict_consequence(ref_codon: str, alt_allele: str,
                        codon_pos: int, codon_number: int,
                        strand: str) -> tuple:
    """
    Predict the functional consequence of substituting alt_allele into
    the reference codon at the given position.

    For minus-strand genes the alt allele is complemented before
    substitution because we work with the coding strand sequence.

    Consequences:
        synonymous  — same amino acid
        missense    — different amino acid
        stop_gained — alternate codon is a stop codon
        stop_lost   — reference codon was a stop, alternate is not
        start_lost  — first codon (Met) changed to something else

    Returns (effect, aa_ref, aa_alt, aa_change_string).
    """
    if len(ref_codon) != 3:
        return "unknown", ".", ".", "."

    # For minus strand, complement the alternate allele
    if strand == "-":
        alt_allele = COMPLEMENT.get(alt_allele.upper(), alt_allele)

    # Build alternate codon by substituting at the correct position
    alt_list        = list(ref_codon.upper())
    alt_list[codon_pos] = alt_allele.upper()
    alt_codon_str   = "".join(alt_list)

    aa_ref = translate(ref_codon.upper())
    aa_alt = translate(alt_codon_str)

    # Classify consequence
    if aa_ref == aa_alt:
        effect = "synonymous"
    elif aa_alt == "*":
        effect = "stop_gained"
    elif aa_ref == "*":
        effect = "stop_lost"
    elif codon_number == 1 and aa_ref == "M" and aa_alt != "M":
        effect = "start_lost"
    else:
        effect = "missense"

    aa_change = f"{aa_ref}{codon_number}{aa_alt}"
    return effect, aa_ref, aa_alt, aa_change


def annotate_site(chrom: str, pos: int, ref_allele: str, alt_allele: str,
                  cds_dict: dict, fasta_path: str, fasta_index: dict) -> tuple:
    """
    Annotate a single candidate site with its functional consequence.

    Finds all CDS features overlapping the position, predicts the
    consequence for each, and returns the most severe one.

    Severity ranking (highest to lowest):
        stop_gained > start_lost > stop_lost > missense > synonymous > unknown

    Returns (effect, gene_name, aa_change).
    """
    SEVERITY = {
        "stop_gained": 6,
        "start_lost":  5,
        "stop_lost":   4,
        "missense":    3,
        "synonymous":  2,
        "unknown":     1,
    }

    cds_hits = find_cds_for_position(chrom, pos, cds_dict)
    if not cds_hits:
        return "intergenic", ".", "."

    best_effect    = "unknown"
    best_gene      = "."
    best_aa_change = "."
    best_severity  = 0

    for cds_start, cds_end, strand, phase, gene_name in cds_hits:
        ref_codon, codon_pos, codon_number = get_codon_for_position(
            fasta_path, fasta_index, chrom, pos,
            cds_start, cds_end, strand, phase
        )

        if ref_codon is None:
            continue

        effect, aa_ref, aa_alt, aa_change = predict_consequence(
            ref_codon, alt_allele, codon_pos, codon_number, strand
        )

        severity = SEVERITY.get(effect, 0)
        if severity > best_severity:
            best_severity  = severity
            best_effect    = effect
            best_gene      = gene_name
            best_aa_change = aa_change

    if best_severity == 0:
        return "intergenic", ".", "."

    return best_effect, best_gene, best_aa_change


# ---------------------------------------------------------------------------
# Reference Allele Lookup
# ---------------------------------------------------------------------------

def get_ref_allele(chrom: str, pos: int,
                   fasta_path: str, fasta_index: dict) -> str:
    """
    Get the reference base at a single position using the FASTA index.
    Returns the base as an uppercase string, or None if lookup fails.
    """
    seq = get_sequence(fasta_path, fasta_index, chrom, pos, pos)
    if seq and seq.upper() in ("A", "C", "G", "T"):
        return seq.upper()
    return None


def get_alleles_from_genotypes(genotypes: list) -> set:
    """
    Extract all unique alleles seen across all genotypes at a site.
    Returns a set of single-character allele strings e.g. {'A', 'C'}.
    """
    alleles = set()
    for g in genotypes:
        if g not in ("NA", "UNKNOWN", "") and g is not None:
            g = g.replace("/", "").replace("|", "").strip()
            if len(g) == 2:
                alleles.add(g[0].upper())
                alleles.add(g[1].upper())
    return alleles


# ---------------------------------------------------------------------------
# Parallelized Annotation
# ---------------------------------------------------------------------------

def annotate_chunk(args_tuple) -> tuple:
    """
    Worker function — annotates one chunk of regional candidate rows.

    Each worker:
        - Reads its assigned chunk file row by row
        - For each site, looks up the reference allele from the FASTA
        - Determines the alternate allele from the genotypes
        - Predicts the functional consequence
        - Writes the annotated row to a temp output file

    Returns (chunk_id, annotated_count, total_count).
    """
    chunk_file, chunk_id, cds_dict, fasta_path, fasta_index = args_tuple
    out_file = os.path.join(TEMP_DIR, f"chunk_{chunk_id}_annotated.tsv")

    annotated = 0
    total     = 0
    buffer    = []

    with open(chunk_file, 'r') as f, open(out_file, 'w') as out:
        for line in f:
            line = line.strip()
            if not line:
                continue

            total += 1
            parts     = line.split('\t')
            site      = parts[0]
            genotypes = parts[1:]

            # Parse site ID: chr_pos format e.g. NC_056679.1_100064397
            last_underscore = site.rfind('_')
            if last_underscore == -1:
                buffer.append(f"{line}\tunknown\t.\t.")
                continue

            chrom = site[:last_underscore]
            try:
                pos = int(site[last_underscore + 1:])
            except ValueError:
                buffer.append(f"{line}\tunknown\t.\t.")
                continue

            # Get reference allele directly from FASTA
            ref = get_ref_allele(chrom, pos, fasta_path, fasta_index)
            if ref is None:
                buffer.append(f"{line}\tref_lookup_failed\t.\t.")
                continue

            # Get alternate allele from genotypes
            alleles = get_alleles_from_genotypes(genotypes)
            alts    = alleles - {ref}
            if not alts:
                buffer.append(f"{line}\tno_alt_allele\t.\t.")
                continue
            alt = sorted(alts)[0]

            # Predict consequence
            effect, gene, aa_change = annotate_site(
                chrom, pos, ref, alt, cds_dict, fasta_path, fasta_index
            )

            buffer.append(f"{line}\t{effect}\t{gene}\t{aa_change}")
            annotated += 1

            if len(buffer) >= CHUNK_SIZE:
                out.write("\n".join(buffer) + "\n")
                buffer = []

        if buffer:
            out.write("\n".join(buffer) + "\n")

    return chunk_id, annotated, total


def split_candidates_into_chunks(header: str) -> list:
    """
    Split the regional candidates file into N chunks for parallel processing.
    Each chunk is written to a temp file without the header.
    Returns list of chunk file paths.
    """
    os.makedirs(TEMP_DIR, exist_ok=True)

    with open(REGIONAL_IN, 'r') as f:
        next(f)  # skip header
        total_rows = sum(1 for _ in f)

    rows_per_chunk = max(1, math.ceil(total_rows / CORES))
    chunk_files    = []
    chunk_id       = 0
    buffer         = []

    with open(REGIONAL_IN, 'r') as f:
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
        if buffer:
            chunk_file = os.path.join(TEMP_DIR, f"chunk_{chunk_id}.tsv")
            with open(chunk_file, 'w') as out:
                out.writelines(buffer)
            chunk_files.append(chunk_file)

    return chunk_files


def merge_chunk_outputs(chunk_files: list, header: str) -> dict:
    """
    Merge all annotated chunk output files into the final annotated file.
    Writes the header with new annotation columns first, then each chunk.
    Cleans up chunk files after merging.
    Returns effect counts for the summary.
    """
    effect_counts = {}

    with open(ANNOTATED_OUT, 'w') as out:
        out.write(header + "\tEffect\tGene\tAA_change\n")

        for i, chunk_file in enumerate(chunk_files):
            ann_file = os.path.join(TEMP_DIR, f"chunk_{i}_annotated.tsv")
            if os.path.exists(ann_file):
                with open(ann_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        out.write(line + "\n")
                        # Count effects from the Effect column (third from last)
                        parts = line.split('\t')
                        if len(parts) >= 3:
                            effect = parts[-3]
                            effect_counts[effect] = effect_counts.get(effect, 0) + 1

    # Cleanup temp files
    for i, chunk_file in enumerate(chunk_files):
        for f in [chunk_file,
                  os.path.join(TEMP_DIR, f"chunk_{i}_annotated.tsv")]:
            if os.path.exists(f):
                os.remove(f)

    if os.path.exists(TEMP_DIR) and not os.listdir(TEMP_DIR):
        os.rmdir(TEMP_DIR)

    return effect_counts


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def main():
    os.makedirs(LOG_DIR, exist_ok=True)

    print("=============================================")
    print(" Step 5 — Annotate Variants")
    print("=============================================")
    print(f"Threshold        : {THRESHOLD}")
    print(f"Regional input   : {REGIONAL_IN}")
    print(f"Regions file     : {REGIONS_FILE}")
    print(f"Reference FASTA  : {FASTA}")
    print(f"Cores            : {CORES}")
    print(f"Output file      : {ANNOTATED_OUT}")
    print(f"Timestamp        : {datetime.datetime.now()}")
    print("")

    # Check format — only GFF3 supported (CDS phase information required)
    fmt = detect_format(REGIONS_FILE)
    if fmt != "gff":
        print(f"[WARNING] Step 5 requires a GFF3 file (CDS features with phase).")
        print(f"          Provided file appears to be {fmt.upper()} format.")
        print(f"          Skipping annotation step.")
        sys.exit(0)

    validate_inputs()

    # ── Phase 1: Parse GFF to build CDS dictionary ───────────────────────────
    log("Phase 1 — Parsing GFF3 for CDS features...")
    cds_dict  = parse_cds_from_gff(REGIONS_FILE)
    total_cds = sum(len(v) for v in cds_dict.values())
    log(f"  -> {total_cds:,} CDS features loaded across {len(cds_dict):,} chromosomes.")
    print("")

    # ── Phase 2: Load FASTA index ─────────────────────────────────────────────
    log("Phase 2 — Loading FASTA index...")
    fasta_index = load_fasta_index(FASTA + ".fai")
    log(f"  -> {len(fasta_index):,} chromosomes indexed.")
    print("")

    # ── Phase 3: Split candidates into chunks ─────────────────────────────────
    log("Phase 3 — Splitting regional candidates into chunks...")
    with open(REGIONAL_IN, 'r') as f:
        header = f.readline().strip()

    chunk_files = split_candidates_into_chunks(header)
    log(f"  -> {len(chunk_files)} chunk(s) created.")
    print("")

    # ── Phase 4: Parallel annotation ─────────────────────────────────────────
    log(f"Phase 4 — Annotating chunks using {CORES} core(s)...")
    worker_args = [
        (chunk_file, i, cds_dict, FASTA, fasta_index)
        for i, chunk_file in enumerate(chunk_files)
    ]

    with Pool(processes=CORES) as pool:
        results = pool.map(annotate_chunk, worker_args)

    total_evaluated = sum(r[2] for r in results)
    total_annotated = sum(r[1] for r in results)

    for chunk_id, annotated, total in sorted(results):
        log(f"  Chunk {chunk_id}: {annotated:,} / {total:,} sites annotated")
    print("")

    # ── Phase 5: Merge chunk outputs ─────────────────────────────────────────
    log("Phase 5 — Merging annotated chunks into final output...")
    effect_counts = merge_chunk_outputs(chunk_files, header)

    print("")
    print("=============================================")
    print(" Step 5 Complete")
    print("=============================================")
    print(f"Sites evaluated      : {total_evaluated:,}")
    print(f"Sites annotated      : {total_annotated:,}")
    print("")
    print("  Effect breakdown:")
    for effect, count in sorted(effect_counts.items(), key=lambda x: -x[1]):
        print(f"    {effect:<35} {count:,}")
    print("")
    print(f"Output saved to      : {ANNOTATED_OUT}")
    print(f"Timestamp            : {datetime.datetime.now()}")
    print("=============================================")


if __name__ == "__main__":
    main()
