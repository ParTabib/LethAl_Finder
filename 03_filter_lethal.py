#!/usr/bin/env python3
# =============================================================================
# Step 3 — Filter Lethal Candidates
# =============================================================================
# Reads the genotype matrix produced by Step 2 and applies the lethal
# recessive selection logic row by row.
#
# For each site, the logic is:
#   1. Collect all non-NA genotypes across all samples
#   2. Find all HET genotypes (e.g. AC, GT) and extract their allele pairs
#   3. For each HET allele pair, check if BOTH corresponding homozygotes
#      exist somewhere in the population
#   4. If either homozygote is missing -> lethal candidate
#
# Example:
#   Site 11383: sample01=AC, sample02=TT, sample03=AA
#   HET pair: A and C -> need AA and CC
#   AA exists (sample03) v  |  CC never appears x  -> lethal candidate
#   TT from sample02 is irrelevant to the AC pair assessment
#
# Output files:
#   1. lethal_candidates_<threshold>.tsv
#      Full matrix rows for lethal candidate sites (same format as input)
#
#   2. lethal_hits_<threshold>.tsv
#      One row per sample carrying a HET genotype at a lethal site:
#      Site | Sample | Genotype
#
# This script is called by the LethAl Finder master script and should not
# be run directly unless you know what you are doing.
#
# Usage:
#   python3 03_filter_lethal.py <threshold> <output_dir>
#   e.g. python3 03_filter_lethal.py 0.25 /results
# =============================================================================

import sys
import os
import datetime
import argparse


# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Filter lethal recessive candidate sites from genotype matrix."
)
parser.add_argument(
    "threshold",
    type=str,
    help="Frequency threshold used in Steps 1 and 2 (e.g. 0.25)"
)
parser.add_argument(
    "output_dir",
    type=str,
    help="Directory containing the genotype matrix and where output will be written"
)
parser.add_argument(
    "--min-carriers",
    type=int,
    default=1,
    help="Minimum number of individuals carrying the HET genotype at a site "
         "for it to be considered a lethal candidate (default: 1)"
)
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

THRESHOLD      = str(float(args.threshold))
OUTPUT_DIR     = args.output_dir
MIN_CARRIERS   = args.min_carriers
LOG_DIR        = os.path.join(OUTPUT_DIR, "logs")
MATRIX_IN      = os.path.join(OUTPUT_DIR, f"master_genotype_matrix_{THRESHOLD}.tsv")
CANDIDATES_OUT = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}.tsv")
HITS_OUT       = os.path.join(OUTPUT_DIR, f"lethal_hits_{THRESHOLD}.tsv")
CHUNK_SIZE     = 100_000  # rows per write chunk


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

def log(msg: str):
    """Print message with timestamp."""
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")
    sys.stdout.flush()


def validate_paths():
    """Check all required input paths exist before processing begins."""
    if not os.path.isfile(MATRIX_IN):
        print(f"[ERROR] Matrix file not found: {MATRIX_IN}")
        print(f"        Make sure Step 2 has completed successfully.")
        sys.exit(1)
    if not os.path.isdir(OUTPUT_DIR):
        print(f"[ERROR] Output directory not found: {OUTPUT_DIR}")
        sys.exit(1)


def is_valid(genotype: str) -> bool:
    """
    Return True if the genotype is a real call — not NA or UNKNOWN.
    """
    return genotype not in ("NA", "UNKNOWN", "") and genotype is not None


def normalize(genotype: str) -> str:
    """
    Normalize genotype string by stripping separators.
    Handles formats: 'AC', 'A/C', 'A|C' -> 'AC'
    """
    return genotype.replace("/", "").replace("|", "").strip()


def is_homozygous(genotype: str) -> bool:
    """Return True if genotype is homozygous (e.g. AA, TT)."""
    g = normalize(genotype)
    return len(g) == 2 and g[0] == g[1]


def is_heterozygous(genotype: str) -> bool:
    """Return True if genotype is heterozygous (e.g. AC, GT)."""
    g = normalize(genotype)
    return len(g) == 2 and g[0] != g[1]


def is_lethal_candidate(genotypes: list) -> bool:
    """
    Apply the lethal recessive selection logic for a single site.

    For each HET genotype, extract its two alleles (e.g. AC -> A, C)
    and check whether BOTH corresponding homozygotes (AA and CC) exist
    somewhere in the population.

    Returns True if the site is a lethal candidate, False if neutral.

    Logic:
        - Build the set of alleles that appear as homozygotes in the population
        - For each HET genotype, check if both alleles have a homozygous form
        - If either is missing -> lethal candidate
    """
    # Build set of alleles that appear as homozygotes somewhere in population
    observed_hom_alleles = set()
    for g in genotypes:
        if is_valid(g) and is_homozygous(g):
            observed_hom_alleles.add(normalize(g)[0])

    # Check each HET genotype's allele pair
    for g in genotypes:
        if is_valid(g) and is_heterozygous(g):
            g_norm  = normalize(g)
            allele1 = g_norm[0]
            allele2 = g_norm[1]
            # If either allele is missing its homozygous form -> lethal
            if allele1 not in observed_hom_alleles or allele2 not in observed_hom_alleles:
                return True

    return False


# ---------------------------------------------------------------------------
# Core Pipeline
# ---------------------------------------------------------------------------

def filter_lethal_candidates():
    """
    Stream through the genotype matrix row by row.
    Apply the lethal candidate logic to each site.

    For each lethal candidate site:
        - Write the full row to lethal_candidates output
        - Write one entry per sample carrying HET at this site to lethal_hits

    Memory usage: one row at a time regardless of matrix size.
    """
    log("Streaming through genotype matrix...")
    print(f"  Input    : {MATRIX_IN}")
    print(f"  Output 1 : {CANDIDATES_OUT}")
    print(f"  Output 2 : {HITS_OUT}")
    print("")

    evaluated  = 0
    candidates = 0
    total_hits = 0

    candidates_buffer = []
    hits_buffer       = []

    with open(MATRIX_IN, 'r') as matrix_f, \
         open(CANDIDATES_OUT, 'w') as cand_f, \
         open(HITS_OUT, 'w') as hits_f:

        # Read and write header
        header_line = matrix_f.readline().strip()
        sample_ids  = header_line.split('\t')[1:]  # skip 'Site' column

        cand_f.write(header_line + "\n")
        hits_f.write("Site\tSample\tGenotype\n")

        for line in matrix_f:
            line = line.strip()
            if not line:
                continue

            parts     = line.split('\t')
            site      = parts[0]
            genotypes = parts[1:]

            evaluated += 1

            # Apply minimum carrier filter before lethal logic
            carrier_count = sum(1 for g in genotypes if is_valid(g) and is_heterozygous(g))
            if carrier_count < MIN_CARRIERS:
                continue

            # Apply lethal candidate logic
            if is_lethal_candidate(genotypes):
                candidates += 1

                # File 1 — full matrix row
                candidates_buffer.append(line)

                # File 2 — one entry per sample carrying HET at this site
                for sample_id, genotype in zip(sample_ids, genotypes):
                    if is_valid(genotype) and is_heterozygous(genotype):
                        hits_buffer.append(f"{site}\t{sample_id}\t{genotype}")
                        total_hits += 1

            # Flush buffers periodically
            if len(candidates_buffer) >= CHUNK_SIZE:
                cand_f.write("\n".join(candidates_buffer) + "\n")
                candidates_buffer = []

            if len(hits_buffer) >= CHUNK_SIZE:
                hits_f.write("\n".join(hits_buffer) + "\n")
                hits_buffer = []

            if evaluated % 1_000_000 == 0:
                log(f"  Evaluated {evaluated:,} sites | "
                    f"Candidates so far: {candidates:,}")

        # Write remaining buffers
        if candidates_buffer:
            cand_f.write("\n".join(candidates_buffer) + "\n")
        if hits_buffer:
            hits_f.write("\n".join(hits_buffer) + "\n")

    return evaluated, candidates, total_hits


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def main():
    # Create log directory
    os.makedirs(LOG_DIR, exist_ok=True)

    print("=============================================")
    print(" Step 3 — Filter Lethal Candidates")
    print("=============================================")
    print(f"Threshold        : {THRESHOLD}")
    print(f"Min carriers     : {MIN_CARRIERS}")
    print(f"Matrix input     : {MATRIX_IN}")
    print(f"Candidates output: {CANDIDATES_OUT}")
    print(f"Hits output      : {HITS_OUT}")
    print(f"Timestamp        : {datetime.datetime.now()}")
    print("")

    validate_paths()

    evaluated, candidates, total_hits = filter_lethal_candidates()

    print("")
    print("=============================================")
    print(" Step 3 Complete")
    print("=============================================")
    print(f"Sites evaluated      : {evaluated:,}")
    print(f"Min carriers filter  : {MIN_CARRIERS}")
    print(f"Lethal candidates    : {candidates:,} "
          f"({100*candidates/evaluated:.2f}% of evaluated sites)")
    print(f"Total HET hits       : {total_hits:,}")
    print(f"Candidates saved to  : {CANDIDATES_OUT}")
    print(f"Hits saved to        : {HITS_OUT}")
    print(f"Timestamp            : {datetime.datetime.now()}")
    print("=============================================")


if __name__ == "__main__":
    main()
