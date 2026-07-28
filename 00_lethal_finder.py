#!/usr/bin/env python3
# =============================================================================
# LethAl Finder — Lethal Recessive Allele Screener
# =============================================================================
# Screens HEARTY genomic output files across a population to identify
# lethal recessive mutations — genomic sites where the heterozygous form
# exists in the population but one or both corresponding homozygotes are
# completely absent, suggesting lethality in double dose.
#
# Pipeline steps:
#   Step 1 — Collect HET sites     (01_collect_het_sites.sh)
#   Step 2 — Build genotype matrix (02_build_matrix.py)
#   Step 3 — Filter lethal candidates (03_filter_lethal.py)
#   Step 4 — Filter by genomic region  (04_filter_regions.py, optional)
#   Step 5 — Annotate variants          (05_annotate_variants.py, optional)
#
# Usage:
#   python3 00_LethAl_Finder.py \
#       --input-dir  /path/to/hearty/files \
#       --output-dir /path/to/results \
#       --threshold  0.25 \
#       --min-coverage 15 \
#       --cores 4
#
# Optional:
#   --sample-sheet samples.tsv   Two-column TSV: filename TAB sample_id
#   --skip-step1                 Skip Step 1 (HET sites already collected)
#   --skip-step2                 Skip Step 2 (matrix already built)
#
# Output files (all in --output-dir):
#   het_sites_<threshold>.txt                   All unique HET site IDs
#   master_genotype_matrix_<threshold>.tsv      Full cross-population matrix
#   lethal_candidates_<threshold>.tsv           Lethal candidate sites
#   lethal_hits_<threshold>.tsv                 Per-sample carrier events
#   lethal_candidates_<threshold>_regional.tsv   Region-filtered candidates (if --regions used)
#   lethal_candidates_<threshold>_annotated.tsv  Annotated candidates (if --fasta used)
#   logs/                                        All step logs
# =============================================================================

import os
import sys
import gzip
import math
import argparse
import datetime
import subprocess


# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="LethAl Finder — Lethal Recessive Allele Screener",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python3 00_LethAl_Finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25
  python3 00_LethAl_Finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --cores 4
  python3 00_LethAl_Finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --skip-step1
  python3 00_LethAl_Finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --sample-sheet samples.tsv
  python3 00_LethAl_Finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --regions genome.gff --feature exon
  python3 00_LethAl_Finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --regions genome.gff --fasta ref.fna
    """
)
parser.add_argument("--input-dir",    required=True,
                    help="Directory containing HEARTY .basecall.txt.gz files")
parser.add_argument("--output-dir",   required=True,
                    help="Directory to write all output files into")
parser.add_argument("--threshold",    required=True,
                    help="Frequency threshold to use (e.g. 0.25 or 0.05)")
parser.add_argument("--min-coverage", type=float, default=15.0,
                    help="Minimum %% of samples with non-NA data for a site "
                         "to be kept (default: 15)")
parser.add_argument("--cores",        type=int, default=1,
                    help="Number of CPU cores for parallelization (default: 1)")
parser.add_argument("--sample-sheet", default=None,
                    help="Optional TSV file mapping filenames to sample IDs "
                         "(columns: filename, sample_id)")
parser.add_argument("--skip-step1",   action="store_true",
                    help="Skip Step 1 — use existing het_sites file in output directory")
parser.add_argument("--skip-step2",   action="store_true",
                    help="Skip Step 2 — use existing genotype matrix in output directory")
parser.add_argument("--min-carriers", type=int, default=1,
                    help="Minimum number of individuals carrying the HET genotype at a site "
                         "for it to be considered a lethal candidate (default: 1)")
parser.add_argument("--skip-step3",   action="store_true",
                    help="Skip Step 3 — use existing lethal candidates file in output directory")
parser.add_argument("--skip-step4",   action="store_true",
                    help="Skip Step 4 — use existing regional candidates file in output directory")
parser.add_argument("--regions",      default=None,
                    help="Optional BED/GFF/GFF3/GTF file to filter lethal candidates "
                         "to specific genomic regions (e.g. exons)")
parser.add_argument("--feature",      type=str, default="exon",
                    help="Feature type to extract from GFF/GTF files "
                         "(default: exon). Ignored for BED files.")
parser.add_argument("--min-depth",    type=int, default=10,
                    help="Minimum read depth (totDepth) required at a site for a sample's "
                         "base call to be included in the matrix. Calls below this threshold "
                         "are written as NA (default: 10)")
parser.add_argument("--max-alleles",  type=int, default=None,
                    help="Maximum number of unique alleles allowed across all individuals "
                         "at a site. Sites exceeding this limit are excluded as potential "
                         "noise or mapping artifacts (default: no limit)")
parser.add_argument("--fasta",        default=None,
                    help="Reference FASTA file (.fna) with .fai index — required to run "
                         "Step 5 variant annotation (requires --regions GFF3)")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

INPUT_DIR    = os.path.abspath(args.input_dir)
OUTPUT_DIR   = os.path.abspath(args.output_dir)
THRESHOLD    = str(float(args.threshold))
MIN_COVERAGE = args.min_coverage
CORES        = args.cores
SAMPLE_SHEET = args.sample_sheet
SKIP_STEP1   = args.skip_step1
SKIP_STEP2   = args.skip_step2
SKIP_STEP3   = args.skip_step3
SKIP_STEP4   = args.skip_step4
MIN_CARRIERS = args.min_carriers
REGIONS      = args.regions
FEATURE      = args.feature
MIN_DEPTH    = args.min_depth
MAX_ALLELES  = args.max_alleles
FASTA        = args.fasta

LOG_DIR        = os.path.join(OUTPUT_DIR, "logs")
HET_SITES_FILE = os.path.join(OUTPUT_DIR, f"het_sites_{THRESHOLD}.txt")
MATRIX_FILE    = os.path.join(OUTPUT_DIR, f"master_genotype_matrix_{THRESHOLD}.tsv")
CANDIDATES_OUT = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}.tsv")
HITS_OUT       = os.path.join(OUTPUT_DIR, f"lethal_hits_{THRESHOLD}.tsv")
REGIONAL_OUT   = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}_regional.tsv")
ANNOTATED_OUT  = os.path.join(OUTPUT_DIR, f"lethal_candidates_{THRESHOLD}_annotated.tsv")

# Locate step scripts relative to this master script
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
STEP1_SCRIPT = os.path.join(SCRIPT_DIR, "01_collect_het_sites.sh")
STEP2_SCRIPT = os.path.join(SCRIPT_DIR, "02_build_matrix.py")
STEP3_SCRIPT = os.path.join(SCRIPT_DIR, "03_filter_lethal.py")
STEP4_SCRIPT = os.path.join(SCRIPT_DIR, "04_filter_regions.py")
STEP5_SCRIPT = os.path.join(SCRIPT_DIR, "05_annotate_variants.py")

# Track timing per step
STEP_TIMES = {}


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log(msg: str):
    """Print message with timestamp."""
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")
    sys.stdout.flush()


def banner(title: str):
    """Print a section banner."""
    print("")
    print("=" * 55)
    print(f" {title}")
    print("=" * 55)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def get_sample_files() -> list:
    """Get sorted list of all HEARTY files in the input directory."""
    return sorted([
        os.path.join(INPUT_DIR, f)
        for f in os.listdir(INPUT_DIR)
        if f.endswith(".basecall.txt.gz")
    ])


def validate_threshold_in_file(file_path: str):
    """Check that threshold columns exist in the HEARTY file header."""
    base_col   = f"Base_{THRESHOLD}"
    status_col = f"Status_{THRESHOLD}"
    with gzip.open(file_path, 'rt') as f:
        header = f.readline().strip().split('\t')
    if base_col not in header or status_col not in header:
        print(f"[ERROR] Columns '{base_col}' or '{status_col}' not found "
              f"in file header.")
        print(f"        Available threshold columns:")
        for col in header:
            if col.startswith("Status_"):
                print(f"          {col.replace('Status_', '')}")
        sys.exit(1)


def validate_inputs():
    """Validate all inputs before starting the pipeline."""
    errors = []

    if not os.path.isdir(INPUT_DIR):
        errors.append(f"Input directory not found: {INPUT_DIR}")

    if not (0 < MIN_COVERAGE <= 100):
        errors.append(f"--min-coverage must be between 0 and 100 "
                      f"(got {MIN_COVERAGE})")

    if CORES < 1:
        errors.append(f"--cores must be at least 1 (got {CORES})")

    if SAMPLE_SHEET and not os.path.isfile(SAMPLE_SHEET):
        errors.append(f"Sample sheet not found: {SAMPLE_SHEET}")

    for script in [STEP1_SCRIPT, STEP2_SCRIPT, STEP3_SCRIPT, STEP4_SCRIPT, STEP5_SCRIPT]:
        if not os.path.isfile(script):
            errors.append(f"Step script not found: {script}")

    sample_files = get_sample_files()
    if not sample_files:
        errors.append(f"No .basecall.txt.gz files found in: {INPUT_DIR}")

    if SKIP_STEP1 and not os.path.isfile(HET_SITES_FILE):
        errors.append(f"--skip-step1 set but HET sites file not found: "
                      f"{HET_SITES_FILE}")

    if SKIP_STEP2 and not os.path.isfile(MATRIX_FILE):
        errors.append(f"--skip-step2 set but matrix file not found: "
                      f"{MATRIX_FILE}")

    if SKIP_STEP3 and not os.path.isfile(CANDIDATES_OUT):
        errors.append(f"--skip-step3 set but candidates file not found: "
                      f"{CANDIDATES_OUT}")

    if SKIP_STEP4 and not os.path.isfile(REGIONAL_OUT):
        errors.append(f"--skip-step4 set but regional candidates file not found: "
                      f"{REGIONAL_OUT}")

    if errors:
        print("\n[ERROR] Validation failed:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)

    # Validate threshold exists in file header
    validate_threshold_in_file(sample_files[0])


# ---------------------------------------------------------------------------
# Resource Estimation
# ---------------------------------------------------------------------------

def estimate_resources():
    """
    Print a resource estimation based on the HET sites file so users
    can verify memory requirements before the heavy work starts.
    """
    banner("Resource Estimation")

    if os.path.isfile(HET_SITES_FILE):
        with open(HET_SITES_FILE) as f:
            n_sites = sum(1 for _ in f)
        # ~67 bytes per site in Python set + hash table overhead (~8 bytes/slot)
        set_size_gb = round((n_sites * 75) / 1e9, 2)
        total_gb    = round(set_size_gb * CORES, 2)

        print(f"  HET sites found      : {n_sites:,}")
        print(f"  RAM per core (set)   : ~{set_size_gb} GB")
        print(f"  Cores requested      : {CORES}")
        print(f"  Total estimated RAM  : ~{total_gb} GB")

        if total_gb > 100:
            print(f"\n  [WARNING] Estimated RAM exceeds 100GB.")
            print(f"            Consider using fewer --cores or a higher "
                  f"--threshold.")
        else:
            print(f"\n  [OK] RAM estimate looks reasonable.")
    else:
        print("  HET sites file not yet available — skipping estimation.")
        print(f"  Cores requested: {CORES}")


# ---------------------------------------------------------------------------
# Step Runners
# ---------------------------------------------------------------------------

def run_step(step_name: str, cmd: list):
    """
    Run a step script as a subprocess.
    Captures timing, checks return code, and exits cleanly on failure.
    """
    log(f"Starting {step_name}...")
    start  = datetime.datetime.now()
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"\n[ERROR] {step_name} failed with exit code {result.returncode}")
        print(f"        Check the output above for details.")
        sys.exit(1)

    elapsed              = datetime.datetime.now() - start
    STEP_TIMES[step_name] = elapsed
    log(f"{step_name} completed in {str(elapsed).split('.')[0]}")


def run_step1():
    """Call 01_collect_het_sites.sh with the required arguments."""
    cmd = [
        "bash", STEP1_SCRIPT,
        THRESHOLD,
        INPUT_DIR,
        OUTPUT_DIR,
        str(CORES)
    ]
    run_step("Step 1", cmd)


def run_step2():
    """Call 02_build_matrix.py with the required arguments."""
    cmd = [
        "python3", STEP2_SCRIPT,
        THRESHOLD,
        INPUT_DIR,
        OUTPUT_DIR,
        "--min-coverage", str(MIN_COVERAGE),
        "--min-depth",    str(MIN_DEPTH),
        "--cores",        str(CORES),
    ]
    if SAMPLE_SHEET:
        cmd += ["--sample-sheet", SAMPLE_SHEET]
    run_step("Step 2", cmd)


def run_step3():
    """Call 03_filter_lethal.py with the required arguments."""
    cmd = [
        "python3", STEP3_SCRIPT,
        THRESHOLD,
        OUTPUT_DIR,
        "--min-carriers", str(MIN_CARRIERS),
    ]
    if MAX_ALLELES is not None:
        cmd += ["--max-alleles", str(MAX_ALLELES)]
    run_step("Step 3", cmd)

def run_step4():
    """Call 04_filter_regions.py with the required arguments."""
    cmd = [
        "python3", STEP4_SCRIPT,
        THRESHOLD,
        OUTPUT_DIR,
        REGIONS,
        "--feature", FEATURE,
        "--cores",   str(CORES),
    ]
    run_step("Step 4", cmd)


def run_step5():
    """Call 05_annotate_variants.py with the required arguments."""
    cmd = [
        "python3", STEP5_SCRIPT,
        THRESHOLD,
        OUTPUT_DIR,
        REGIONS,
        "--fasta",  FASTA,
        "--cores",  str(CORES),
    ]
    run_step("Step 5", cmd)


# ---------------------------------------------------------------------------
# Final Summary
# ---------------------------------------------------------------------------

def print_summary():
    """Print a final summary of all outputs, file sizes, and timing."""
    banner("Pipeline Complete — Summary")

    sample_files = get_sample_files()
    print(f"  Threshold        : {THRESHOLD}")
    print(f"  Min coverage     : {MIN_COVERAGE}%")
    print(f"  Cores used       : {CORES}")
    print(f"  Samples          : {len(sample_files)}")
    print("")

    print("  Output files:")
    for f in [HET_SITES_FILE, MATRIX_FILE, CANDIDATES_OUT, HITS_OUT, REGIONAL_OUT, ANNOTATED_OUT]:
        if os.path.isfile(f):
            size = os.path.getsize(f) / 1e9
            print(f"    {os.path.basename(f):<45} {size:.2f} GB")
    print("")

    print("  Step timing:")
    total = datetime.timedelta()
    for step, elapsed in STEP_TIMES.items():
        print(f"    {step:<10} {str(elapsed).split('.')[0]}")
        total += elapsed
    print(f"    {'Total':<10} {str(total).split('.')[0]}")
    print("")

    print(f"  Logs directory   : {LOG_DIR}")
    print("")


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def main():
    # Create output and log directories
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    banner("LethAl Finder — Lethal Recessive Allele Screener")
    print(f"  Input directory  : {INPUT_DIR}")
    print(f"  Output directory : {OUTPUT_DIR}")
    print(f"  Threshold        : {THRESHOLD}")
    print(f"  Min coverage     : {MIN_COVERAGE}%")
    print(f"  Cores            : {CORES}")
    print(f"  Sample sheet     : "
          f"{SAMPLE_SHEET if SAMPLE_SHEET else 'Auto-detect'}")
    print(f"  Min carriers     : {MIN_CARRIERS}")
    print(f"  Skip Step 1      : {SKIP_STEP1}")
    print(f"  Skip Step 2      : {SKIP_STEP2}")
    print(f"  Skip Step 3      : {SKIP_STEP3}")
    print(f"  Skip Step 4      : {SKIP_STEP4}")
    print(f"  Min depth        : {MIN_DEPTH}")
    print(f"  Max alleles      : {MAX_ALLELES if MAX_ALLELES is not None else 'No limit'}")
    print(f"  Regions file     : {REGIONS if REGIONS else 'Not provided'}")
    print(f"  Feature type     : {FEATURE}")
    print(f"  Reference FASTA  : {FASTA if FASTA else 'Not provided (Step 5 will be skipped)'}")
    print(f"  Timestamp        : {datetime.datetime.now()}")

    # Validate all inputs
    log("Validating inputs...")
    validate_inputs()
    log("  -> All inputs validated.")

    # Step 1 — Collect HET sites
    if SKIP_STEP1:
        banner("Step 1 — Collect Heterozygous Sites")
        log(f"Skipping — using existing file: {HET_SITES_FILE}")
        site_count = sum(1 for _ in open(HET_SITES_FILE))
        log(f"  -> {site_count:,} HET sites in existing file.")
        STEP_TIMES["Step 1"] = datetime.timedelta(0)
    else:
        run_step1()

    # Resource estimation after Step 1
    estimate_resources()

    # Step 2 — Build genotype matrix
    if SKIP_STEP2:
        banner("Step 2 — Build Genotype Matrix")
        log(f"Skipping — using existing file: {MATRIX_FILE}")
        STEP_TIMES["Step 2"] = datetime.timedelta(0)
    else:
        run_step2()

    # Step 3 — Filter lethal candidates
    if SKIP_STEP3:
        banner("Step 3 — Filter Lethal Candidates")
        log(f"Skipping — using existing file: {CANDIDATES_OUT}")
        STEP_TIMES["Step 3"] = datetime.timedelta(0)
    else:
        run_step3()

    # Step 4 — Filter by genomic region (optional)
    if REGIONS:
        if SKIP_STEP4:
            banner("Step 4 — Filter by Genomic Region")
            log(f"Skipping — using existing file: {REGIONAL_OUT}")
            STEP_TIMES["Step 4"] = datetime.timedelta(0)
        else:
            run_step4()

    # Step 5 — Annotate variants (optional — requires --fasta and --regions GFF3)
    if FASTA and REGIONS:
        run_step5()

    # Final summary
    print_summary()


if __name__ == "__main__":
    main()
