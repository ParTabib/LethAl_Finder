#!/bin/bash
# =============================================================================
# Step 1 — Collect Heterozygous Sites
# =============================================================================
# Scans all gzipped HEARTY output files for sites where at least one sample
# is heterozygous at the given threshold. Outputs a sorted unique list of
# site IDs (chr_pos) to het_sites_<threshold>.txt in the output directory.
#
# Parallelizes across --cores processes — each core processes one sample
# file simultaneously, then results are merged and deduplicated.
#
# This script is called by the LethAl Finder master script and should not
# be run directly unless you know what you are doing.
#
# Usage:
#   bash 01_collect_het_sites.sh <threshold> <input_dir> <output_dir> <cores>
#   e.g. bash 01_collect_het_sites.sh 0.25 /data/hearty /results 4
# =============================================================================

# --- Input validation --------------------------------------------------------
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    echo "[ERROR] Missing arguments."
    echo "        Usage: bash 01_collect_het_sites.sh <threshold> <input_dir> <output_dir> <cores>"
    exit 1
fi

THRESHOLD=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
CORES=$4
LOG_DIR="${OUTPUT_DIR}/logs"
OUTPUT_FILE="${OUTPUT_DIR}/het_sites_${THRESHOLD}.txt"
STATUS_COL="Status_${THRESHOLD}"
TEMP_DIR="${OUTPUT_DIR}/temp_het_sites"

echo "============================================="
echo " Step 1 — Collect Heterozygous Sites"
echo "============================================="
echo "Threshold        : $THRESHOLD"
echo "Target column    : $STATUS_COL"
echo "Input directory  : $INPUT_DIR"
echo "Output directory : $OUTPUT_DIR"
echo "Output file      : $OUTPUT_FILE"
echo "Cores            : $CORES"
echo "Timestamp        : $(date)"
echo ""

# --- Validate directories ----------------------------------------------------
if [ ! -d "$INPUT_DIR" ]; then
    echo "[ERROR] Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "[ERROR] Output directory not found: $OUTPUT_DIR"
    exit 1
fi

# --- Create log and temp directories -----------------------------------------
mkdir -p "$LOG_DIR"
mkdir -p "$TEMP_DIR"

# --- Validate HEARTY files exist ---------------------------------------------
FILES=("${INPUT_DIR}"/*.basecall.txt.gz)
if [ ${#FILES[@]} -eq 0 ]; then
    echo "[ERROR] No .basecall.txt.gz files found in: $INPUT_DIR"
    exit 1
fi
echo "Sample files detected: ${#FILES[@]}"
echo ""

# --- Detect Status column number from header of first file -------------------
echo "Detecting column position of '${STATUS_COL}' from file header..."

COL_NUM=$(zcat "${FILES[0]}" | head -1 | tr '\t' '\n' | grep -n "^${STATUS_COL}$" | cut -d: -f1)

if [ -z "$COL_NUM" ]; then
    echo "[ERROR] Column '${STATUS_COL}' not found in file header."
    echo "        Available threshold columns:"
    zcat "${FILES[0]}" | head -1 | tr '\t' '\n' | grep "Status_"
    echo ""
    echo "        Re-run with one of the above threshold values."
    exit 1
fi

echo "  -> '${STATUS_COL}' is column number: $COL_NUM"
echo ""

# --- Worker function — process one sample file -------------------------------
process_sample() {
    local FILE=$1
    local COL_NUM=$2
    local TEMP_DIR=$3
    local SAMPLE
    SAMPLE=$(basename "$FILE" .basecall.txt.gz)
    local TMP_OUT="${TEMP_DIR}/${SAMPLE}_hetsites.txt"

    echo "  Processing: $SAMPLE"

    # Decompress and scan:
    # - Skip the header line (NR > 1)
    # - Check if the status column starts with HET
    # - Print chr_pos as the site ID
    # - Deduplicate per file immediately to keep temp file small
    zcat "$FILE" | awk -v col="$COL_NUM" '
        NR > 1 && $col ~ /^HET/ {
            print $1"_"$2
        }
    ' | sort -u > "$TMP_OUT"
}

export -f process_sample

# --- Scan all files in parallel ----------------------------------------------
echo "Scanning all sample files for HET sites using $CORES core(s)..."
echo ""

# Use GNU parallel if available, otherwise fall back to xargs
if command -v parallel &> /dev/null; then
    printf '%s\n' "${FILES[@]}" | \
        parallel -j "$CORES" process_sample {} "$COL_NUM" "$TEMP_DIR"
else
    printf '%s\n' "${FILES[@]}" | \
        xargs -P "$CORES" -I{} bash -c 'process_sample "$@"' _ {} "$COL_NUM" "$TEMP_DIR"
fi

echo ""

# --- Final deduplication across all samples ----------------------------------
echo "Merging and deduplicating final site list across all samples..."

TMP_FILES=("${TEMP_DIR}"/*_hetsites.txt)
sort -u "${TMP_FILES[@]}" > "$OUTPUT_FILE"

# Cleanup temp files
rm -rf "$TEMP_DIR"

SITE_COUNT=$(wc -l < "$OUTPUT_FILE")

echo ""
echo "============================================="
echo " Step 1 Complete"
echo "============================================="
echo "Total unique HET sites found : $SITE_COUNT"
echo "Output saved to              : $OUTPUT_FILE"
echo "Timestamp                    : $(date)"
echo "============================================="
