# LethAl_Finder

A computational pipeline for screening population-scale genomic data to identify lethal recessive alleles — genomic sites where a heterozygous call exists in the population but one or both corresponding homozygotes are completely absent, suggesting lethality in double dose.

Developed as part of a special course project under the **Biodiversity and Extinction Research Group** at the **Technical University of Denmark (DTU)**, supervised by **Mick Westbury**.

---

## How It Works

The pipeline operates in three steps:

**Step 1 — Collect Heterozygous Sites**
Scans all HEARTY output files across the population and collects every genomic site where at least one individual is heterozygous at the given frequency threshold. Outputs a sorted unique list of site IDs.

**Step 2 — Build Genotype Matrix**
Using the HET site list from Step 1, extracts the base call for every individual at every HET site and assembles a full cross-population genotype matrix. Sites where fewer than a user-defined percentage of individuals have data are filtered out.

**Step 3 — Filter Lethal Candidates**
Applies the lethal recessive selection logic row by row. For each HET site, the allele pair is extracted (e.g. `AC` → `A` and `C`) and both corresponding homozygotes (`AA` and `CC`) are looked up across the population. If either homozygote is absent, the site is flagged as a lethal candidate.

---

## Input

LethAl_Finder is designed to work with the output of [HEARTY](https://github.com/BiodiversityExtinction/HEARTY) — a toolkit for profiling heterozygosity patterns from sequencing data. HEARTY generates per-sample `.basecall.txt.gz` files containing base calls and heterozygosity status at each genomic site across one or more frequency thresholds.

---

## Important: Pre-filtering Recommendations

Before running LethAl_Finder, we strongly recommend filtering your input data to remove genomic regions that can produce false lethal candidates:

### 1. Sex Chromosomes

Sex chromosomes behave fundamentally differently from autosomes. Males are hemizygous on the X chromosome — they carry only one copy — so X-linked sites will always appear to have a missing homozygote and will be incorrectly flagged as lethal candidates. The same applies to Y-linked and Z/W-linked scaffolds in non-mammalian species.

**Recommendation:** Use [SATC (Sex Assignment Through Coverage)](https://github.com/popgenDK/SATC) to identify sex-linked scaffolds in your population before running LethAl_Finder. SATC uses depth-of-coverage from BAM files and PCA-based clustering to identify sex chromosomes and sex-linked scaffolds in non-model organisms. Once identified, exclude those scaffolds from your HEARTY input using the `-r` or `-f` regions argument in HEARTY.

### 2. Mitochondrial DNA

Mitochondrial DNA is present in many copies per cell, is maternally inherited, and is effectively haploid — all sites will appear homozygous or as false heterozygotes due to heteroplasmy. Mitochondrial scaffolds (typically a single small circular chromosome) should be excluded.

**Recommendation:** Ensure the mitochondrial scaffold is excluded from your input files before running LethAl_Finder.

### 3. Unplaced Scaffolds (NW_ / Unanchored Contigs)

Unplaced scaffolds are genomic sequences that could not be confidently assigned to a chromosome during genome assembly. They may produce false lethal candidates due to lower assembly confidence.

**Recommendation:** Ensure unplaced scaffolds (e.g. `NW_*` in NCBI assemblies) are excluded from your input files before running LethAl_Finder, and include only chromosome-level scaffolds.

---

## Requirements

- Python 3.6+
- bash, awk, zcat (standard on Linux/macOS)
- GNU parallel (optional but recommended for Step 1 parallelization — falls back to `xargs` if not available)

No external Python packages are required — the pipeline uses only the standard library.

---

## Installation

```bash
git clone https://github.com/ParTabib/LethAl_Finder.git
cd LethAl_Finder
```

---

## Usage

```bash
python3 00_lethal_finder.py \
    --input-dir  /path/to/hearty/files \
    --output-dir /path/to/results \
    --threshold  0.25
```

### All arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--input-dir` | Yes | — | Directory containing HEARTY `.basecall.txt.gz` files |
| `--output-dir` | Yes | — | Directory to write all output files into |
| `--threshold` | Yes | — | Frequency threshold (e.g. `0.05` or `0.25`) |
| `--min-coverage` | No | `15` | Min % of individuals with non-NA data for a site to be kept |
| `--cores` | No | `1` | Number of CPU cores for parallelization |
| `--sample-sheet` | No | Auto-detect | TSV file mapping filenames to sample IDs |
| `--skip-step1` | No | False | Skip Step 1 — use existing HET sites file |
| `--skip-step2` | No | False | Skip Step 2 — use existing genotype matrix |

### Examples

```bash
# Basic run with default settings
python3 00_lethal_finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25

# Use 4 cores and 20% minimum coverage
python3 00_lethal_finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --cores 4 --min-coverage 20

# Resume from Step 3 (Steps 1 and 2 already completed)
python3 00_lethal_finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --skip-step1 --skip-step2

# Provide custom sample names via sample sheet
python3 00_lethal_finder.py --input-dir /data/hearty --output-dir /results --threshold 0.25 --sample-sheet samples.tsv
```

### Sample sheet format

If your filenames are not meaningful, provide a two-column TSV mapping filenames to sample IDs:

```
sample01.basecall.txt.gz    Simba
sample02.basecall.txt.gz    Nala
sample03.basecall.txt.gz    Mufasa
```

If no sample sheet is provided, sample IDs are derived automatically by stripping the common suffix from all filenames.

---

## Output Files

All output files are written to `--output-dir`:

| File | Description |
|---|---|
| `het_sites_<threshold>.txt` | Sorted unique list of all HET site IDs across the population |
| `master_genotype_matrix_<threshold>.tsv` | Full cross-population genotype matrix (rows = sites, columns = individuals) |
| `lethal_candidates_<threshold>.tsv` | Subset of the matrix containing only lethal candidate sites |
| `lethal_hits_<threshold>.tsv` | One row per individual carrying a HET genotype at a lethal site |
| `logs/` | Log files for each step |

### lethal_hits format

```
Site                    Sample    Genotype
NC_056679.1_11381       Simba     AT
NC_056679.1_11381       Nala      AT
NC_056679.1_11390       Simba     AC
```

---

## Resource Requirements

Memory and runtime scale with the number of HET sites and the number of individuals:

| Threshold | HET sites (example) | RAM per core | Recommended RAM |
|---|---|---|---|
| 0.25 | ~11M | ~1 GB | 16G |
| 0.05 | ~368M | ~30 GB | 128G |

**Runtime estimate with 4 cores (29 individuals, 7.8 GB compressed files):**

| Step | Time |
|---|---|
| Step 1 | ~4 hours |
| Step 2 | ~8 hours |
| Step 3 | ~30 minutes |
| **Total** | **~12.5 hours** |

For HPC clusters, the pipeline is compatible with any workload manager (SLURM, PBS, LSF). Simply wrap the run command in your cluster submission script and request resources accordingly.

---

## Parallelization

Steps 1 and 2 support multi-core parallelization via `--cores`:

- **Step 1** — sample files scanned simultaneously
- **Step 2 Phase 2** — sample files processed simultaneously (each core holds its own copy of the HET set in memory)
- **Step 2 Phase 3** — temporary files sorted simultaneously
- **Step 2 Phase 4** — merge is inherently sequential

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Acknowledgements

Developed as a special course project under the Biodiversity and Extinction Research Group, Technical University of Denmark (DTU).
Supervised by Mick Westbury.
