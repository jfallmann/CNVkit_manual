# CNVkit_manual
Run CNVkit after Sarek for specified contrasts, then summarize the called
copy-number states (tables, PDF/PNG plots and interactive HTML).

# 1. Copy files into WorkDir/
# profiles/, scripts/ and envs/ must come along: run_pipeline.sh defaults to
# profiles/slurm (per-rule resources), rule genelist invokes
# scripts/gene_calls.py, and the summary rule uses the managed conda env
# envs/cnv_summary.yaml.
cp -r 00_prepare.sh config.yaml Snakefile run_pipeline.sh profiles scripts envs ${WorkDir}/
cd ${WorkDir}/

# 2. One-time setup (symlinks + directory tree)
# 00_prepare.sh takes three arguments:
#   bash 00_prepare.sh <PROJECT_DIR> <OUTPUT_DIR> <BAM_DIR>
#     PROJECT_DIR  Sarek results root; BAMs are expected at
#                  {PROJECT_DIR}/preprocessing/recalibrated/{sample}/{sample}.recal.bam
#     OUTPUT_DIR   pipeline output root; everything lands in {OUTPUT_DIR}/cnvkit_manual
#     BAM_DIR      Sarek recalibrated BAM root (source of the symlinks)
# It symlinks the sample BAMs into cnvkit_manual/bams/ as {sample}.bam for
# inspection; the pipeline itself reads the Sarek layout directly (see bam_dir
# below) and does not use those symlinks.
bash 00_prepare.sh ${PROJECT_DIR} ${OUTPUT_DIR} ${BAM_DIR}

# 3. Export the config paths (or replace them with literals in config.yaml)
# config.yaml uses ${PROJECT_DIR} and ${OUTPUT_DIR} environment variables.
export PROJECT_DIR=/path/to/sarek/results
export OUTPUT_DIR=/path/to/outputs
# bam_dir (optional): set it in config.yaml to read BAMs from an alternative
# Sarek root instead of {project_dir}/preprocessing/recalibrated. It must use
# the same {bam_dir}/{sample}/{sample}.recal.bam layout.

# 4. Download refFlat if missing (see output of step 2)
# Must match cnvkit_manual/reference/ — the path config.yaml's refflat: uses.
# (00_prepare.sh fetches the ENCODE blacklist itself.)
wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz \
    | gunzip > cnvkit_manual/reference/hg38.refFlat.txt

# 5. Dry-run to check the DAG, then launch
conda activate snakemake   # (or whichever env has snakemake)
bash run_pipeline.sh --dry-run
bash run_pipeline.sh
# For a single local run (no SLURM): bash run_pipeline.sh --executor local

# ── Environments ────────────────────────────────────────────────────────────
# CNVkit jobs run in the existing conda env named in config.yaml (conda_env,
# default "cnvkit") via `conda run -n cnvkit`. The summary rule instead uses a
# Snakemake-managed conda env (envs/cnv_summary.yaml, activated with
# --use-conda) containing R and the plotting packages. Keep both available.

# ── Outputs ─────────────────────────────────────────────────────────────────
# Per-comparison results: {outdir}/<mode>/<sample>/<sample>.* with comparison
# IDs "<mode>__<sample>" (e.g. vs_reference__BC139_naive,
# pairwise_bc139_vs_bc1__BC139_naive). Pairwise modes live in their own
# pairwise_<name>/ directories.
# Summary analysis: {outdir}/summary/ — scored segment/region tables (CSV),
# gene matrices (TSV), Manhattan/heatmap/trajectory plots (PDF/PNG) and
# interactive scatters (HTML) under cn_scatters/, temporal/ and conserved/.

# ── Configuration knobs ─────────────────────────────────────────────────────
# ploidy/purity precedence (highest first):
#   comparison_overrides.<sample>.<mode> > comparison_overrides.<sample>
#   > mode_defaults.<mode> > global cnvkit.ploidy/purity
# analysis.* controls the summary: thresholds, plots, oncogenes,
# resistance_samples, temporal_clones, conserved_comparisons,
# scatter_comparisons.

# ── Interpretation notes ────────────────────────────────────────────────────
# - Missing data is NA, never treated as neutral.
# - The confidence column is a heuristic score, not a probability.
# - overlaps_vs_reference in pairwise regions means genomic overlap with the
#   vs_reference regions, not ancestry or resistance specificity.
# - temporal_clones timeline order is the analysis order, not chronological.
# - Called CN is relative to each contrast's baseline (reference/ploidy).

# ── Testing status ──────────────────────────────────────────────────────────
# The BAM-based steps (coverage, references) and SLURM execution have not been
# tested locally; only the summary workflow was exercised end-to-end.
