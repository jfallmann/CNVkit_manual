#!/usr/bin/env bash
# =============================================================================
# 00_prepare.sh — CNVkit WGS manual pipeline: symlink setup & validation
# =============================================================================
# Run this script ONCE before launching the Snakemake pipeline.
# It creates the output directory tree, symlinks the 16 unique biological
# sample BAMs into cnvkit_manual/bams/, and validates key inputs.
#
# Safe to run from any directory — all paths are absolute.
#
# Usage:
#   bash 00_prepare.sh
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# PATHS — edit here if the project moves
# ─────────────────────────────────────────────────────────────────────────────
PROJECT_DIR=$1
OUTPUT_DIR=$2
BAM_DIR=$3

# All pipeline output lands here
OUTDIR="${OUTPUT_DIR}/cnvkit_manual"

FASTA="/resources/references/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
REFFLAT="${OUTDIR}/reference/hg38.refFlat.txt"

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLES
# 16 unique biological samples.
# Bulk_sensitive_1 … Bulk_sensitive_7 are Sarek somatic-pairing aliases —
# all map to the same FASTQ and produce identical BAMs; they are intentionally
# excluded here.
# ─────────────────────────────────────────────────────────────────────────────
SAMPLES=(
    Bulk_sensitive
    Bulk_resistant
    BC1_naive       BC1_naive_1
    BC3_naive       BC3_naive_1
    BC217_naive     BC217_naive_1
    BC217_resistant BC217_resistant_1
    BC139_naive     BC139_naive_1
    BC139_resistant BC139_resistant_1
    BC1444_naive    BC1444_naive_1
)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Directory structure
# ─────────────────────────────────────────────────────────────────────────────
echo "==> Creating output directory structure under: ${OUTDIR}"
mkdir -p \
    "${OUTDIR}/bams" \
    "${OUTDIR}/bins" \
    "${OUTDIR}/coverage" \
    "${OUTDIR}/references" \
    "${OUTDIR}/vs_reference" \
    "${OUTDIR}/vs_normal" \
    "${OUTDIR}/resources" \
    "${OUTDIR}/logs/coverage" \
    "${OUTDIR}/logs/vs_reference" \
    "${OUTDIR}/logs/vs_normal"
echo "    OK"

# ─────────────────────────────────────────────────────────────────────────────
# 2. Sanity-check FASTA
# ─────────────────────────────────────────────────────────────────────────────
if [[ ! -f "${FASTA}" ]]; then
    echo "ERROR: FASTA not found: ${FASTA}" >&2
    exit 1
fi
echo "==> FASTA verified: ${FASTA}"

# ─────────────────────────────────────────────────────────────────────────────
# 3. BAM / BAI symlinks
# ─────────────────────────────────────────────────────────────────────────────
echo "==> Symlinking recalibrated BAMs into ${OUTDIR}/bams/"
N_MISSING=0

for SAMPLE in "${SAMPLES[@]}"; do
    SRC_BAM="${BAM_DIR}/${SAMPLE}/${SAMPLE}.recal.bam"
    SRC_BAI="${SRC_BAM}.bai"
    DST_BAM="${OUTDIR}/bams/${SAMPLE}.bam"
    DST_BAI="${OUTDIR}/bams/${SAMPLE}.bam.bai"

    if [[ ! -f "${SRC_BAM}" ]]; then
        echo "  [MISSING BAM] ${SRC_BAM}" >&2
        N_MISSING=$((N_MISSING + 1))
        continue
    fi

    ln -sfn "${SRC_BAM}" "${DST_BAM}"

    if [[ -f "${SRC_BAI}" ]]; then
        ln -sfn "${SRC_BAI}" "${DST_BAI}"
        echo "    ${SAMPLE}  ✓"
    else
        echo "  [MISSING BAI] ${SRC_BAI}" >&2
        echo "  --> Run: samtools index ${SRC_BAM}" >&2
        N_MISSING=$((N_MISSING + 1))
    fi
done

if (( N_MISSING > 0 )); then
    echo ""
    echo "  WARNING: ${N_MISSING} file(s) missing — see messages above before continuing."
else
    echo "    All BAM/BAI symlinks OK."
fi

# ─────────────────────────────────────────────────────────────────────────────
# 4. refFlat gene annotation
# ─────────────────────────────────────────────────────────────────────────────
# The pipeline requires a chr-prefixed hg38 refFlat.txt for gene annotation
# (--annotate in cnvkit reference).  This file uses UCSC contig naming and
# matches Homo_sapiens_assembly38.fasta.
#
# If not already present, download with:
#
#   wget -qO- \
#     https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz \
#     | gunzip > "${REFFLAT}"
#
if [[ -f "${REFFLAT}" ]]; then
    LINES=$(wc -l < "${REFFLAT}")
    echo "==> refFlat found: ${REFFLAT}  (${LINES} lines)"
else
    echo ""
    echo "==> [ACTION REQUIRED] refFlat not found at:"
    echo "      ${REFFLAT}"
    echo ""
    echo "    Download it before running the pipeline:"
    echo ""
    echo "      wget -qO- \\"
    echo "        https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz \\"
    echo "        | gunzip > \"${REFFLAT}\""
    echo ""
fi

# ─────────────────────────────────────────────────────────────────────────────
# 5. Summary
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " Setup complete.  Output directory: ${OUTDIR}"
echo "============================================================"
echo ""
echo "Next steps:"
echo "  1. Download refFlat.txt if not already present (see above)."
echo "  2. Review config.yaml — adjust SLURM resource limits if needed."
echo "  3. Dry-run to verify the DAG:"
echo "       bash run_pipeline.sh --dry-run"
echo "  4. Launch:"
echo "       bash run_pipeline.sh"
echo ""
