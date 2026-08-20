#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Launch the CNVkit Snakemake pipeline on SLURM
# =============================================================================
# Prerequisites
#   1. 00_prepare.sh has been run successfully (symlinks, refFlat in place).
#   2. Snakemake is available in the current shell.
#      If snakemake is installed inside the cnvkit conda environment:
#        conda activate cnvkit && bash run_pipeline.sh
#      If snakemake has its own environment:
#        conda activate snakemake && bash run_pipeline.sh
#      Individual CNVkit jobs are activated via 'conda run -n cnvkit' inside
#      the Snakefile, so the snakemake environment does NOT need CNVkit.
#   3. Your SLURM profile is at ~/.config/snakemake/<PROFILE>/config.yaml
#      and maps Snakemake's resources (mem_mb, runtime) to SBATCH directives.
#
# Usage
#   bash run_pipeline.sh [options]
#
# Options
#   --dry-run  | -n        Print jobs without submitting (DAG check)
#   --profile  NAME        SLURM profile name  (default: slurm)
#   --jobs     N           Max concurrent cluster jobs  (default: 100)
#   --until    RULE        Run up to (and including) a specific rule
#   --forcerun RULE        Force-rerun a rule even if outputs are up to date
#   --touch                Touch all output files without running (mark done)
#   --unlock               Unlock a stale Snakemake lock
# =============================================================================

set -euo pipefail

# ─── Locate this script's directory so paths work from any cwd ───────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEFILE="${SCRIPT_DIR}/Snakefile"
CONFIG="${SCRIPT_DIR}/config.yaml"

# ─── Defaults ────────────────────────────────────────────────────────────────
PROFILE="slurm"
MAX_JOBS=100
DRY_RUN=""
EXTRA_ARGS=()

# ─── Parse arguments ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run|-n)
            DRY_RUN="--dry-run --quiet"
            ;;
        --profile)
            PROFILE="$2"; shift
            ;;
        --jobs|-j)
            MAX_JOBS="$2"; shift
            ;;
        --until|--forcerun|--touch|--unlock)
            # Pass through to snakemake as-is (with optional next arg)
            EXTRA_ARGS+=("$1")
            if [[ $# -gt 1 && ! "$2" == --* ]]; then
                EXTRA_ARGS+=("$2"); shift
            fi
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Usage: bash run_pipeline.sh [--dry-run] [--profile PROFILE] [--jobs N]" >&2
            exit 1
            ;;
    esac
    shift
done

# ─── Validate snakemake is available ─────────────────────────────────────────
if ! command -v snakemake &>/dev/null; then
    echo "ERROR: snakemake not found in PATH." >&2
    echo "  Activate the environment that contains snakemake first:" >&2
    echo "    conda activate snakemake # if snakemake has its own env" >&2
    exit 1
fi

SM_VERSION=$(snakemake --version 2>/dev/null || echo "unknown")
echo "==> Snakemake ${SM_VERSION}"
echo "==> Snakefile:  ${SNAKEFILE}"
echo "==> Config:     ${CONFIG}"
echo "==> Profile:    ${PROFILE}"
echo "==> Max jobs:   ${MAX_JOBS}"
[[ -n "${DRY_RUN}" ]] && echo "==> MODE:       DRY RUN (no jobs submitted)"
echo ""

# ─── Run ─────────────────────────────────────────────────────────────────────
snakemake \
    --snakefile  "${SNAKEFILE}" \
    --configfile "${CONFIG}" \
    --profile    "${PROFILE}" \
    --jobs       "${MAX_JOBS}" \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \    
    ${DRY_RUN} \
    "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"

# ─── Exit message ────────────────────────────────────────────────────────────
if [[ -z "${DRY_RUN}" ]]; then
    echo ""
    echo "==> Pipeline launched.  Monitor with:"
    echo "      squeue -u \$USER"
    echo "      tail -f <OUTDIR>/logs/<rule>/<sample>.log"
fi
