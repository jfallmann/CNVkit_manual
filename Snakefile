# =============================================================================
# Snakefile — CNVkit WGS copy-number pipeline
# =============================================================================
#
# Two parallel analyses per sample:
#
#   vs_reference  — sample vs a flat (theoretical, log2 = 0) genome reference.
#                   Reveals CNVs relative to a perfect diploid genome.
#                   Includes Bulk_sensitive itself as a QC check.
#
#   vs_normal     — sample vs Bulk_sensitive as the matched normal.
#                   Removes germline CNVs shared with the baseline.
#
# Pipeline steps (one SLURM job each, all parallelised across samples):
#
#   0. access             — accessible genome regions from the FASTA (wgs -g),
#                           minus blacklisted regions and dropped chromosomes
#   1. autobin            — build genome-wide WGS bins once from the normal BAM
#   2. coverage           — per-sample read-depth in target + antitarget bins
#   3a. flat_reference    — flat reference (from bins only, no coverage)
#   3b. normal_reference  — matched normal reference (from Bulk_sensitive coverage)
#   4. fix                — subtract reference + GC/bias correction  → .cnr
#   5. segment            — circular binary segmentation (CBS)        → .cns
#   6. call               — integer copy-number calls                 → .call.cns
#   7. genemetrics        — per-gene statistics                       → .genemetrics.tsv
#  7b. genelist           — amplified / deleted gene lists            → .amplified.tsv,
#                                                                       .deleted.tsv
#   8. scatter            — genome-wide scatter plot                  → .scatter.png
#   9. diagram            — chromosome arm diagram                    → .diagram.pdf
#  10. heatmap            — multi-sample heatmap per mode             → heatmap.pdf
#
# Run:
#   bash run_pipeline.sh            # full run via SLURM
#   bash run_pipeline.sh --dry-run  # check DAG only
# =============================================================================

import os
import sys
import json
import re

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from workflow_config import (comparisons, expand_path, resolve_call,
                             summary_manifest, summary_outputs)

if not workflow.overwrite_configfiles:
    configfile: os.path.join(workflow.basedir, "config.yaml")

for key in ("project_dir", "outdir", "fasta", "refflat"):
    config[key] = expand_path(config[key])
if config.get("bam_dir"):
    config["bam_dir"] = expand_path(config["bam_dir"])
config["cnvkit"]["access_exclude"] = [
    expand_path(p) for p in config["cnvkit"].get("access_exclude", [])
]
COMPARISONS = comparisons(config)
COMPARISON_MAP = {(r["mode"], r["sample"]): r for r in COMPARISONS}
MANIFEST = summary_manifest(config, COMPARISONS)
REPORT_OUTPUTS = summary_outputs(MANIFEST)
REFERENCE_SAMPLES = sorted({r["reference"] for r in COMPARISONS
                            if r["reference"] not in ("human_reference", config["normal_sample"])})

# ─── Derived constants ────────────────────────────────────────────────────────
PROJECT  = config["project_dir"]
OUTDIR   = config["outdir"]
FASTA    = config["fasta"]
REFFLAT  = config["refflat"]
CONDA    = config["conda_env"]
NORMAL   = config["normal_sample"]
BAM_DIR  = config.get("bam_dir")

ALL_SAMPLES    = config["all_samples"]
VS_REF_SAMPLES = config["vs_ref_samples"]
VS_NOR_SAMPLES = config["vs_normal_samples"]

PAIRWISE_MODES = ["pairwise_" + name for name in (config.get("pairwise_comparisons") or {})]
COMPARISON_MODES   = [r["mode"] for r in COMPARISONS]
COMPARISON_SAMPLES = [r["sample"] for r in COMPARISONS]

# Prefix for every CNVkit call — activates the named conda env without
# requiring the SLURM job script to run in an interactive shell.
CNVKIT = f"conda run --no-capture-output -n {CONDA} cnvkit.py"

# Regions to drop from the accessible genome (ENCODE blacklist, centromeres,
# ...).  Optional — leave the config list empty to keep the whole assembly
# minus its N-runs.
ACCESS_EXCLUDE  = config["cnvkit"].get("access_exclude") or []
ACCESS_MIN_GAP  = config["cnvkit"].get("access_min_gap", 5000)

METHOD      = config["cnvkit"]["method"]
MAPQ        = config["cnvkit"]["min_mapq"]
SEG_METH    = config["cnvkit"]["segment_method"]
COV_THREADS = config["cnvkit"]["coverage_threads"]

# Chromosomes kept in the accessible genome: the canonical set minus
# cnvkit.drop_chr.  Everything else (alt/random/unplaced contigs, chrM,
# chrEBV) is dropped, so no downstream step ever sees those bins.
DROP_CHR   = set(config["cnvkit"].get("drop_chr") or [])
CANON_CHR  = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
KEEP_CHR   = [c for c in CANON_CHR if c not in DROP_CHR]

GENELIST   = config["cnvkit"].get("genelist") or {}
if "amp_log2" in GENELIST or "del_log2" in GENELIST:
    sys.stderr.write(
        "WARNING: cnvkit.genelist.amp_log2/del_log2 are deprecated and "
        "ignored; gene amp/del calling now uses absolute copy number "
        "(cnvkit.genelist.amp_offset/del_offset).\n"
    )
AMP_OFFSET = GENELIST.get("amp_offset", 1)
DEL_OFFSET = GENELIST.get("del_offset", 1)
MIN_PROBES = GENELIST.get("min_probes", 5)

# ─── Helper functions ─────────────────────────────────────────────────────────

def bam_path(sample):
    """Absolute path to the recalibrated BAM for *sample*."""
    if BAM_DIR:
        return f"{BAM_DIR}/{sample}/{sample}.recal.bam"
    return f"{PROJECT}/preprocessing/recalibrated/{sample}/{sample}.recal.bam"

def bai_path(sample):
    return bam_path(sample) + ".bai"

def get_bam(wildcards):
    return bam_path(wildcards.sample)

def get_bai(wildcards):
    return bai_path(wildcards.sample)

def get_reference(wildcards):
    """Return the CNN reference matching the requested analysis mode."""
    row = COMPARISON_MAP.get((wildcards.mode, wildcards.sample))
    if row is None:
        raise ValueError(f"Unregistered comparison: {wildcards.mode}/{wildcards.sample}")
    if row["reference"] == "human_reference":
        return f"{OUTDIR}/references/flat_reference.cnn"
    if row["reference"] == NORMAL:
        return f"{OUTDIR}/references/normal_reference.cnn"
    return f"{OUTDIR}/references/{row['reference']}.cnn"

def resolve_ploidy_purity(mode, sample):
    """Resolve (ploidy, purity) for one (mode, sample) comparison."""
    return resolve_call(config, mode, sample)

def call_opts(wildcards):
    """--ploidy/--purity flags for `cnvkit call`, resolved per comparison."""
    ploidy, purity = resolve_ploidy_purity(wildcards.mode, wildcards.sample)
    opts = [f"--ploidy {ploidy}"]
    if purity:
        opts.append(f"--purity {purity}")
    return " ".join(opts)

def get_ploidy(wildcards):
    """Resolved ploidy for a (mode, sample) comparison, for genelist."""
    ploidy, _ = resolve_ploidy_purity(wildcards.mode, wildcards.sample)
    return ploidy

# ─── Wildcard constraints ─────────────────────────────────────────────────────
wildcard_constraints:
    mode      = "|".join(re.escape(m) for m in ["vs_reference", "vs_normal"] + PAIRWISE_MODES),
    sample    = "|".join(re.escape(s) for s in ALL_SAMPLES),
    reference = "|".join(re.escape(r) for r in REFERENCE_SAMPLES) or "(?!)",

# =============================================================================
# Target rule — collect all final outputs
# =============================================================================
rule all:
    input:
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.genemetrics.tsv",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.amplified.tsv",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.deleted.tsv",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.call.cns",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.scatter.png",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.diagram.pdf",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.amplified.genes.txt",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.deleted.genes.txt",
            zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
        ),
        *([f"{OUTDIR}/vs_reference/heatmap.pdf"] if VS_REF_SAMPLES else []),
        *([f"{OUTDIR}/vs_normal/heatmap.pdf"] if VS_NOR_SAMPLES else []),
        *REPORT_OUTPUTS,


# =============================================================================
# STEP 0 — Accessible regions
# Sequencing-accessible regions of the genome (assembly minus long N runs).
# CNVkit's 'wgs' binning method requires these; derived once from the FASTA.
# Any BED files in cnvkit.access_exclude (e.g. the ENCODE blacklist) are
# additionally subtracted via -x, so blacklisted regions are never binned.
# =============================================================================
rule access:
    input:
        fasta   = FASTA,
        exclude = ACCESS_EXCLUDE,
    output:
        bed = f"{OUTDIR}/bins/access.bed",
    params:
        cnvkit  = CNVKIT,
        min_gap = ACCESS_MIN_GAP,
        exclude = lambda w, input: " ".join(f"-x {bed}" for bed in input.exclude),
        keep    = ",".join(KEEP_CHR),
    threads: 1
    log:
        f"{OUTDIR}/logs/access.log"
    shell:
        """
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})"

        {params.cnvkit} access {input.fasta} \
            --min-gap-size {params.min_gap} \
            {params.exclude} \
            -o {output.bed}.all \
        2>&1 | tee {log}

        # Keep only the wanted chromosomes; drops alt/random/unplaced contigs
        # and anything listed in cnvkit.drop_chr.
        awk -v keep="{params.keep}" \
            'BEGIN {{ n = split(keep, k, ","); for (i = 1; i <= n; i++) ok[k[i]] = 1 }}
             ok[$1]' \
            {output.bed}.all > {output.bed}
        rm -f {output.bed}.all

        echo "kept: $(cut -f1 {output.bed} | sort -u | paste -sd' ' -)" | tee -a {log}
        """


# =============================================================================
# STEP 1 — Autobin
# Generate genome-wide WGS target and antitarget BED files.
# Done once using the normal BAM; the resulting bins are shared by all samples.
# Method 'wgs' requires the accessible-regions BED, passed via -g.
# =============================================================================
rule autobin:
    input:
        bam     = bam_path(NORMAL),
        bai     = bai_path(NORMAL),
        access  = f"{OUTDIR}/bins/access.bed",
        fasta   = FASTA,
        refflat = REFFLAT,
    output:
        target     = f"{OUTDIR}/bins/cnvkit_targets.bed",
        antitarget = f"{OUTDIR}/bins/cnvkit_antitargets.bed",
    params:
        cnvkit = CNVKIT,
        method = METHOD,
    threads: 1
    log:
        f"{OUTDIR}/logs/autobin.log"
    shell:
        """
        mkdir -p "$(dirname {output.target})" "$(dirname {log})"

        {params.cnvkit} autobin {input.bam} \
            --method {params.method} \
            --fasta {input.fasta} \
            --access {input.access} \
            --annotate {input.refflat} \
            --target-output-bed {output.target} \
            --antitarget-output-bed {output.antitarget} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 2 — Coverage
# Compute per-sample read depth for target and antitarget bins in a single job
# (splitting into two jobs would double the BAM I/O cost).
# Parallelised across samples on the cluster.
# =============================================================================
rule coverage:
    input:
        bam        = get_bam,
        bai        = get_bai,
        target     = f"{OUTDIR}/bins/cnvkit_targets.bed",
        antitarget = f"{OUTDIR}/bins/cnvkit_antitargets.bed",
    output:
        target_cov     = f"{OUTDIR}/coverage/{{sample}}.targetcoverage.cnn",
        antitarget_cov = f"{OUTDIR}/coverage/{{sample}}.antitargetcoverage.cnn",
    params:
        cnvkit = CNVKIT,
        mapq   = MAPQ,
    threads: COV_THREADS
    log:
        f"{OUTDIR}/logs/coverage/{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname {log})"

        {{
            {params.cnvkit} coverage {input.bam} {input.target} \
                -p {threads} \
                -q {params.mapq} \
                -o {output.target_cov}

            {params.cnvkit} coverage {input.bam} {input.antitarget} \
                -p {threads} \
                -q {params.mapq} \
                -o {output.antitarget_cov}
        }} 2>&1 | tee {log}
        """


# =============================================================================
# STEP 3a — Flat reference
# Builds a theoretical (log2 = 0) reference from the bin BED files alone.
# GC content and mappability bias factors are computed from the FASTA.
# This is the baseline for the "vs_reference" analysis mode.
# =============================================================================
rule build_flat_reference:
    input:
        target     = f"{OUTDIR}/bins/cnvkit_targets.bed",
        antitarget = f"{OUTDIR}/bins/cnvkit_antitargets.bed",
        fasta      = FASTA,
    output:
        ref = f"{OUTDIR}/references/flat_reference.cnn",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/build_flat_reference.log"
    shell:
        """
        {params.cnvkit} reference \
            -t {input.target} \
            -a {input.antitarget} \
            --fasta {input.fasta} \
            -o {output.ref} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 3b — Normal reference
# Builds a matched-normal reference from Bulk_sensitive coverage files.
# All Bulk_sensitive_1..7 aliases share the same biological sample (same
# FASTQ); a single reference from the canonical Bulk_sensitive BAM is correct.
# This is the baseline for the "vs_normal" analysis mode.
# =============================================================================
rule build_normal_reference:
    input:
        target_cov     = f"{OUTDIR}/coverage/{NORMAL}.targetcoverage.cnn",
        antitarget_cov = f"{OUTDIR}/coverage/{NORMAL}.antitargetcoverage.cnn",
        fasta          = FASTA,
    output:
        ref = f"{OUTDIR}/references/normal_reference.cnn",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/build_normal_reference.log"
    shell:
        """
        {params.cnvkit} reference \
            {input.target_cov} {input.antitarget_cov} \
            --fasta {input.fasta} \
            -o {output.ref} \
        2>&1 | tee {log}
        """


rule build_sample_reference:
    input:
        target_cov     = f"{OUTDIR}/coverage/{{reference}}.targetcoverage.cnn",
        antitarget_cov = f"{OUTDIR}/coverage/{{reference}}.antitargetcoverage.cnn",
        fasta          = FASTA,
    output:
        ref = f"{OUTDIR}/references/{{reference}}.cnn",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/build_sample_reference.{{reference}}.log"
    shell:
        """
        {params.cnvkit} reference \
            {input.target_cov} {input.antitarget_cov} \
            --fasta {input.fasta} \
            -o {output.ref} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 4 — Fix
# Subtract the reference log2 values, apply GC/edge/repeat-mask bias
# corrections, and produce the per-bin ratio file (.cnr).
# Parallelised over (mode, sample) combinations.
# =============================================================================
rule fix:
    input:
        target_cov     = f"{OUTDIR}/coverage/{{sample}}.targetcoverage.cnn",
        antitarget_cov = f"{OUTDIR}/coverage/{{sample}}.antitargetcoverage.cnn",
        reference      = get_reference,
    output:
        cnr = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/fix.{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname {output.cnr})" "$(dirname {log})"

        {params.cnvkit} fix \
            {input.target_cov} {input.antitarget_cov} {input.reference} \
            -o {output.cnr} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 5 — Segment
# Circular binary segmentation (CBS by default) to partition the genome into
# copy-number segments.
# =============================================================================
rule segment:
    input:
        cnr = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
    output:
        cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cns",
    params:
        cnvkit = CNVKIT,
        method = SEG_METH,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/segment.{{sample}}.log"
    shell:
        """
        {params.cnvkit} segment \
            {input.cnr} \
            --method {params.method} \
            -o {output.cns} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 6 — Call
# Convert segment log2 ratios to integer absolute copy-number calls.
# --ploidy / --purity are resolved per (mode, sample) comparison (see
# resolve_ploidy_purity / call_opts) from cnvkit.ploidy/purity,
# cnvkit.mode_defaults and cnvkit.comparison_overrides in the config; CNVkit
# does not estimate them and neither does this pipeline. Runs for every
# comparison listed in rule all (both vs_reference and vs_normal).
# =============================================================================
rule call:
    input:
        cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cns",
    output:
        call_cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.call.cns",
    params:
        cnvkit = CNVKIT,
        opts   = call_opts,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/call.{{sample}}.log"
    shell:
        """
        {params.cnvkit} call \
            {input.cns} \
            {params.opts} \
            -o {output.call_cns} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 7 — Genemetrics
# Per-gene copy-number statistics — mean log2 ratio, p-value, etc.
# Requires gene labels in the bins (provided by autobin --annotate refflat).
# Segments come from the ploidy/purity-aware .call.cns (not the plain .cns),
# so every gene row also carries the absolute copy-number (cn) column that
# rule genelist thresholds on. -t 0 -m 1 disable genemetrics' own log2/probe
# pre-filtering, since that filtering is now done downstream on cn.
# =============================================================================
rule genemetrics:
    input:
        cnr      = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
        call_cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.call.cns",
    output:
        tsv = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.genemetrics.tsv",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/genemetrics.{{sample}}.log"
    shell:
        """
        {params.cnvkit} genemetrics \
            {input.cnr} \
            -s {input.call_cns} \
            -t 0 -m 1 \
            -o {output.tsv} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 7b — Amplified / deleted gene lists
# Split the genemetrics table (which now carries an absolute copy-number `cn`
# column, from the .call.cns segments used in rule genemetrics) into an
# amplified and a deleted table, plus bare gene-symbol lists for pasting into
# enrichment tools. Cutoffs are absolute copy number, derived from this
# comparison's resolved ploidy (cnvkit.ploidy / mode_defaults /
# comparison_overrides) and cnvkit.genelist.amp_offset/del_offset — see
# scripts/gene_calls.py. Sorted by |cn - ploidy|, strongest first.
# =============================================================================
rule genelist:
    input:
        tsv    = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.genemetrics.tsv",
        script = os.path.join(workflow.basedir, "scripts", "gene_calls.py"),
    output:
        amp       = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.amplified.tsv",
        dele      = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.deleted.tsv",
        amp_genes = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.amplified.genes.txt",
        del_genes = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.deleted.genes.txt",
    params:
        ploidy     = get_ploidy,
        amp_offset = AMP_OFFSET,
        del_offset = DEL_OFFSET,
        min_probes = MIN_PROBES,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/genelist.{{sample}}.log"
    shell:
        """
        python3 {input.script:q} \
            --genemetrics {input.tsv} \
            --ploidy {params.ploidy} \
            --amp-offset {params.amp_offset} \
            --del-offset {params.del_offset} \
            --min-probes {params.min_probes} \
            --amp-out {output.amp} \
            --del-out {output.dele} \
            --amp-genes-out {output.amp_genes} \
            --del-genes-out {output.del_genes} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 8 — Scatter plot
# Genome-wide scatter of per-bin log2 ratios with segments overlaid.
# =============================================================================
rule scatter:
    input:
        cnr = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
        cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cns",
    output:
        png = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.scatter.png",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/scatter.{{sample}}.log"
    shell:
        """
        {params.cnvkit} scatter \
            {input.cnr} \
            -s {input.cns} \
            -o {output.png} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 9 — Chromosome arm diagram
# Ideogram-style view of amplifications and deletions per chromosome arm.
# =============================================================================
rule diagram:
    input:
        cnr = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
        cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cns",
    output:
        pdf = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.diagram.pdf",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/diagram.{{sample}}.log"
    shell:
        """
        {params.cnvkit} diagram \
            {input.cnr} \
            -s {input.cns} \
            -o {output.pdf} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 10 — Multi-sample heatmap (aggregate, runs after all samples finish)
# Uses segmented (.cns) files for a clean cross-sample CNV overview.
# =============================================================================
rule heatmap_vs_reference:
    input:
        cns = expand(
            f"{OUTDIR}/vs_reference/{{sample}}/{{sample}}.cns",
            sample=VS_REF_SAMPLES,
        ),
    output:
        pdf = f"{OUTDIR}/vs_reference/heatmap.pdf",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/heatmap_vs_reference.log"
    shell:
        """
        {params.cnvkit} heatmap \
            {input.cns} \
            -o {output.pdf} \
        2>&1 | tee {log}
        """

rule heatmap_vs_normal:
    input:
        cns = expand(
            f"{OUTDIR}/vs_normal/{{sample}}/{{sample}}.cns",
            sample=VS_NOR_SAMPLES,
        ),
    output:
        pdf = f"{OUTDIR}/vs_normal/heatmap.pdf",
    params:
        cnvkit = CNVKIT,
    threads: 1
    log:
        f"{OUTDIR}/logs/heatmap_vs_normal.log"
    shell:
        """
        {params.cnvkit} heatmap \
            {input.cns} \
            -o {output.pdf} \
        2>&1 | tee {log}
        """


if MANIFEST is not None:
    rule manifest:
        input:
            config = workflow.configfiles,
            wfcfg  = os.path.join(workflow.basedir, "scripts", "workflow_config.py"),
        output:
            json = f"{OUTDIR}/summary/manifest.json",
        params:
            manifest = lambda w: json.dumps({**MANIFEST, "expected_outputs": REPORT_OUTPUTS}, indent=2),
        threads: 1
        log:
            f"{OUTDIR}/logs/manifest.log"
        shell:
            """
            mkdir -p "$(dirname "{output.json}")" "$(dirname "{log}")"
            printf '%s\n' {params.manifest:q} > "{output.json}"
            """

    rule summary:
        input:
            cnr = expand(
                f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
                zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
            ),
            call_cns = expand(
                f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.call.cns",
                zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
            ),
            genemetrics = expand(
                f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.genemetrics.tsv",
                zip, mode=COMPARISON_MODES, sample=COMPARISON_SAMPLES,
            ),
            manifest = f"{OUTDIR}/summary/manifest.json",
            script   = os.path.join(workflow.basedir, "scripts", "cnv_summary.R"),
        output:
            REPORT_OUTPUTS,
        conda:
            "envs/cnv_summary.yaml",
        threads: 1
        log:
            f"{OUTDIR}/logs/summary.log"
        shell:
            """
            mkdir -p "$(dirname "{log}")"
            Rscript --vanilla {input.script:q} {input.manifest:q} 2>&1 | tee "{log}"
            """