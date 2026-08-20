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

configfile: "config.yaml"

# ─── Derived constants ────────────────────────────────────────────────────────
PROJECT  = config["project_dir"]
OUTDIR   = config["outdir"]
FASTA    = config["fasta"]
REFFLAT  = config["refflat"]
CONDA    = config["conda_env"]
NORMAL   = config["normal_sample"]

ALL_SAMPLES    = config["all_samples"]
VS_REF_SAMPLES = config["vs_ref_samples"]
VS_NOR_SAMPLES = config["vs_normal_samples"]

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

# Ploidy / purity are never estimated here: pass them through only if the
# config sets them, otherwise leave CNVkit's own defaults in place.
PLOIDY = config["cnvkit"].get("ploidy")
PURITY = config["cnvkit"].get("purity")
CALL_OPTS = " ".join(
    o for o in (
        f"--ploidy {PLOIDY}" if PLOIDY else "",
        f"--purity {PURITY}" if PURITY else "",
    ) if o
)

# Chromosomes kept in the accessible genome: the canonical set minus
# cnvkit.drop_chr.  Everything else (alt/random/unplaced contigs, chrM,
# chrEBV) is dropped, so no downstream step ever sees those bins.
DROP_CHR   = set(config["cnvkit"].get("drop_chr") or [])
CANON_CHR  = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
KEEP_CHR   = [c for c in CANON_CHR if c not in DROP_CHR]

GENELIST   = config["cnvkit"].get("genelist") or {}
AMP_LOG2   = GENELIST.get("amp_log2", 0.585)
DEL_LOG2   = GENELIST.get("del_log2", -0.585)
MIN_PROBES = GENELIST.get("min_probes", 5)

# ─── Helper functions ─────────────────────────────────────────────────────────

def bam_path(sample):
    """Absolute path to the recalibrated BAM for *sample*."""
    return f"{PROJECT}/preprocessing/recalibrated/{sample}/{sample}.recal.bam"

def bai_path(sample):
    return bam_path(sample) + ".bai"

def get_bam(wildcards):
    return bam_path(wildcards.sample)

def get_bai(wildcards):
    return bai_path(wildcards.sample)

def get_reference(wildcards):
    """Return the CNN reference matching the requested analysis mode."""
    if wildcards.mode == "vs_reference":
        return f"{OUTDIR}/references/flat_reference.cnn"
    elif wildcards.mode == "vs_normal":
        return f"{OUTDIR}/references/normal_reference.cnn"
    raise ValueError(f"Unrecognised mode wildcard: {wildcards.mode}")

# ─── Wildcard constraints ─────────────────────────────────────────────────────
wildcard_constraints:
    mode   = "vs_reference|vs_normal",
    sample = "|".join(ALL_SAMPLES),

# =============================================================================
# Target rule — collect all final outputs
# =============================================================================
rule all:
    input:
        # ── vs_reference ──────────────────────────────────────────────────────
        expand(
            f"{OUTDIR}/vs_reference/{{sample}}/{{sample}}.genemetrics.tsv",
            sample=VS_REF_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_reference/{{sample}}/{{sample}}.amplified.tsv",
            sample=VS_REF_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_reference/{{sample}}/{{sample}}.deleted.tsv",
            sample=VS_REF_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_reference/{{sample}}/{{sample}}.call.cns",
            sample=VS_REF_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_reference/{{sample}}/{{sample}}.scatter.png",
            sample=VS_REF_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_reference/{{sample}}/{{sample}}.diagram.pdf",
            sample=VS_REF_SAMPLES,
        ),
        f"{OUTDIR}/vs_reference/heatmap.pdf",

        # ── vs_normal ─────────────────────────────────────────────────────────
        expand(
            f"{OUTDIR}/vs_normal/{{sample}}/{{sample}}.genemetrics.tsv",
            sample=VS_NOR_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_normal/{{sample}}/{{sample}}.amplified.tsv",
            sample=VS_NOR_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_normal/{{sample}}/{{sample}}.deleted.tsv",
            sample=VS_NOR_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_normal/{{sample}}/{{sample}}.call.cns",
            sample=VS_NOR_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_normal/{{sample}}/{{sample}}.scatter.png",
            sample=VS_NOR_SAMPLES,
        ),
        expand(
            f"{OUTDIR}/vs_normal/{{sample}}/{{sample}}.diagram.pdf",
            sample=VS_NOR_SAMPLES,
        ),
        f"{OUTDIR}/vs_normal/heatmap.pdf",


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
        bam    = bam_path(NORMAL),
        bai    = bai_path(NORMAL),
        access = f"{OUTDIR}/bins/access.bed",
    output:
        target     = f"{OUTDIR}/bins/cnvkit_targets.bed",
        antitarget = f"{OUTDIR}/bins/cnvkit_antitargets.bed",
    params:
        cnvkit  = CNVKIT,
        method  = METHOD,
        fasta   = FASTA,
        refflat = REFFLAT,
    threads: 1
    log:
        f"{OUTDIR}/logs/autobin.log"
    shell:
        """
        mkdir -p "$(dirname {output.target})" "$(dirname {log})"

        {params.cnvkit} autobin {input.bam} \
            --method {params.method} \
            --fasta {params.fasta} \
            --access {input.access} \
            --annotate {params.refflat} \
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
    output:
        ref = f"{OUTDIR}/references/flat_reference.cnn",
    params:
        cnvkit = CNVKIT,
        fasta  = FASTA,
    threads: 1
    log:
        f"{OUTDIR}/logs/build_flat_reference.log"
    shell:
        """
        {params.cnvkit} reference \
            -t {input.target} \
            -a {input.antitarget} \
            --fasta {params.fasta} \
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
    output:
        ref = f"{OUTDIR}/references/normal_reference.cnn",
    params:
        cnvkit = CNVKIT,
        fasta  = FASTA,
    threads: 1
    log:
        f"{OUTDIR}/logs/build_normal_reference.log"
    shell:
        """
        {params.cnvkit} reference \
            {input.target_cov} {input.antitarget_cov} \
            --fasta {params.fasta} \
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
# --ploidy / --purity are passed only when set in the config; CNVkit does not
# estimate them and neither does this pipeline.
# =============================================================================
rule call:
    input:
        cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cns",
    output:
        call_cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.call.cns",
    params:
        cnvkit = CNVKIT,
        opts   = CALL_OPTS,
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
# =============================================================================
rule genemetrics:
    input:
        cnr = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cnr",
        cns = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.cns",
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
            -s {input.cns} \
            -o {output.tsv} \
        2>&1 | tee {log}
        """


# =============================================================================
# STEP 7b — Amplified / deleted gene lists
# Split the genemetrics table at fixed log2 cutoffs (cnvkit.genelist in the
# config) into an amplified and a deleted table, plus bare gene-symbol lists
# for pasting into enrichment tools.  Sorted by |log2|, strongest first.
# =============================================================================
rule genelist:
    input:
        tsv = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.genemetrics.tsv",
    output:
        amp       = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.amplified.tsv",
        dele      = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.deleted.tsv",
        amp_genes = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.amplified.genes.txt",
        del_genes = f"{OUTDIR}/{{mode}}/{{sample}}/{{sample}}.deleted.genes.txt",
    params:
        amp_log2   = AMP_LOG2,
        del_log2   = DEL_LOG2,
        min_probes = MIN_PROBES,
    threads: 1
    log:
        f"{OUTDIR}/logs/{{mode}}/genelist.{{sample}}.log"
    shell:
        r"""
        # Column positions are read from the header, not hardcoded.
        # Rows are prefixed with |log2| for sorting, then the key is cut off.
        split_genes() {{
            awk -F'\t' -v OFS='\t' -v want="$1" \
                -v amp={params.amp_log2} -v del={params.del_log2} \
                -v minp={params.min_probes} '
                NR == 1 {{
                    for (i = 1; i <= NF; i++) col[$i] = i
                    L = col["log2"]; P = col["probes"]
                    print "key", $0
                    next
                }}
                $P < minp {{ next }}
                (want == "amp" && $L >= amp) || (want == "del" && $L <= del) {{
                    print ($L < 0 ? -$L : $L), $0
                }}' {input.tsv} \
            | {{ IFS= read -r hdr; echo "$hdr"; sort -t"$(printf '\t')" -k1,1gr; }} \
            | cut -f2-
        }}

        split_genes amp > {output.amp}
        split_genes del > {output.dele}

        # bare gene symbols (column 1 of the genemetrics table), header dropped
        cut -f1 {output.amp}  | tail -n +2 > {output.amp_genes}
        cut -f1 {output.dele} | tail -n +2 > {output.del_genes}

        {{
            echo "amplified (log2 >= {params.amp_log2}): $(wc -l < {output.amp_genes}) genes"
            echo "deleted   (log2 <= {params.del_log2}): $(wc -l < {output.del_genes}) genes"
        }} | tee {log}
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