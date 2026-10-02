"""Validated comparison registry and explicit report targets for the workflow."""
import math
import os
import re
from pathlib import Path


DEFAULT_ONCOGENES = (
    "KRAS EGFR MYC MYCL MYCN MET ERBB2 BRAF PIK3CA CDKN2A CDKN2B TP53 RB1 "
    "PTEN MDM2 MDM4 CCND1 CCNE1 CDK4 CDK6 SOX2 TERT NKX2-1 ALK ROS1 RET "
    "FGFR1 KEAP1 STK11 SMARCA4 YAP1 NFE2L2 AURKA TITF1"
).split()
DEFAULT_THRESHOLDS = dict(log2_high=0.4, log2_low=0.2, probe_high=100,
                          probe_med=30, min_size_kb=100, gap_merge=1e6,
                          gene_log2=0.3, gene_weight=0.5)
RESERVED_REFERENCES = ("flat_reference", "normal_reference", "human_reference")


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError("Sample and analysis identifiers must use letters, digits, _, . or -")
    return value


def expand_path(value):
    value = os.path.expanduser(os.path.expandvars(value))
    if "$" in value:
        raise ValueError("Unresolved environment variable in a configured path; export PROJECT_DIR/OUTPUT_DIR or use literal paths")
    return os.path.abspath(value)


def resolve_call(config, mode, sample):
    cnv = config["cnvkit"]
    values = dict(ploidy=cnv.get("ploidy", 2), purity=cnv.get("purity"))
    values.update((cnv.get("mode_defaults") or {}).get(mode, {}))
    override = (cnv.get("comparison_overrides") or {}).get(sample) or {}
    values.update({k: override[k] for k in ("ploidy", "purity") if k in override})
    values.update(override.get(mode, {}))
    if mode.startswith("pairwise_"):
        spec = (config.get("pairwise_comparisons") or {})[mode[len("pairwise_"):]]
        values.update({k: spec[k] for k in ("ploidy", "purity") if k in spec})
    ploidy, purity = values["ploidy"], values["purity"]
    if isinstance(ploidy, bool) or not isinstance(ploidy, (int, float)) or not math.isfinite(ploidy) or ploidy < 1 or int(ploidy) != ploidy:
        raise ValueError(f"{mode}/{sample}: ploidy must be a positive integer")
    if purity is not None and (isinstance(purity, bool) or not isinstance(purity, (int, float)) or not math.isfinite(purity) or not 0 < purity <= 1):
        raise ValueError(f"{mode}/{sample}: purity must be null or in (0, 1]")
    return int(ploidy), purity


def comparisons(config):
    samples = config["all_samples"]
    if not samples or len(set(samples)) != len(samples):
        raise ValueError("all_samples must be nonempty and unique")
    for sample in samples:
        identifier(sample)
        if sample in RESERVED_REFERENCES:
            raise ValueError(f"Sample name {sample!r} collides with a reserved reference name")
    if config["normal_sample"] not in samples:
        raise ValueError("normal_sample must be in all_samples")
    rows = []

    def add(mode, sample, reference):
        if sample not in samples or (reference != "human_reference" and reference not in samples):
            raise ValueError(f"Unknown sample/reference in {mode}")
        ploidy, purity = resolve_call(config, mode, sample)
        prefix = f'{config["outdir"]}/{mode}/{sample}/{sample}'
        rows.append(dict(id=f"{mode}__{sample}", mode=mode, sample=sample,
                         reference=reference, comparison_type="vs_reference" if mode == "vs_reference" else "pairwise",
                         ploidy=ploidy, purity=purity, cnr=prefix + ".cnr",
                         cns=prefix + ".call.cns", genemetrics=prefix + ".genemetrics.tsv"))

    for mode, key in (("vs_reference", "vs_ref_samples"), ("vs_normal", "vs_normal_samples")):
        selected = config.get(key) or []
        if len(set(selected)) != len(selected):
            raise ValueError(f"Duplicate samples in {key}")
        for sample in selected:
            add(mode, sample, "human_reference" if mode == "vs_reference" else config["normal_sample"])
    for name, spec in (config.get("pairwise_comparisons") or {}).items():
        identifier(name)
        if spec["sample"] == spec["reference"]:
            raise ValueError(f"Self-comparison is not allowed: {name}")
        add("pairwise_" + name, spec["sample"], spec["reference"])
    if not rows:
        raise ValueError("Configure at least one comparison")
    if len({r["id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate comparison IDs (check pairwise names and sample names for '__' collisions)")
    return rows


def summary_manifest(config, rows):
    settings = config.get("analysis") or {}
    if not settings.get("enabled", False):
        return None
    registry = {row["id"]: row for row in rows}
    thresholds = {**DEFAULT_THRESHOLDS, **(settings.get("thresholds") or {})}
    for key, value in thresholds.items():
        if key not in DEFAULT_THRESHOLDS or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError(f"Invalid analysis threshold: {key}")
    if thresholds["log2_low"] > thresholds["log2_high"] or thresholds["probe_med"] > thresholds["probe_high"]:
        raise ValueError("Analysis low/medium thresholds must not exceed high thresholds")
    resistant = settings.get("resistance_samples") or []
    if not set(resistant) <= set(config["all_samples"]):
        raise ValueError("Unknown analysis.resistance_samples entry")
    temporal, conserved, scatters = [], [], []
    for clone, spec in (settings.get("temporal_clones") or {}).items():
        identifier(clone)
        timeline = []
        for n, entry in enumerate(spec.get("timeline") or [], 1):
            key = entry["comparison"]
            if key not in registry:
                raise ValueError(f"Unknown comparison in {clone}: {key}")
            label = entry.get("label", key)
            if not isinstance(label, str) or not label.strip():
                raise ValueError(f"{clone}: timeline labels must be nonempty strings")
            timeline.append(dict(timepoint=n, label=label, dir=key))
        if len(timeline) < 2 or len({x["label"] for x in timeline}) != len(timeline) or len({x["dir"] for x in timeline}) != len(timeline):
            raise ValueError(f"{clone}: series needs >=2 unique comparisons and labels")
        temporal.append(dict(clone=clone, out_subdir="temporal/" + clone, timeline=timeline))
    clones = {x["clone"]: x for x in temporal}
    for name, spec in (settings.get("conserved_comparisons") or {}).items():
        identifier(name)
        pair = spec.get("clones") or []
        if len(pair) != 2 or len(set(pair)) != 2 or not set(pair) <= set(clones):
            raise ValueError(f"{name}: specify two configured temporal clones")
        for clone in pair:
            if registry[clones[clone]["timeline"][0]["dir"]]["mode"] != "vs_reference":
                raise ValueError(f"{name}: conserved analysis requires a flat-reference first comparison")
        conserved.append(dict(clones=pair, out_subdir="conserved/" + name))
    for name, spec in (settings.get("scatter_comparisons") or {}).items():
        identifier(name)
        pair = spec.get("comparisons") or []
        if len(pair) != 2 or len(set(pair)) != 2 or not set(pair) <= set(registry):
            raise ValueError(f"{name}: specify two distinct registered comparisons")
        title, xlab, ylab = (spec.get("title", name), spec.get("xlab", pair[0]),
                             spec.get("ylab", pair[1]))
        for lbl in (title, xlab, ylab):
            if not isinstance(lbl, str) or not lbl.strip():
                raise ValueError(f"{name}: scatter title/xlab/ylab must be nonempty strings")
        scatter = dict(out_prefix=name, title=title, xlab=xlab, ylab=ylab,
                       timeline=[dict(timepoint=i, label=key, dir=key) for i, key in enumerate(pair, 1)])
        gm = spec.get("genemetrics_comparisons")
        if gm is not None:
            if len(gm) != 2 or not set(gm) <= set(registry):
                raise ValueError(f"{name}: invalid genemetrics_comparisons")
            scatter["gm"] = dict(x=gm[0], y=gm[1])
        scatters.append(scatter)
    plots = dict(top_n=40, heatmap_max=60, cn_cap=10, heatmap_cap=6, label_top=25, log2_cap=3)
    plots.update(settings.get("plots") or {})
    for key, value in plots.items():
        if key not in {"top_n", "heatmap_max", "cn_cap", "heatmap_cap", "label_top", "log2_cap"} or isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"Invalid positive integer plotting setting: {key}")
    oncogenes = settings.get("oncogenes")
    if oncogenes is None:
        oncogenes = DEFAULT_ONCOGENES
    return dict(comparisons=rows, out_dir=f'{config["outdir"]}/summary', thresholds=thresholds,
                oncogenes=oncogenes, plots=plots,
                resistance_samples=resistant, temporal_clones=temporal,
                conserved_comparisons=conserved, scatter_comparisons=scatters)


def summary_outputs(manifest):
    if manifest is None:
        return []
    outputs = ["all_cnv_segments_scored.csv", "segments_scored_vs_reference.csv",
               "segments_scored_pairwise.csv", "cnv_segments_high_moderate.csv",
               "all_genemetrics_combined.csv", "cnv_regions_vs_reference.csv",
               "cnv_regions_pairwise.csv", "cnv_regions_shared_all_comparisons.csv",
               "gene_log2_matrix.tsv", "resistant_group_genes.tsv", "gene_log2_matrix_cleaned.tsv"]

    def plot(prefix, html=False):
        outputs.extend(prefix + ext for ext in ([".pdf", ".png", ".html"] if html else [".pdf", ".png"]))

    outputs.extend("by_comparison/" + x["id"] + ".csv" for x in manifest["comparisons"])
    plot("manhattan_cnv_all_comparisons")
    for spec in manifest["temporal_clones"]:
        prefix = spec["out_subdir"] + "/" + spec["clone"].lower()
        outputs.append(prefix + "_gene_cn_drift.csv")
        for suffix in ("manhattan_log2_relative", "manhattan_cn_absolute", "temporal_heatmap_genomic",
                       "temporal_heatmap_clustered", "trajectories", "cn_state_alluvial"):
            plot(prefix + "_" + suffix)
        plot("Genomic_drift_" + spec["clone"])
    for spec in manifest["conserved_comparisons"]:
        prefix = spec["out_subdir"] + "/" + "_".join(spec["clones"]).lower()
        outputs.extend(prefix + suffix for suffix in ("_gene_cn_all.csv", "_conserved_cnv.csv"))
        plot(prefix + "_cn_scatter", html=True)
        plot(prefix + "_conserved_heatmap")
    for spec in manifest["scatter_comparisons"]:
        prefix = "cn_scatters/" + spec["out_prefix"]
        plot(prefix, html=True)
        if "gm" in spec:
            plot(prefix + "_abs", html=True)
    paths = [str(Path(manifest["out_dir"]) / p) for p in outputs]
    if len(set(paths)) != len(paths):
        raise ValueError("Duplicate summary output paths (check scatter/conserved/temporal names)")
    return paths
