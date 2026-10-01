import csv
import json
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from workflow_config import comparisons, summary_manifest, summary_outputs

_default_rscript = shutil.which("Rscript")
RSCRIPT = shlex.split(os.environ.get(
    "CNV_SUMMARY_RSCRIPT",
    f"{_default_rscript} --vanilla" if _default_rscript else "Rscript --vanilla",
))
SUMMARY_R = os.path.join(REPO, "scripts", "cnv_summary.R")

CHROMS = ["chr1", "chr2", "chr3"]

SEG_BASE = [
    ("chr1", 100000, 200000, "KRAS;ENST00000311936,EGFR", 0.8, 4),
    ("chr1", 300000, 400000, "TP53,RB1", -0.6, 1),
    ("chr2", 100000, 200000, "MYC,CDK4", 0.5, 3),
    ("chr2", 300000, 400000, "PTEN,MDM2", 0.1, 2),
    ("chr3", 100000, 200000, "CCND1,ERBB2", 0.6, 3),
    ("chr3", 300000, 400000, "CDKN2A,SMARCA4", -0.3, 1),
]

GM_BASE = [
    ("KRAS", "chr1", 100000, 200000, 50, 0.8, 1.0, 4),
    ("EGFR", "chr1", 100000, 200000, 50, 0.7, 1.0, 4),
    ("TP53", "chr1", 300000, 400000, 40, -0.6, 1.0, 1),
    ("RB1", "chr1", 300000, 400000, 40, -0.5, 1.0, 1),
    ("MYC", "chr2", 100000, 200000, 45, 0.5, 1.0, 3),
    ("CDK4", "chr2", 100000, 200000, 45, 0.4, 1.0, 3),
    ("PTEN", "chr2", 300000, 400000, 35, 0.1, 1.0, 2),
    ("MDM2", "chr2", 300000, 400000, 35, 0.2, 1.0, 2),
    ("CCND1", "chr3", 100000, 200000, 30, 0.6, 1.0, 3),
    ("ERBB2", "chr3", 100000, 200000, 30, 0.5, 1.0, 3),
    ("CDKN2A", "chr3", 300000, 400000, 25, -0.3, 1.0, 1),
    ("SMARCA4", "chr3", 300000, 400000, 25, -0.2, 1.0, 1),
    ("TERT", "chr1", 500000, 600000, 20, 0.05, 0.4, 2),
    ("ALK", "chr2", 500000, 600000, 20, 0.1, 1.0, 2),
    ("RET", "chr3", 500000, 600000, 20, 0.15, 1.0, 2),
]

EMPTY_IDS = ("vs_reference__N", "vs_normal__T2")


def _config(tmp_path):
    return {
        "project_dir": str(tmp_path / "proj"),
        "outdir": str(tmp_path / "out"),
        "fasta": str(tmp_path / "genome.fa"),
        "refflat": str(tmp_path / "refFlat.txt"),
        "conda_env": "cnvkit",
        "normal_sample": "N",
        "all_samples": ["N", "T1", "T2"],
        "vs_ref_samples": ["N", "T1", "T2"],
        "vs_normal_samples": ["T1", "T2"],
        "cnvkit": {
            "method": "wgs",
            "min_mapq": 20,
            "segment_method": "cbs",
            "coverage_threads": 1,
            "ploidy": 2,
            "purity": None,
            "mode_defaults": {"vs_normal": {"ploidy": 3}},
            "comparison_overrides": {},
            "drop_chr": ["chrY"],
            "genelist": {"amp_offset": 1, "del_offset": 1, "min_probes": 5},
            "access_min_gap": 5000,
            "access_exclude": [],
        },
        "pairwise_comparisons": {"t1_vs_t2": {"sample": "T1", "reference": "T2"}},
        "analysis": {
            "enabled": True,
            "resistance_samples": ["T1"],
            "thresholds": {
                "log2_high": 0.4,
                "log2_low": 0.2,
                "probe_high": 100,
                "probe_med": 30,
                "min_size_kb": 100,
                "gap_merge": 1000000,
                "gene_log2": 0.3,
                "gene_weight": 0.5,
            },
            "temporal_clones": {
                "C1": {"timeline": [
                    {"comparison": "vs_reference__T1", "label": "T1 vs reference"},
                    {"comparison": "vs_normal__T1", "label": "T1 vs normal"},
                ]},
                "C2": {"timeline": [
                    {"comparison": "vs_reference__T2", "label": "T2 vs reference"},
                    {"comparison": "vs_normal__T2", "label": "T2 vs normal"},
                ]},
                "C3": {"timeline": [
                    {"comparison": "vs_reference__N", "label": "N vs reference"},
                    {"comparison": "vs_normal__T2", "label": "T2 vs normal"},
                ]},
            },
            "conserved_comparisons": {"c1_c2": {"clones": ["C1", "C2"]},
                                      "c1_c3": {"clones": ["C1", "C3"]}},
            "scatter_comparisons": {
                "s1": {
                    "comparisons": ["vs_reference__T1", "vs_reference__T2"],
                    "title": "T1 vs T2 called CN",
                    "xlab": "T1 called CN (vs reference)",
                    "ylab": "T2 called CN (vs reference)",
                    "genemetrics_comparisons": ["vs_reference__T1", "vs_reference__T2"],
                },
                "s2": {
                    "comparisons": ["vs_normal__T1", "vs_normal__T2"],
                    "title": "T1 vs T2 vs normal",
                    "xlab": "T1 called CN (vs normal)",
                    "ylab": "T2 called CN (vs normal)",
                },
            },
        },
    }


def _segments_for(cid):
    if cid == "vs_normal__T1":
        return [
            ("chr1", 100000, 200000, "KRAS,EGFR", 0.7, 5),
            ("chr1", 300000, 400000, "TP53,RB1", -0.5, 2),
            ("chr2", 100000, 200000, "MYC,CDK4", 0.4, 4),
            ("chr2", 300000, 400000, "PTEN,MDM2", 0.0, 3),
        ]
    return SEG_BASE


def _write_cnr(path):
    rows = []
    for c in CHROMS:
        for i in range(6):
            rows.append([c, i * 100000 + 1, (i + 1) * 100000, ".",
                         1.0, round((i % 5 - 2) * 0.1, 4), 1.0])
    with open(path, "w") as fh:
        fh.write("chromosome\tstart\tend\tgene\tdepth\tlog2\tweight\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


def _write_cns(path, segments):
    with open(path, "w") as fh:
        fh.write("chromosome\tstart\tend\tgene\tlog2\tcn\tprobes\tweight\tdepth\n")
        for chrom, start, end, genes, log2, cn in segments:
            fh.write("\t".join(map(str, [chrom, start, end, genes, log2, cn,
                                         100, 1.0, 1.0])) + "\n")


def _write_genemetrics(path, rows=None):
    rows = GM_BASE if rows is None else rows
    with open(path, "w") as fh:
        fh.write("gene\tchromosome\tstart\tend\tprobes\tlog2\tweight\tcn\n")
        for g in rows:
            fh.write("\t".join(map(str, g)) + "\n")


def _write_cnr_header(path):
    with open(path, "w") as fh:
        fh.write("chromosome\tstart\tend\tgene\tdepth\tlog2\tweight\n")


def _write_cns_header(path):
    with open(path, "w") as fh:
        fh.write("chromosome\tstart\tend\tgene\tlog2\tcn\tprobes\tweight\tdepth\n")


def _write_genemetrics_header(path):
    with open(path, "w") as fh:
        fh.write("gene\tchromosome\tstart\tend\tprobes\tlog2\tweight\tcn\n")


def _make_inputs(rows, seg_fn=None, gm_fn=None):
    for row in rows:
        for key in ("cnr", "cns", "genemetrics"):
            os.makedirs(os.path.dirname(row[key]), exist_ok=True)
        _write_cnr(row["cnr"])
        segs = [] if row["id"] in EMPTY_IDS else (
            seg_fn(row["id"]) if seg_fn else _segments_for(row["id"]))
        _write_cns(row["cns"], segs)
        _write_genemetrics(row["genemetrics"], None if gm_fn is None else gm_fn(row["id"]))


def _make_empty_inputs(rows):
    for row in rows:
        for key in ("cnr", "cns", "genemetrics"):
            os.makedirs(os.path.dirname(row[key]), exist_ok=True)
        _write_cnr_header(row["cnr"])
        _write_cns_header(row["cns"])
        _write_genemetrics_header(row["genemetrics"])


def _write_manifest(manifest):
    path = os.path.join(manifest["out_dir"], "manifest.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(manifest, fh)
    return path


def _run_rscript(manifest_path):
    return subprocess.run(RSCRIPT + [SUMMARY_R, manifest_path],
                          capture_output=True, text=True)


def _read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_tsv_header(path):
    with open(path, newline="") as fh:
        return next(csv.reader(fh, delimiter="\t"))


def _read_csv_header(path):
    with open(path, newline="") as fh:
        return next(csv.reader(fh))


def _read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _pdf_text(path):
    res = subprocess.run(["pdftotext", path, "-"], capture_output=True, text=True)
    assert res.returncode == 0, path
    return res.stdout


def _html_widget_json(path):
    with open(path, "r", errors="ignore") as fh:
        html = fh.read()
    m = re.search(r'<script type="application/json" data-for="htmlwidget-[^"]+">(.*?)</script>',
                  html, re.S)
    assert m, path
    return json.loads(m.group(1))


def _html_trace_names(path):
    data = _html_widget_json(path)
    return [tr.get("name") for tr in data["x"]["data"]]


def _png_size(path):
    with open(path, "rb") as fh:
        data = fh.read(26)
    assert data[:8] == b"\x89PNG\r\n\x1a\n", path
    return struct.unpack(">I", data[16:20])[0], struct.unpack(">I", data[20:24])[0]


def _pdf_ok(path):
    with open(path, "rb") as fh:
        head = fh.read(5)
    return head == b"%PDF-" and os.path.getsize(path) > 100


def _html_has_plotly(path):
    _html_widget_json(path)
    return True


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_end_to_end_all_output_families(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_inputs(rows)
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode == 0, res.stderr[-4000:]

    expected = manifest["expected_outputs"]
    for p in expected:
        assert os.path.exists(p), p

    pdfs = [p for p in expected if p.endswith(".pdf")]
    pngs = [p for p in expected if p.endswith(".png")]
    htmls = [p for p in expected if p.endswith(".html")]
    for p in pdfs:
        assert _pdf_ok(p), p
    for p in pngs:
        w, h = _png_size(p)
        assert w > 100 and h > 100, p
    for p in htmls:
        assert _html_has_plotly(p), p

    out = manifest["out_dir"]

    seg = _read_csv(os.path.join(out, "all_cnv_segments_scored.csv"))
    assert seg and {"comparison", "id", "mode", "sample", "reference",
                    "comparison_type", "ploidy", "purity", "confidence",
                    "direction"} <= set(seg[0])
    assert any(r["id"] == "vs_normal__T1" and r["ploidy"] == "3" for r in seg)
    assert any(r["id"] == "vs_reference__T1" and r["purity"] in ("", "NA") for r in seg)

    empty_by = _read_csv(os.path.join(out, "by_comparison", "vs_reference__N.csv"))
    assert empty_by == [] and os.path.getsize(
        os.path.join(out, "by_comparison", "vs_reference__N.csv")) > 0

    pair_regions = _read_csv(os.path.join(out, "cnv_regions_pairwise.csv"))
    assert "overlaps_vs_reference" in pair_regions[0]
    assert "origin" not in pair_regions[0]

    shared_all = _read_csv(os.path.join(out, "cnv_regions_shared_all_comparisons.csv"))
    assert shared_all == []

    ids = {r["id"] for r in rows}
    gm_matrix = _read_tsv(os.path.join(out, "gene_log2_matrix.tsv"))
    assert gm_matrix and {"gene", "chromosome"} | ids <= set(gm_matrix[0])
    assert _read_tsv(os.path.join(out, "resistant_group_genes.tsv"))
    assert _read_tsv(os.path.join(out, "gene_log2_matrix_cleaned.tsv"))

    drift_c1 = _read_csv(os.path.join(out, "temporal", "C1", "c1_gene_cn_drift.csv"))
    assert drift_c1 and {"gene", "drift", "is_oncogene"} <= set(drift_c1[0])
    assert "T1 vs normal" in drift_c1[0]
    assert any(r["T1 vs normal"] == "3" for r in drift_c1)

    drift_c3 = _read_csv(os.path.join(out, "temporal", "C3", "c3_gene_cn_drift.csv"))
    assert drift_c3 == []

    cons = _read_csv(os.path.join(out, "conserved", "c1_c2", "c1_c2_gene_cn_all.csv"))
    assert cons and {"gene", "cn_C1", "cn_C2", "conserved_state"} <= set(cons[0])
    assert _read_csv(os.path.join(out, "conserved", "c1_c2", "c1_c2_conserved_cnv.csv"))

    # HTML widget JSON carries actual trace data
    names = _html_trace_names(os.path.join(out, "cn_scatters", "s1.html"))
    assert names and any(n in ("gain in both", "loss in both", "divergent",
                               "neutral in both") for n in names)
    assert _pdf_text(os.path.join(out, "cn_scatters", "s1.pdf")).strip()


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_missing_input_fails(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_inputs(rows)
    os.remove(rows[0]["cnr"])
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode != 0
    assert "missing" in (res.stderr + res.stdout).lower()


def _neutral_segments(ploidy):
    return [
        ("chr1", 100000, 200000, "KRAS,EGFR", 0.0, ploidy),
        ("chr1", 300000, 400000, "TP53,RB1", 0.0, ploidy),
        ("chr2", 100000, 200000, "MYC,CDK4", 0.0, ploidy),
        ("chr2", 300000, 400000, "PTEN,MDM2", 0.0, ploidy),
        ("chr3", 100000, 200000, "CCND1,ERBB2", 0.0, ploidy),
        ("chr3", 300000, 400000, "CDKN2A,SMARCA4", 0.0, ploidy),
    ]


def _segments_ploidy_shift(cid):
    if cid == "vs_reference__T1":
        return [("chr1", 100000, 200000, "KRAS,EGFR", 0.6, 3),
                ("chr1", 300000, 400000, "TP53,RB1", 0.6, 3),
                ("chr2", 100000, 200000, "MYC,CDK4", 0.6, 3),
                ("chr2", 300000, 400000, "PTEN,MDM2", 0.6, 3),
                ("chr3", 100000, 200000, "CCND1,ERBB2", 0.6, 3),
                ("chr3", 300000, 400000, "CDKN2A,SMARCA4", 0.6, 3)]
    if cid in ("vs_reference__T2", "vs_normal__T1"):
        return [("chr1", 100000, 200000, "KRAS,EGFR", -0.6, 2),
                ("chr1", 300000, 400000, "TP53,RB1", -0.6, 2),
                ("chr2", 100000, 200000, "MYC,CDK4", -0.6, 2),
                ("chr2", 300000, 400000, "PTEN,MDM2", -0.6, 2),
                ("chr3", 100000, 200000, "CCND1,ERBB2", -0.6, 2),
                ("chr3", 300000, 400000, "CDKN2A,SMARCA4", -0.6, 2)]
    return _segments_for(cid)


def _gm_ploidy_shift(cid):
    if cid == "vs_reference__T1":
        return [(g[0], g[1], g[2], g[3], g[4], 0.6, 1.0, 3) for g in GM_BASE]
    if cid in ("vs_reference__T2", "vs_normal__T1"):
        return [(g[0], g[1], g[2], g[3], g[4], -0.6, 1.0, 2) for g in GM_BASE]
    return GM_BASE


def _segments_missing_gene(cid):
    if cid == "vs_normal__T1":
        return [s for s in _segments_for(cid) if "KRAS" not in s[3]]
    return _segments_for(cid)


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_all_empty_inputs(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_empty_inputs(rows)
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode == 0, res.stderr[-4000:]

    out = manifest["out_dir"]
    for p in manifest["expected_outputs"]:
        assert os.path.exists(p) and os.path.getsize(p) > 0, p

    ids = {r["id"] for r in rows}
    gm_header = _read_tsv_header(os.path.join(out, "gene_log2_matrix.tsv"))
    assert ids <= set(gm_header)
    assert _read_tsv(os.path.join(out, "gene_log2_matrix.tsv")) == []

    cons_header = _read_csv_header(
        os.path.join(out, "conserved", "c1_c2", "c1_c2_gene_cn_all.csv"))
    assert {"gene", "cn_C1", "cn_C2", "conserved_state"} <= set(cons_header)
    assert _read_csv(os.path.join(out, "conserved", "c1_c2", "c1_c2_conserved_cnv.csv")) == []

    drift_header = _read_csv_header(
        os.path.join(out, "temporal", "C1", "c1_gene_cn_drift.csv"))
    assert {"gene", "drift", "is_oncogene", "T1 vs reference", "T1 vs normal"} <= set(drift_header)

    txt = _pdf_text(os.path.join(out, "cn_scatters", "s1.pdf"))
    assert "No shared genes" in txt


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_one_empty_scatter(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_inputs(rows)
    for row in rows:
        if row["id"] == "vs_reference__T2":
            _write_cns_header(row["cns"])
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode == 0, res.stderr[-4000:]

    out = manifest["out_dir"]
    for ext in (".pdf", ".png", ".html"):
        p = os.path.join(out, "cn_scatters", "s1" + ext)
        assert os.path.exists(p) and os.path.getsize(p) > 0, p
    assert "No shared genes" in _pdf_text(os.path.join(out, "cn_scatters", "s1.pdf"))


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_all_neutral_no_conserved_alterations(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_inputs(rows, seg_fn=lambda cid: _neutral_segments(
        2 if cid.startswith("vs_reference") else 3))
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode == 0, res.stderr[-4000:]

    out = manifest["out_dir"]
    cons_alt = _read_csv(os.path.join(out, "conserved", "c1_c2", "c1_c2_conserved_cnv.csv"))
    assert cons_alt == []
    hm = os.path.join(out, "conserved", "c1_c2", "c1_c2_conserved_heatmap.pdf")
    assert os.path.exists(hm) and os.path.getsize(hm) > 0
    assert "No conserved altered genes" in _pdf_text(hm)


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_missing_cn_column_fails(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_inputs(rows)
    row = rows[0]
    with open(row["cns"], "w") as fh:
        fh.write("chromosome\tstart\tend\tgene\tlog2\tprobes\tweight\tdepth\n")
        for chrom, start, end, genes, log2, cn in _segments_for(row["id"]):
            fh.write("\t".join(map(str, [chrom, start, end, genes, log2, 100, 1.0, 1.0])) + "\n")
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode != 0
    assert "missing required column" in (res.stderr + res.stdout).lower()


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_ploidy_baselines_not_medians(tmp_path):
    cfg = _config(tmp_path)
    cfg["cnvkit"]["comparison_overrides"] = {"T2": {"ploidy": 3}}
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_inputs(rows, seg_fn=_segments_ploidy_shift, gm_fn=_gm_ploidy_shift)
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode == 0, res.stderr[-4000:]

    out = manifest["out_dir"]
    cons = _read_csv(os.path.join(out, "conserved", "c1_c2", "c1_c2_gene_cn_all.csv"))
    assert cons
    states = {r["conserved_state"] for r in cons}
    assert states == {"divergent"}

    cons_names = _html_trace_names(
        os.path.join(out, "conserved", "c1_c2", "c1_c2_cn_scatter.html"))
    assert "divergent" in cons_names and "neutral in both" not in cons_names

    s1_names = _html_trace_names(os.path.join(out, "cn_scatters", "s1.html"))
    assert "divergent" in s1_names and "neutral in both" not in s1_names

    abs_names = _html_trace_names(os.path.join(out, "cn_scatters", "s1_abs.html"))
    assert "divergent" in abs_names and "neutral in both" not in abs_names

    assert "configured ploidy (2, 3)" in _pdf_text(os.path.join(out, "cn_scatters", "s1.pdf"))


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_missing_calls_remain_na(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    manifest["expected_outputs"] = summary_outputs(manifest)
    _make_inputs(rows, seg_fn=_segments_missing_gene)
    manifest_path = _write_manifest(manifest)

    res = _run_rscript(manifest_path)
    assert res.returncode == 0, res.stderr[-4000:]

    out = manifest["out_dir"]
    drift_c1 = _read_csv(os.path.join(out, "temporal", "C1", "c1_gene_cn_drift.csv"))
    kras = [r for r in drift_c1 if r["gene"] == "KRAS"]
    assert kras and kras[0]["T1 vs normal"] in ("", "NA")
