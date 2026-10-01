import json
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from workflow_config import (comparisons, resolve_call, summary_manifest,
                             summary_outputs)

SNAKEFILE = os.path.join(REPO, "Snakefile")


def base_config(**overrides):
    cfg = {
        "project_dir": "/tmp/cnvproj",
        "outdir": "/tmp/cnvout",
        "fasta": "/tmp/genome.fa",
        "refflat": "/tmp/refFlat.txt",
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
        "analysis": {"enabled": False},
    }
    cfg.update(overrides)
    return cfg


def test_registry():
    rows = comparisons(base_config())
    by_id = {r["id"]: r for r in rows}
    assert set(by_id) == {
        "vs_reference__N", "vs_reference__T1", "vs_reference__T2",
        "vs_normal__T1", "vs_normal__T2", "pairwise_t1_vs_t2__T1",
    }
    assert by_id["vs_reference__T1"]["reference"] == "human_reference"
    assert by_id["vs_normal__T1"]["reference"] == "N"
    assert by_id["pairwise_t1_vs_t2__T1"]["reference"] == "T2"
    assert by_id["pairwise_t1_vs_t2__T1"]["comparison_type"] == "pairwise"
    assert by_id["vs_reference__T1"]["cnr"].endswith("vs_reference/T1/T1.cnr")


def test_ploidy_purity_precedence():
    cfg = base_config()
    cfg["cnvkit"]["comparison_overrides"] = {
        "T1": {"ploidy": 4, "vs_normal": {"purity": 0.7}},
    }
    assert resolve_call(cfg, "vs_reference", "T1") == (4, None)
    assert resolve_call(cfg, "vs_normal", "T1") == (4, 0.7)
    assert resolve_call(cfg, "vs_reference", "T2") == (2, None)
    assert resolve_call(cfg, "vs_normal", "T2") == (3, None)
    cfg["pairwise_comparisons"]["t1_vs_t2"]["ploidy"] = 5
    assert resolve_call(cfg, "pairwise_t1_vs_t2", "T1") == (5, None)


def test_unknowns():
    cfg = base_config()
    cfg["analysis"] = {
        "enabled": True,
        "temporal_clones": {"C1": {"timeline": [
            {"comparison": "vs_reference__T1", "label": "a"},
            {"comparison": "ghost__T1", "label": "b"},
        ]}},
    }
    with pytest.raises(ValueError):
        summary_manifest(cfg, comparisons(cfg))
    cfg2 = base_config()
    cfg2["pairwise_comparisons"] = {"x": {"sample": "T1", "reference": "GHOST"}}
    with pytest.raises(ValueError):
        comparisons(cfg2)


def test_disabled():
    manifest = summary_manifest(base_config(), comparisons(base_config()))
    assert manifest is None
    assert summary_outputs(None) == []


def test_null_mappings():
    cfg = base_config()
    cfg["cnvkit"]["mode_defaults"] = None
    cfg["cnvkit"]["comparison_overrides"] = None
    cfg["analysis"] = {
        "enabled": True,
        "resistance_samples": None,
        "temporal_clones": None,
        "conserved_comparisons": None,
        "scatter_comparisons": None,
    }
    manifest = summary_manifest(cfg, comparisons(cfg))
    assert manifest is not None
    assert manifest["temporal_clones"] == []
    assert manifest["conserved_comparisons"] == []
    assert manifest["scatter_comparisons"] == []


def test_invalid_ploidy_purity():
    for bad in (0, 2.5, True, "2"):
        cfg = base_config()
        cfg["cnvkit"]["ploidy"] = bad
        with pytest.raises(ValueError):
            comparisons(cfg)
    for bad in (0, -0.1, 1.5, True):
        cfg = base_config()
        cfg["cnvkit"]["purity"] = bad
        with pytest.raises(ValueError):
            comparisons(cfg)


def test_invalid_mode_overrides():
    cfg = base_config()
    cfg["cnvkit"]["comparison_overrides"] = {"T1": {"vs_normal": {"purity": 1.5}}}
    with pytest.raises(ValueError):
        comparisons(cfg)
    cfg = base_config()
    cfg["cnvkit"]["comparison_overrides"] = {"T1": {"vs_normal": {"ploidy": 2.5}}}
    with pytest.raises(ValueError):
        comparisons(cfg)


def test_reserved_reference_sample_rejected():
    cfg = base_config()
    cfg["all_samples"] = ["flat_reference", "T1", "T2"]
    cfg["normal_sample"] = "flat_reference"
    with pytest.raises(ValueError):
        comparisons(cfg)


def test_duplicate_comparison_ids_rejected():
    cfg = base_config()
    cfg["all_samples"] = ["N", "T1", "T2", "y__T1"]
    cfg["pairwise_comparisons"] = {
        "x__y": {"sample": "T1", "reference": "T2"},
        "x": {"sample": "y__T1", "reference": "T2"},
    }
    with pytest.raises(ValueError):
        comparisons(cfg)


def test_duplicate_summary_outputs_rejected():
    cfg = base_config()
    cfg["analysis"] = {
        "enabled": True,
        "scatter_comparisons": {
            "foo": {"comparisons": ["vs_reference__T1", "vs_reference__T2"],
                    "genemetrics_comparisons": ["vs_reference__T1", "vs_reference__T2"]},
            "foo_abs": {"comparisons": ["vs_normal__T1", "vs_normal__T2"]},
        },
    }
    with pytest.raises(ValueError):
        summary_outputs(summary_manifest(cfg, comparisons(cfg)))


def test_empty_label_rejected():
    cfg = base_config()
    cfg["analysis"] = {
        "enabled": True,
        "temporal_clones": {"C1": {"timeline": [
            {"comparison": "vs_reference__T1", "label": ""},
            {"comparison": "vs_normal__T1", "label": "b"},
        ]}},
    }
    with pytest.raises(ValueError):
        summary_manifest(cfg, comparisons(cfg))


def _write_config(tmp_path, cfg):
    cfgfile = tmp_path / "config.json"
    cfgfile.write_text(json.dumps(cfg))
    return str(cfgfile)


def test_synthetic_dag(tmp_path):
    bamdir = tmp_path / "bams"
    fasta = tmp_path / "genome.fa"
    refflat = tmp_path / "refFlat.txt"
    for s in ("N", "T1", "T2"):
        d = bamdir / s
        d.mkdir(parents=True)
        (d / f"{s}.recal.bam").touch()
        (d / f"{s}.recal.bam.bai").touch()
    fasta.touch()
    refflat.touch()
    cfg = base_config(
        project_dir=str(tmp_path / "proj"),
        outdir=str(tmp_path / "out"),
        fasta=str(fasta),
        refflat=str(refflat),
        bam_dir=str(bamdir),
    )
    cfg["analysis"] = {
        "enabled": True,
        "temporal_clones": {"C1": {"timeline": [
            {"comparison": "vs_reference__T1", "label": "a"},
            {"comparison": "pairwise_t1_vs_t2__T1", "label": "b"},
        ]}},
    }
    cfgfile = _write_config(tmp_path, cfg)
    res = subprocess.run(
        ["snakemake", "-n", "-s", SNAKEFILE, "--configfile", cfgfile],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr
    assert "build_sample_reference" in res.stdout
    assert "manifest" in res.stdout
    assert "fix" in res.stdout


def test_manifest_execution_and_rerun(tmp_path):
    outdir = tmp_path / "out"
    cfg = base_config(outdir=str(outdir))
    cfg["analysis"] = {
        "enabled": True,
        "temporal_clones": {"C1": {"timeline": [
            {"comparison": "vs_reference__T1", "label": "a"},
            {"comparison": "pairwise_t1_vs_t2__T1", "label": "b"},
        ]}},
    }
    cfgfile = _write_config(tmp_path, cfg)
    target = str(outdir / "summary" / "manifest.json")
    cmd = ["snakemake", "--cores", "1", "-s", SNAKEFILE, target, "--configfile", cfgfile]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    with open(target) as fh:
        manifest = json.load(fh)
    ids = {c["id"] for c in manifest["comparisons"]}
    assert "pairwise_t1_vs_t2__T1" in ids
    assert "expected_outputs" in manifest
    assert manifest["temporal_clones"][0]["timeline"][1]["dir"] == "pairwise_t1_vs_t2__T1"
    res2 = subprocess.run(cmd, capture_output=True, text=True)
    assert res2.returncode == 0, res2.stderr
    assert "Nothing to be done" in res2.stdout + res2.stderr
