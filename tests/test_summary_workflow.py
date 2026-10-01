import json
import os
import re
import subprocess
import sys
from collections import Counter

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))
sys.path.insert(0, os.path.join(REPO, "tests"))
from workflow_config import comparisons, summary_manifest, summary_outputs
from test_cnv_summary import RSCRIPT, _config, _html_trace_names, _make_inputs

SNAKEFILE = os.path.join(REPO, "Snakefile")


def _write_config(tmp_path, cfg):
    cfgfile = tmp_path / "config.json"
    cfgfile.write_text(json.dumps(cfg))
    return str(cfgfile)


def _r_env():
    env = dict(os.environ)
    rbin = os.path.dirname(RSCRIPT[0])
    env["PATH"] = rbin + os.pathsep + env["PATH"]
    env["R_ENVIRON_USER"] = os.devnull
    env["R_PROFILE_USER"] = os.devnull
    return env


def _run_snakemake(cfgfile, targets, cwd, env):
    cmd = ["snakemake", "--cores", "1", "-s", SNAKEFILE] + list(targets) + \
          ["--allowed-rules", "manifest", "summary", "--configfile", cfgfile]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, env=env)


def _ran_rule(combined, name):
    return bool(re.search(rf"^(?:local)?rule {name}:", combined, re.M))


def _job_counts(stdout):
    return Counter(re.findall(r"^(?:local)?rule (\w+):", stdout, re.M))


@pytest.mark.skipif(not os.path.exists(RSCRIPT[0]), reason="R interpreter not found")
def test_summary_workflow_execution(tmp_path):
    cfg = _config(tmp_path)
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    expected = summary_outputs(manifest)
    _make_inputs(rows)
    cfgfile = _write_config(tmp_path, cfg)
    env = _r_env()

    res = _run_snakemake(cfgfile, expected, tmp_path, env)
    assert res.returncode == 0, (res.stdout + res.stderr)[-4000:]
    for p in expected:
        assert os.path.exists(p), p
    out = cfg["outdir"]
    assert os.path.exists(os.path.join(out, "logs", "summary.log"))

    res2 = _run_snakemake(cfgfile, expected, tmp_path, env)
    assert res2.returncode == 0, (res2.stdout + res2.stderr)[-4000:]
    assert "Nothing to be done" in res2.stdout + res2.stderr

    html = os.path.join(out, "summary", "cn_scatters", "s1.html")
    os.remove(html)
    res3 = _run_snakemake(cfgfile, expected, tmp_path, env)
    assert res3.returncode == 0, (res3.stdout + res3.stderr)[-4000:]
    assert os.path.exists(html)
    assert _ran_rule(res3.stdout + res3.stderr, "summary")

    cfg["analysis"]["thresholds"]["log2_high"] = 0.5
    _write_config(tmp_path, cfg)
    res4 = _run_snakemake(cfgfile, expected, tmp_path, env)
    assert res4.returncode == 0, (res4.stdout + res4.stderr)[-4000:]
    with open(os.path.join(out, "summary", "manifest.json")) as fh:
        m2 = json.load(fh)
    assert m2["thresholds"]["log2_high"] == 0.5
    assert _ran_rule(res4.stdout + res4.stderr, "summary")

    s1 = _html_trace_names(os.path.join(out, "summary", "cn_scatters", "s1.html"))
    assert any(n in ("gain in both", "loss in both", "divergent", "neutral in both") for n in s1)
    assert os.path.exists(os.path.join(out, "summary", "cn_scatters", "s1_abs.html"))
    s2 = _html_trace_names(os.path.join(out, "summary", "cn_scatters", "s2.html"))
    assert s2
    assert not os.path.exists(os.path.join(out, "summary", "cn_scatters", "s2_abs.html"))
    for c in ("c1_c2", "c1_c3"):
        assert os.path.exists(os.path.join(out, "summary", "conserved", c, c + "_gene_cn_all.csv"))


def test_default_config_dag(tmp_path):
    import yaml
    with open(os.path.join(REPO, "config.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    proj = tmp_path / "proj"
    out = tmp_path / "out"
    cfg["project_dir"] = str(proj)
    cfg["outdir"] = str(out)
    cfg["fasta"] = str(tmp_path / "genome.fa")
    cfg["refflat"] = str(out / "reference" / "hg38.refFlat.txt")
    cfg["cnvkit"]["access_exclude"] = [str(out / "reference" / "hg38-blacklist.v2.bed")]
    for p in (cfg["fasta"], cfg["refflat"], cfg["cnvkit"]["access_exclude"][0]):
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, "w").close()
    for s in cfg["all_samples"]:
        d = proj / "preprocessing" / "recalibrated" / s
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{s}.recal.bam").touch()
        (d / f"{s}.recal.bam.bai").touch()
    cfgfile = _write_config(tmp_path, cfg)
    res = subprocess.run(
        ["snakemake", "-n", "-s", SNAKEFILE, "--configfile", cfgfile],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert res.returncode == 0, res.stderr[-4000:]
    rows = comparisons(cfg)
    manifest = summary_manifest(cfg, rows)
    outputs = summary_outputs(manifest)
    counts = _job_counts(res.stdout)
    assert counts["fix"] == len(rows)
    assert counts["summary"] == 1
    assert counts["manifest"] == 1
    assert counts["coverage"] == len(cfg["all_samples"])
    refs = sorted({r["reference"] for r in rows
                   if r["reference"] not in ("human_reference", cfg["normal_sample"])})
    assert counts["build_sample_reference"] == len(refs)
    for p in outputs:
        assert p in res.stdout, p
