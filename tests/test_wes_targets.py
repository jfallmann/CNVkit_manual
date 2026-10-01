import json
import os
import shlex
import subprocess
import time

import pytest

from test_workflow_config import SNAKEFILE, base_config


@pytest.fixture
def binning(tmp_path):
    cfg = base_config(
        project_dir=str(tmp_path / "project"),
        outdir=str(tmp_path / "out"),
        fasta=str(tmp_path / "genome.fa"),
        refflat=str(tmp_path / "refFlat.txt"),
        bam_dir=str(tmp_path / "bams"),
    )
    bam = tmp_path / "bams" / "N" / "N.recal.bam"
    access = tmp_path / "out" / "bins" / "access.bed"
    for path in (bam, bam.with_suffix(".bam.bai"), access,
                 tmp_path / "genome.fa", tmp_path / "refFlat.txt"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    cfg["cnvkit"]["method"] = "hybrid"
    return cfg


def dry_run(tmp_path, cfg):
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(cfg))
    target = str(tmp_path / "out" / "bins" / "cnvkit_targets.bed")
    result = subprocess.run(
        ["snakemake", "-n", "-p", "-s", SNAKEFILE, target,
         "--configfile", str(config_path)],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
    )
    return result.returncode, result.stdout + result.stderr


def autobin_args(output):
    command = next(line.strip() for line in output.splitlines()
                   if "cnvkit.py autobin " in line)
    return shlex.split(command)


def test_hybrid_targets_expanded_quoted_and_tracked(tmp_path, binning, monkeypatch):
    bed = tmp_path / "capture targets.bed"
    bed.write_text("chr1\t0\t1000\n")
    monkeypatch.setenv("WES_CAPTURE_BED", str(bed))
    binning["cnvkit"]["target_bed"] = "${WES_CAPTURE_BED}"
    code, output = dry_run(tmp_path, binning)
    assert code == 0, output
    args = autobin_args(output)
    assert args[args.index("--method") + 1] == "hybrid"
    assert args[args.index("--targets") + 1] == str(bed)
    inputs = next(line for line in output.splitlines() if line.strip().startswith("input:"))
    assert str(bed) in inputs

    timestamp = time.time() + 10
    for name in ("cnvkit_targets.bed", "cnvkit_antitargets.bed"):
        path = tmp_path / "out" / "bins" / name
        path.write_text("chr1\t0\t1000\n")
        os.utime(path, (timestamp, timestamp))
    code, output = dry_run(tmp_path, binning)
    assert code == 0, output
    assert "Nothing to be done" in output
    os.utime(bed, (timestamp + 10, timestamp + 10))
    code, output = dry_run(tmp_path, binning)
    assert code == 0, output
    assert "rule autobin:" in output
    assert "Updated input files" in output


@pytest.mark.parametrize("target_bed", [None, "", "   ", []])
def test_hybrid_requires_target_bed(tmp_path, binning, target_bed):
    if target_bed is not None:
        binning["cnvkit"]["target_bed"] = target_bed
    code, output = dry_run(tmp_path, binning)
    assert code != 0
    assert "requires cnvkit.target_bed" in output


def test_hybrid_missing_bed_fails_dag(tmp_path, binning):
    bed = tmp_path / "missing.bed"
    binning["cnvkit"]["target_bed"] = str(bed)
    code, output = dry_run(tmp_path, binning)
    assert code != 0
    assert "MissingInputException" in output
    assert str(bed) in output


@pytest.mark.parametrize("configured_bed", [False, True])
def test_wgs_does_not_use_capture_targets(tmp_path, binning, configured_bed):
    binning["cnvkit"]["method"] = "wgs"
    if configured_bed:
        binning["cnvkit"]["target_bed"] = str(tmp_path / "unused.bed")
    code, output = dry_run(tmp_path, binning)
    assert code == 0, output
    args = autobin_args(output)
    assert args[args.index("--method") + 1] == "wgs"
    assert "--targets" not in args
