"""End-to-end data, tables and public-document consistency."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import re
import shutil

import pytest

from hea_md.campaign import load_campaign
from hea_md.cli import main

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("field,value", [("strain_unit", "percent"), ("stress_unit", "MPa"),
                                         ("stress_convention", "tension positive")])
def test_manifest_units_fail_closed(tmp_path, field, value):
    shutil.copytree(ROOT / "configs", tmp_path / "configs")
    p = tmp_path / "configs/campaign.json"
    manifest = json.loads(p.read_text(encoding="utf-8"))
    manifest[field] = value
    p.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="unit or convention"):
        load_campaign(tmp_path)


@pytest.mark.parametrize("file", ["README.md"])
def test_readme_primary_numbers(file):
    rows = json.loads((ROOT / "results/analysis.json").read_text(encoding="utf-8"))["summary"]
    text = (ROOT / file).read_text(encoding="utf-8")
    for row in rows:
        if row["case_id"] in {"eq", "nb28", "ta12"}:
            for key in ("slope_GPa", "proof_stress_GPa", "peak_common_GPa"):
                assert f"{row[key]:.2f}" in text
        if row["case_id"] == "nb28":
            for key in ("slope_change_percent", "peak_change_percent"):
                assert f"{row[key]:.2f}%" in text


def test_markdown_local_links_exist():
    for file in ROOT.rglob("*.md"):
        if any(part.startswith(".") for part in file.relative_to(ROOT).parts):
            continue
        for href in re.findall(r"\]\(([^)\s]+)\)", file.read_text(encoding="utf-8")):
            if "://" in href or href.startswith(("#", "mailto:")):
                continue
            assert (file.parent / href.split("#", 1)[0]).exists(), (file, href)


def test_cli_verify(capsys):
    assert main(["--root", str(ROOT), "verify"]) == 0
    assert "PASS" in capsys.readouterr().out


def test_cli_output_isolated_and_checks_match(tmp_path, capsys):
    assert main(["--root", str(ROOT), "analyze", "--output", str(tmp_path / "result")]) == 0
    assert main(["--root", str(ROOT), "verify", "--results", str(tmp_path / "result/analysis.json")]) == 0
    assert "PASS" in capsys.readouterr().out
