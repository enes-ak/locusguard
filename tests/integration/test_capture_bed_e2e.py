"""End-to-end: capture-bed aware WES and Panel pipelines."""
from __future__ import annotations

import json
import subprocess
import sys

import pytest

from locusguard.api import Annotator
from locusguard.config import load_config

_SMN1_CFG_YAML = """
schema_version: "1.0"
locusguard_compat: ">=0.1.0,<1.0.0"
locus: {id: SMN1, name: SMN1, gene_family: SMN, paralogs: [SMN2]}
reference: grch38
coordinates:
  primary: {chrom: chr5, start: 12000, end: 15000}
  paralogs:
    SMN2: {chrom: chr5, start: 2000, end: 5000}
psvs:
  - name: "c.840C>T"
    chrom: chr5
    pos: 14001
    alleles: {SMN1: C, SMN2: T}
evidence_weights:
  default:
    psv_match: 0.40
    haplotype_consistency: 0.25
    mapq_pattern: 0.10
    softclip: 0.05
    unique_kmer: 0.10
    coverage_ratio: 0.10
  profile_overrides:
    ont_panel:
      disable: [coverage_ratio]
    ont_wes:
      disable: [coverage_ratio]
confidence_thresholds: {resolved: 0.80, probable: 0.50}
"""


@pytest.mark.integration
def test_panel_with_partial_coverage_reports_missing_psv(
    tmp_path, smn_like_bam, smn_like_fasta, mini_vcf, sma_panel_bed_missing_psv,
):
    """Panel mode + bed missing the PSV position → missing PSV in summary,
    manifest warning."""
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=smn_like_fasta,
        data_type="panel",
        capture_bed=sma_panel_bed_missing_psv,
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    manifest_path = tmp_path / "out.manifest.json"
    annotator.annotate_vcf(
        bam=smn_like_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
        manifest_path=manifest_path,
    )
    summary = json.loads(summary_path.read_text())
    assert summary["loci"]["SMN1"]["psv_coverage"]["missing"] == ["c.840C>T"]
    assert summary["loci"]["SMN1"]["psv_coverage"]["covered"] == []
    assert summary["loci"]["SMN1"]["psv_coverage"]["fraction_covered"] == 0.0

    manifest = json.loads(manifest_path.read_text())
    assert any(
        "SMN1: PSV c.840C>T" in w and "not in capture bed" in w
        for w in manifest["warnings"]
    )


@pytest.mark.integration
def test_wes_with_full_coverage_reports_all_covered(
    tmp_path, smn_like_bam, smn_like_fasta, mini_vcf, sma_panel_bed_full_coverage,
):
    """WES mode + bed covering the PSV → fraction_covered=1.0, no missing-PSV warnings."""
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=smn_like_fasta,
        data_type="wes",
        capture_bed=sma_panel_bed_full_coverage,
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    manifest_path = tmp_path / "out.manifest.json"
    annotator.annotate_vcf(
        bam=smn_like_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
        manifest_path=manifest_path,
    )
    summary = json.loads(summary_path.read_text())
    assert summary["loci"]["SMN1"]["psv_coverage"]["fraction_covered"] == 1.0
    assert summary["loci"]["SMN1"]["psv_coverage"]["missing"] == []

    manifest = json.loads(manifest_path.read_text())
    assert not any(
        "not in capture bed" in w for w in manifest["warnings"]
    )


@pytest.mark.integration
def test_wes_without_capture_bed_emits_informational_warning(
    tmp_path, smn_like_bam, smn_like_fasta, mini_vcf,
):
    """--data-type wes without --capture-bed → pipeline runs, summary lacks
    psv_coverage field, manifest has 'cannot be verified' notice."""
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=smn_like_fasta,
        data_type="wes",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    manifest_path = tmp_path / "out.manifest.json"
    annotator.annotate_vcf(
        bam=smn_like_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
        manifest_path=manifest_path,
    )
    summary = json.loads(summary_path.read_text())
    assert "psv_coverage" not in summary["loci"]["SMN1"]

    manifest = json.loads(manifest_path.read_text())
    assert any(
        "without --capture-bed" in w and "cannot be verified" in w
        for w in manifest["warnings"]
    )


@pytest.mark.integration
def test_wgs_rejects_capture_bed_via_cli(
    tmp_path, smn_like_bam, smn_like_fasta, mini_vcf, sma_panel_bed_full_coverage,
):
    """CLI invocation with --data-type wgs + --capture-bed → non-zero exit code."""
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    # The CLI subprocess path uses the conda env's python (has numpy).
    result = subprocess.run(
        [
            sys.executable, "-m", "locusguard.cli.main",
            "annotate",
            "--bam", str(smn_like_bam),
            "--vcf", str(mini_vcf),
            "--output", str(tmp_path / "out.vcf.gz"),
            "--reference-fasta", str(smn_like_fasta),
            "--data-type", "wgs",
            "--capture-bed", str(sma_panel_bed_full_coverage),
            "--config", str(cfg_path),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    stderr_stdout = (result.stderr + result.stdout).lower()
    assert "capture-bed" in stderr_stdout
    assert "wes or panel" in stderr_stdout
