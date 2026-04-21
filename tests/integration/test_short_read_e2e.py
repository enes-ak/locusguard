"""End-to-end: short-read BAM produces valid output with graceful degradation."""
from __future__ import annotations

import json

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
  - name: PSV1
    chrom: chr5
    pos: 13000
    alleles: {SMN1: C, SMN2: T}
  - name: PSV2
    chrom: chr5
    pos: 13500
    alleles: {SMN1: A, SMN2: G}
evidence_weights:
  default:
    psv_match: 0.40
    haplotype_consistency: 0.25
    mapq_pattern: 0.10
    softclip: 0.05
    unique_kmer: 0.10
    coverage_ratio: 0.10
  profile_overrides:
    short_read_wgs:
      psv_match: 0.50
      cap_confidence: 0.90
      disable: [haplotype_consistency]
confidence_thresholds: {resolved: 0.80, probable: 0.50}
"""


@pytest.mark.integration
def test_short_read_produces_valid_output_with_cap(
    tmp_path, short_read_bam, multi_psv_fasta, mini_vcf,
):
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=multi_psv_fasta,
        tech="short-read",
        data_type="wgs",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    manifest_path = tmp_path / "out.manifest.json"
    annotator.annotate_vcf(
        bam=short_read_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
        manifest_path=manifest_path,
    )

    summary = json.loads(summary_path.read_text())
    assert summary["tech"] == "short-read"
    # Cap confidence at 0.90
    assert summary["loci"]["SMN1"]["mean_confidence"] <= 0.90

    manifest = json.loads(manifest_path.read_text())
    # Degradations should list haplotype_consistency as disabled
    disabled = {d["evidence"] for d in manifest["degradations"]}
    assert "haplotype_consistency" in disabled

    # Phase 2.5: CN estimation should be explicitly not-supported for short-read
    cn_block = summary["loci"]["SMN1"].get("cn_estimate")
    if cn_block is not None:
        assert cn_block["status"] == "not_supported_for_tech"
        assert cn_block["absolute_cn"] is None
