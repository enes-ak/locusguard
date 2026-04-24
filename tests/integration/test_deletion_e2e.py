"""End-to-end: homozygous-deletion call on a synthetic SMN1-deleted fixture."""
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
  profile_overrides: {}
confidence_thresholds: {resolved: 0.80, probable: 0.50}
"""


@pytest.mark.integration
def test_smn1_homozygous_deletion_detected(
    tmp_path, smn1_deleted_bam, smn_like_fasta, mini_vcf,
):
    """SMN1 region populated by SMN2-pattern reads → HOMOZYGOUS_DELETION."""
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=smn_like_fasta,
        data_type="wgs",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    annotator.annotate_vcf(
        bam=smn1_deleted_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
    )
    doc = json.loads(summary_path.read_text())
    assert doc["loci"]["SMN1"]["deletion_status"] == "HOMOZYGOUS_DELETION"


@pytest.mark.integration
def test_smn_like_bam_not_deleted(
    tmp_path, smn_like_bam, smn_like_fasta, mini_vcf,
):
    """The standard smn_like fixture has SMN1-pattern reads → PRESENT or
    INDETERMINATE (depending on how many reach the assigner). Must NOT
    be HOMOZYGOUS_DELETION."""
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=smn_like_fasta,
        data_type="wgs",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    annotator.annotate_vcf(
        bam=smn_like_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
    )
    doc = json.loads(summary_path.read_text())
    assert doc["loci"]["SMN1"]["deletion_status"] in ("PRESENT", "INDETERMINATE")
