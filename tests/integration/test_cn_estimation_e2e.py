"""End-to-end: CN estimation on a synthetic fixture with controlled depths."""
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
control_regions:
  - {name: ctrl_A, chrom: chr5, start: 6000, end: 9000}
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
def test_cn_estimation_on_synthetic_fixture(
    tmp_path, smn_like_bam, smn_like_fasta, mini_vcf,
):
    """The smn_like fixture has 10 reads in SMN1 region (12000-15000) and
    10 reads in SMN2 region (2000-5000). A control region at 6000-9000 has
    ZERO reads, so CN estimation will report insufficient_depth — but the
    pipeline must produce a valid CN block gracefully."""
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=smn_like_fasta,
        tech="ont",
        data_type="wgs",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    result = annotator.annotate_vcf(
        bam=smn_like_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
    )

    # AnnotationResult carries CN
    assert "SMN1" in result.cn_by_locus
    cn = result.cn_by_locus["SMN1"]
    # Control region (6000-9000) is empty → insufficient_depth graceful fail
    assert cn.status == "insufficient_depth"

    # Summary JSON contains the cn_estimate block
    summary = json.loads(summary_path.read_text())
    cn_block = summary["loci"]["SMN1"]["cn_estimate"]
    assert cn_block["status"] == "insufficient_depth"
    assert cn_block["absolute_cn"] is None
