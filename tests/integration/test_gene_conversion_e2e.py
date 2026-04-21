"""End-to-end: synthetic gene-conversion BAM triggers GENE_CONVERSION_FLAG."""
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
  - name: PSV3
    chrom: chr5
    pos: 14000
    alleles: {SMN1: T, SMN2: A}
gene_conversion:
  known_hotspots:
    - {name: region_A, chrom: chr5, start: 12900, end: 13100}
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
def test_gene_conversion_flagged(
    tmp_path, gene_conversion_bam, multi_psv_fasta, mini_vcf,
):
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text(_SMN1_CFG_YAML)
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=multi_psv_fasta,
        tech="ont",
        data_type="wgs",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    summary_path = tmp_path / "out.summary.json"
    result = annotator.annotate_vcf(
        bam=gene_conversion_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=summary_path,
    )

    # Gene conversion flag should be raised on SMN1 locus
    summary = json.loads(summary_path.read_text())
    assert summary["loci"]["SMN1"]["gene_conv_flag"] is True

    # And the haplotype cluster should carry the note
    clusters = result.haplotype_clusters_by_locus["SMN1"]
    assert any("gene_conversion_suspected" in c.notes for c in clusters)
