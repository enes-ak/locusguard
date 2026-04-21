import json

from locusguard.api import Annotator
from locusguard.config import load_config


def test_annotator_runs_end_to_end_on_fixtures(
    smn_like_bam, smn_like_fasta, mini_vcf, tmp_path, monkeypatch
):
    # Use the SMN1 fixture config (not the bundled, since bundled coords don't
    # match our synthetic FASTA). We define it inline for the test.
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text("""
schema_version: "1.0"
locusguard_compat: ">=0.1.0,<1.0.0"
locus:
  id: SMN1
  name: SMN1
  gene_family: SMN
  paralogs: [SMN2]
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
confidence_thresholds:
  resolved: 0.80
  probable: 0.50
""")
    smn1 = load_config(cfg_path)

    cfg_path2 = tmp_path / "smn2.yaml"
    cfg_path2.write_text("""
schema_version: "1.0"
locusguard_compat: ">=0.1.0,<1.0.0"
locus:
  id: SMN2
  name: SMN2
  gene_family: SMN
  paralogs: [SMN1]
reference: grch38
coordinates:
  primary: {chrom: chr5, start: 2000, end: 5000}
  paralogs:
    SMN1: {chrom: chr5, start: 12000, end: 15000}
psvs:
  - name: "c.840C>T"
    chrom: chr5
    pos: 4001
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
confidence_thresholds:
  resolved: 0.80
  probable: 0.50
""")
    smn2 = load_config(cfg_path2)

    annotator = Annotator(
        configs=[smn1, smn2],
        reference_fasta=smn_like_fasta,
        tech="ont",
        data_type="wgs",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    annotator.annotate_vcf(
        bam=smn_like_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
        summary_path=tmp_path / "summary.json",
        manifest_path=tmp_path / "manifest.json",
    )

    assert out_vcf.exists()
    assert (tmp_path / "summary.json").exists()
    assert (tmp_path / "manifest.json").exists()

    summary = json.loads((tmp_path / "summary.json").read_text())
    assert "SMN1" in summary["loci"]
    assert "SMN2" in summary["loci"]
    assert summary["loci"]["SMN1"]["status"] == "RESOLVED"


def test_annotation_result_has_cn_by_locus(
    smn_like_bam, smn_like_fasta, mini_vcf, tmp_path,
):
    from locusguard.api import Annotator
    from locusguard.config import load_config

    # Minimal SMN1 config with a control region (will be insufficient_depth
    # since the fixture FASTA is only 20kb — but graceful path is valid).
    cfg_path = tmp_path / "smn1.yaml"
    cfg_path.write_text("""
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
  - {name: ctrl_A, chrom: chr5, start: 17000, end: 18000}
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
""")
    smn1 = load_config(cfg_path)

    annotator = Annotator(
        configs=[smn1],
        reference_fasta=smn_like_fasta,
        tech="ont",
        data_type="wgs",
    )
    out_vcf = tmp_path / "out.vcf.gz"
    result = annotator.annotate_vcf(
        bam=smn_like_bam,
        vcf_in=mini_vcf,
        vcf_out=out_vcf,
    )
    assert "SMN1" in result.cn_by_locus
    cn = result.cn_by_locus["SMN1"]
    assert cn.locus_id == "SMN1"
    # Control region is empty → insufficient_depth or similar graceful fail.
    assert cn.status in ("ok", "insufficient_depth", "no_control_regions")
