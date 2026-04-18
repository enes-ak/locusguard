import pytest
import yaml

from locusguard.config import load_config


MINIMAL_YAML = """
schema_version: "1.0"
locusguard_compat: ">=0.1.0,<1.0.0"
locus:
  id: SMN1
  name: "Survival of Motor Neuron 1"
  gene_family: SMN
  paralogs: [SMN2]
  clinical_relevance: SMA
reference: grch38
coordinates:
  primary: {chrom: chr5, start: 70924941, end: 70953012}
  paralogs:
    SMN2: {chrom: chr5, start: 69345350, end: 69373421}
psvs:
  - name: "c.840C>T"
    chrom: chr5
    pos: 70951946
    alleles: {SMN1: C, SMN2: T}
    exon: 7
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
"""


def test_load_config_reads_yaml(tmp_path):
    path = tmp_path / "SMN1.yaml"
    path.write_text(MINIMAL_YAML)
    cfg = load_config(path)
    assert cfg.locus.id == "SMN1"
    assert cfg.psvs[0].pos == 70951946


def test_load_config_raises_on_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "does_not_exist.yaml")


def test_load_config_raises_on_malformed_yaml(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("not: valid: yaml: [")
    with pytest.raises(yaml.YAMLError):
        load_config(path)


def test_load_config_raises_on_schema_violation(tmp_path):
    doc = yaml.safe_load(MINIMAL_YAML)
    del doc["psvs"]
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(doc))
    from pydantic import ValidationError
    with pytest.raises(ValidationError):
        load_config(path)
