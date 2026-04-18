"""End-to-end test: CLI invocation on synthetic fixtures produces a valid annotated VCF."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from locusguard.io.vcf import VcfReader


@pytest.mark.integration
def test_annotate_cli_end_to_end(tmp_path, smn_like_bam, smn_like_fasta, mini_vcf):
    # Write two config YAMLs matching the smn_like_fasta layout
    smn1_yaml = tmp_path / "smn1.yaml"
    smn1_yaml.write_text("""
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
    pos: 14000
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
""")

    smn2_yaml = tmp_path / "smn2.yaml"
    smn2_yaml.write_text("""
schema_version: "1.0"
locusguard_compat: ">=0.1.0,<1.0.0"
locus: {id: SMN2, name: SMN2, gene_family: SMN, paralogs: [SMN1]}
reference: grch38
coordinates:
  primary: {chrom: chr5, start: 2000, end: 5000}
  paralogs:
    SMN1: {chrom: chr5, start: 12000, end: 15000}
psvs:
  - name: "c.840C>T"
    chrom: chr5
    pos: 4000
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
""")

    output = tmp_path / "out.vcf.gz"

    env_args = {
        "LOCUSGUARD_GRCH38_FASTA": str(smn_like_fasta),
    }

    result = subprocess.run(
        [
            "locusguard", "annotate",
            "--bam", str(smn_like_bam),
            "--vcf", str(mini_vcf),
            "--output", str(output),
            "--reference", "grch38",
            "--tech", "ont",
            "--data-type", "wgs",
            "--config", str(smn1_yaml),
            "--config", str(smn2_yaml),
        ],
        capture_output=True,
        text=True,
        env={**_os_environ(), **env_args},
    )
    assert result.returncode == 0, (
        f"CLI exited non-zero. stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    )
    assert "done:" in result.stdout

    # Validate output artifacts
    assert output.exists()
    summary = tmp_path / "out.summary.json"
    manifest = tmp_path / "out.manifest.json"
    assert summary.exists()
    assert manifest.exists()

    reader = VcfReader(output)
    variants = list(reader.iter_variants())
    assert len(variants) == 2

    # The SMN1-region variant (POS=13950 in 1-based VCF) should have TRUE_LOCUS=SMN1
    smn1_variant = next(v for v in variants if v.POS == 13950)
    assert smn1_variant.INFO.get("TRUE_LOCUS") == "SMN1"
    assert smn1_variant.INFO.get("LOCUS_STATUS") == "RESOLVED"
    assert float(smn1_variant.INFO.get("LOCUS_CONF")) >= 0.80

    # SMN2-region variant (POS=3950)
    smn2_variant = next(v for v in variants if v.POS == 3950)
    assert smn2_variant.INFO.get("TRUE_LOCUS") == "SMN2"

    # Summary JSON sanity
    summary_doc = json.loads(summary.read_text())
    assert summary_doc["loci"]["SMN1"]["status"] == "RESOLVED"
    assert summary_doc["loci"]["SMN2"]["status"] == "RESOLVED"


def _os_environ() -> dict[str, str]:
    import os
    return dict(os.environ)
