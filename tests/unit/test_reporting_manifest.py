import json

from locusguard.reporting.manifest import write_manifest


def test_manifest_contains_run_metadata(tmp_path):
    out = tmp_path / "manifest.json"
    write_manifest(
        output_path=out,
        locusguard_version="0.1.0",
        command_line="locusguard annotate --bam X.bam --vcf Y.vcf.gz",
        reference_fasta_path="/refs/grch38.fa",
        reference_fasta_md5="abc123",
        config_hashes={"SMN1": "sha256:aaaa", "SMN2": "sha256:bbbb"},
        data_type="wgs",
        profile_used="ont_wgs",
        runtime_seconds=42.5,
        warnings=["profile 'foo' not defined; falling back to default"],
    )
    doc = json.loads(out.read_text())
    assert doc["locusguard_version"] == "0.1.0"
    assert "locusguard annotate" in doc["command_line"]
    assert doc["profile_used"] == "ont_wgs"
    assert doc["reference_fasta"]["md5"] == "abc123"
    assert doc["config_hashes"]["SMN1"] == "sha256:aaaa"
    assert doc["warnings"] == ["profile 'foo' not defined; falling back to default"]


