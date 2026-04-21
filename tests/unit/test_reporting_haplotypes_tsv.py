from locusguard.reporting.haplotypes_tsv import write_haplotypes_tsv
from locusguard.types import HaplotypeCluster


def _c(hap_id: str, n_reads: int, pattern: dict[str, str], notes: list[str]) -> HaplotypeCluster:
    return HaplotypeCluster(
        hap_id=hap_id,
        supporting_reads=[f"r{i}" for i in range(n_reads)],
        psv_pattern=pattern,
        assigned_locus="SMN1",
        confidence=0.9,
        notes=notes,
    )


def test_writes_haplotype_rows(tmp_path):
    out = tmp_path / "haplotypes.tsv"
    write_haplotypes_tsv(
        output_path=out,
        clusters_by_locus={
            "SMN1": [
                _c("H1", 15, {"P1": "C", "P2": "A"}, []),
                _c("H2", 1, {"P1": "?", "P2": "?"}, []),
            ],
        },
    )
    lines = out.read_text().splitlines()
    assert lines[0].startswith("hap_id\t")
    assert len(lines) == 3
    # H2 has 1 support read -> low_support=true
    h2_row = next(line for line in lines[1:] if line.startswith("H2\t"))
    assert h2_row.endswith("true") or "\ttrue" in h2_row


def test_notes_are_compacted(tmp_path):
    out = tmp_path / "haplotypes.tsv"
    cluster = _c("H1", 5, {"P1": "C"}, ["gene_conversion_suspected", "hotspot_match:exon_7"])
    write_haplotypes_tsv(output_path=out, clusters_by_locus={"SMN1": [cluster]})
    lines = out.read_text().splitlines()
    assert "gene_conversion_suspected" in lines[1]
    assert "hotspot_match:exon_7" in lines[1]
