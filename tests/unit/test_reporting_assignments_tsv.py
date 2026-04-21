from locusguard.reporting.assignments_tsv import write_assignments_tsv
from locusguard.types import Assignment, EvidenceScore


def _a(read_id: str, locus: str) -> Assignment:
    return Assignment(
        read_id=read_id,
        assigned_locus=locus,
        confidence=0.9,
        status="RESOLVED",
        evidence_scores=[
            EvidenceScore(source="psv_match", normalized=1.0, raw={}, available=True),
            EvidenceScore(source="softclip", normalized=0.95, raw={}, available=True),
        ],
        locus_key=f"{locus}:abc",
        flags=set(),
        cluster_id="H1",
    )


def test_writes_tsv_with_header_and_rows(tmp_path):
    out = tmp_path / "assignments.tsv"
    write_assignments_tsv(
        output_path=out,
        assignments_by_locus={
            "SMN1": [_a("r1", "SMN1"), _a("r2", "SMN1")],
        },
    )
    lines = out.read_text().splitlines()
    assert lines[0].startswith("read_id\t")
    assert len(lines) == 3  # header + 2 rows
    assert "r1" in lines[1]
    assert "H1" in lines[1]
    assert "psv_match:1.00" in lines[1]


def test_handles_empty_evidence_list(tmp_path):
    out = tmp_path / "empty.tsv"
    empty = Assignment(
        read_id="r0",
        assigned_locus=None,
        confidence=0.0,
        status="UNASSIGNED",
        evidence_scores=[],
        locus_key="X:y",
        flags={"insufficient_coverage"},
    )
    write_assignments_tsv(output_path=out, assignments_by_locus={"X": [empty]})
    lines = out.read_text().splitlines()
    assert lines[1].split("\t")[0] == "r0"
