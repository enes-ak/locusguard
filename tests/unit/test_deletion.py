"""Unit tests for deletion status classification."""
from __future__ import annotations

import pytest

from locusguard.deletion import classify_deletion
from locusguard.types import Assignment


def _make_assignments(statuses: list[str]) -> list[Assignment]:
    """Build a list of Assignment stubs with the given status values."""
    return [
        Assignment(
            read_id=f"r{i}",
            assigned_locus="SMN1" if s != "UNASSIGNED" else None,
            confidence=0.9 if s == "RESOLVED" else 0.6 if s == "PROBABLE" else 0.3,
            status=s,  # type: ignore[arg-type]
            evidence_scores=[],
            locus_key="SMN1:abc",
        )
        for i, s in enumerate(statuses)
    ]


def test_empty_assignments_is_indeterminate():
    assert classify_deletion([]) == "INDETERMINATE"


def test_too_few_assignments_is_indeterminate():
    # 9 reads, all RESOLVED — below the 10-read floor
    assert classify_deletion(_make_assignments(["RESOLVED"] * 9)) == "INDETERMINATE"


def test_all_resolved_is_present():
    assert classify_deletion(_make_assignments(["RESOLVED"] * 20)) == "PRESENT"


def test_mixed_confident_is_present():
    statuses = ["RESOLVED"] * 5 + ["PROBABLE"] * 5 + ["AMBIGUOUS"] * 5
    assert classify_deletion(_make_assignments(statuses)) == "PRESENT"


def test_all_ambiguous_with_enough_reads_is_deletion():
    # 15 reads, none confident → SMN1 not detected despite coverage
    assert classify_deletion(_make_assignments(["AMBIGUOUS"] * 15)) == "HOMOZYGOUS_DELETION"


def test_all_unassigned_with_enough_reads_is_deletion():
    assert classify_deletion(_make_assignments(["UNASSIGNED"] * 15)) == "HOMOZYGOUS_DELETION"


def test_mixed_ambiguous_and_unassigned_is_deletion():
    statuses = ["AMBIGUOUS"] * 8 + ["UNASSIGNED"] * 7
    assert classify_deletion(_make_assignments(statuses)) == "HOMOZYGOUS_DELETION"


def test_single_probable_read_flips_to_present():
    # 14 AMBIGUOUS + 1 PROBABLE → enough evidence
    statuses = ["AMBIGUOUS"] * 14 + ["PROBABLE"]
    assert classify_deletion(_make_assignments(statuses)) == "PRESENT"


@pytest.mark.parametrize("count", [10, 11, 50, 500])
def test_threshold_boundary_all_ambiguous(count: int):
    assert classify_deletion(_make_assignments(["AMBIGUOUS"] * count)) == "HOMOZYGOUS_DELETION"
