"""Deletion-status classification from per-read locus assignments.

A locus is ``HOMOZYGOUS_DELETION`` when enough reads were evaluated at that
locus but none showed confident evidence for it (no RESOLVED or PROBABLE
assignments). This is chemistry- and aligner-independent: it relies on the
assigner's evidence-based attribution, not on raw depth ratios.
"""
from __future__ import annotations

from typing import Literal

from locusguard.types import Assignment

DeletionStatus = Literal["HOMOZYGOUS_DELETION", "PRESENT", "INDETERMINATE"]

_MIN_READS_FOR_CALL = 10


def classify_deletion(assignments: list[Assignment]) -> DeletionStatus:
    """Classify a locus's deletion status from its per-read assignment list.

    Returns:
      - ``INDETERMINATE`` when fewer than 10 reads were evaluated.
      - ``HOMOZYGOUS_DELETION`` when ≥10 reads were evaluated and none are
        RESOLVED or PROBABLE (reads overlap the region but no confident
        PSV evidence supports the target locus).
      - ``PRESENT`` when ≥1 assignment is RESOLVED or PROBABLE.
    """
    if len(assignments) < _MIN_READS_FOR_CALL:
        return "INDETERMINATE"
    for a in assignments:
        if a.status in ("RESOLVED", "PROBABLE"):
            return "PRESENT"
    return "HOMOZYGOUS_DELETION"
