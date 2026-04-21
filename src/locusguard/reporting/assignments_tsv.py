"""Per-read assignments TSV writer."""
from __future__ import annotations

import csv
from pathlib import Path

from locusguard.types import Assignment


_HEADER = [
    "read_id",
    "locus_config",
    "assigned_locus",
    "confidence",
    "status",
    "cluster_id",
    "locus_key",
    "evidence_scores_compact",
    "flags",
]


def write_assignments_tsv(
    output_path: Path,
    assignments_by_locus: dict[str, list[Assignment]],
) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(_HEADER)
        for locus_cfg_id, assignments in assignments_by_locus.items():
            for a in assignments:
                writer.writerow([
                    a.read_id,
                    locus_cfg_id,
                    a.assigned_locus or "",
                    f"{a.confidence:.4f}",
                    a.status,
                    a.cluster_id or "",
                    a.locus_key,
                    _compact_evidence(a.evidence_scores),
                    ",".join(sorted(a.flags)),
                ])


def _compact_evidence(scores) -> str:
    parts = []
    for ev in scores:
        if ev.available:
            parts.append(f"{ev.source}:{ev.normalized:.2f}")
        else:
            parts.append(f"{ev.source}:na")
    return "|".join(parts)
