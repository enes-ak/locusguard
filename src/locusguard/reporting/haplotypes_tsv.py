"""Per-cluster haplotypes TSV writer."""
from __future__ import annotations

import csv
from pathlib import Path

from locusguard.types import HaplotypeCluster


_HEADER = [
    "hap_id",
    "locus_config",
    "psv_pattern",
    "supporting_read_count",
    "assigned_locus",
    "cluster_confidence",
    "notes",
    "low_support",
]

_LOW_SUPPORT_THRESHOLD = 2


def write_haplotypes_tsv(
    output_path: Path,
    clusters_by_locus: dict[str, list[HaplotypeCluster]],
) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(_HEADER)
        for locus_cfg_id, clusters in clusters_by_locus.items():
            for c in clusters:
                pattern = ";".join(f"{k}={v}" for k, v in sorted(c.psv_pattern.items()))
                low_support = (
                    "true" if len(c.supporting_reads) < _LOW_SUPPORT_THRESHOLD else "false"
                )
                writer.writerow([
                    c.hap_id,
                    locus_cfg_id,
                    pattern,
                    len(c.supporting_reads),
                    c.assigned_locus or "",
                    f"{c.confidence:.4f}",
                    ",".join(c.notes),
                    low_support,
                ])
