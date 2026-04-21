"""Write per-locus summary JSON."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from locusguard.types import Assignment


def write_summary(
    output_path: Path,
    sample_name: str,
    reference: str,
    tech: str,
    data_type: str,
    runtime_seconds: float,
    assignments_by_locus: dict[str, list[Assignment]],
    variant_counts_by_locus: dict[str, int],
    gene_conv_flags_by_locus: dict[str, bool] | None = None,
    gene_conv_evidence_by_locus: dict[str, str] | None = None,
) -> None:
    gene_conv_flags_by_locus = gene_conv_flags_by_locus or {}
    gene_conv_evidence_by_locus = gene_conv_evidence_by_locus or {}

    loci_block = {}
    for locus_id, assignments in assignments_by_locus.items():
        block = _per_locus_block(
            locus_id=locus_id,
            assignments=assignments,
            variants_annotated=variant_counts_by_locus.get(locus_id, 0),
        )
        block["gene_conv_flag"] = gene_conv_flags_by_locus.get(locus_id, False)
        if locus_id in gene_conv_evidence_by_locus:
            block["gene_conv_evidence"] = gene_conv_evidence_by_locus[locus_id]
        loci_block[locus_id] = block

    doc = {
        "sample": sample_name,
        "reference": reference,
        "tech": tech,
        "data_type": data_type,
        "runtime_sec": runtime_seconds,
        "loci": loci_block,
    }
    output_path.write_text(json.dumps(doc, indent=2, sort_keys=True))


def _per_locus_block(
    locus_id: str,
    assignments: list[Assignment],
    variants_annotated: int,
) -> dict[str, object]:
    if not assignments:
        return {
            "status": "UNASSIGNED",
            "read_count": 0,
            "variants_annotated": variants_annotated,
            "mean_confidence": 0.0,
        }

    status_counts = Counter(a.status for a in assignments)
    dominant_status = status_counts.most_common(1)[0][0]
    mean_conf = sum(a.confidence for a in assignments) / len(assignments)

    flag_counts: Counter[str] = Counter()
    for a in assignments:
        for f in a.flags:
            flag_counts[f] += 1

    return {
        "status": dominant_status,
        "read_count": len(assignments),
        "variants_annotated": variants_annotated,
        "mean_confidence": round(mean_conf, 4),
        "status_counts": dict(status_counts),
        "flag_counts": dict(flag_counts),
    }
