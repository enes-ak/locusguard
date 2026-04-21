"""Project per-locus assignments onto a VCF by adding INFO fields."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import cyvcf2

from locusguard.io.vcf import VcfReader, VcfWriter
from locusguard.types import Assignment

LocusRegion = tuple[str, str, int, int]  # (locus_id, chrom, start_1based_incl, end_1based_incl)


_INFO_FIELDS = [
    ("TRUE_LOCUS", "1", "String", "LocusGuard assigned locus"),
    ("LOCUS_CONF", "1", "Float", "LocusGuard confidence in [0,1]"),
    ("LOCUS_STATUS", "1", "String",
     "LocusGuard status: RESOLVED|PROBABLE|AMBIGUOUS|UNASSIGNED"),
    ("LOCUS_EVIDENCE", "1", "String", "LocusGuard evidence decomposition"),
    ("LOCUS_KEY", "1", "String", "LocusGuard locus key for cross-output linkage"),
    ("GENE_CONVERSION_FLAG", "1", "Integer",
     "1 if gene conversion suspected at this locus, 0 otherwise"),
    ("CN_CONTEXT", ".", "String",
     "Comma-separated locus:cn pairs for all configured loci in this run"),
]


class VcfProjector:
    """Reads input VCF, annotates variants inside configured locus regions,
    writes annotated VCF. Variants outside any region are passed through
    unchanged (A mode / pure annotator contract).
    """

    def __init__(
        self,
        input_vcf: Path,
        output_vcf: Path,
        locus_regions: list[LocusRegion],
        cn_by_locus: dict[str, float | None] | None = None,
    ) -> None:
        self._in_path = input_vcf
        self._out_path = output_vcf
        self._regions = locus_regions
        self._cn_by_locus = cn_by_locus or {}

    def run(self, assignments_by_locus: dict[str, list[Assignment]]) -> None:
        per_locus_summary = {
            locus_id: _summarize(a) for locus_id, a in assignments_by_locus.items()
        }
        cn_context_value = self._format_cn_context()

        reader = VcfReader(self._in_path)
        writer = VcfWriter(
            path=self._out_path,
            template_reader=reader,
            extra_info_fields=_INFO_FIELDS,
        )

        for variant in reader.iter_variants():
            info_updates = self._info_for_variant(variant, per_locus_summary)
            if cn_context_value and info_updates:
                info_updates["CN_CONTEXT"] = cn_context_value
            if info_updates:
                writer.write_annotated(variant, info_updates)
            else:
                writer.write_annotated(variant, {})

        writer.close()

    def _format_cn_context(self) -> str:
        if not self._cn_by_locus:
            return ""
        parts = []
        for locus_id, cn in sorted(self._cn_by_locus.items()):
            if cn is None:
                parts.append(f"{locus_id}:na")
            else:
                parts.append(f"{locus_id}:{cn:.2f}")
        return ",".join(parts)

    def _info_for_variant(
        self,
        variant: cyvcf2.Variant,
        per_locus_summary: dict[str, dict[str, Any]],
    ) -> dict[str, object]:
        for locus_id, chrom, start, end in self._regions:
            if chrom != variant.CHROM:
                continue
            if start <= variant.POS <= end:
                summary = per_locus_summary.get(locus_id)
                if summary is None:
                    return {}
                return {
                    "TRUE_LOCUS": summary["best_locus"],
                    "LOCUS_CONF": summary["mean_conf"],
                    "LOCUS_STATUS": summary["dominant_status"],
                    "LOCUS_EVIDENCE": summary["evidence_summary"],
                    "LOCUS_KEY": summary["locus_key"],
                    "GENE_CONVERSION_FLAG": summary["gene_conv_flag"],
                }
        return {}


def _summarize(assignments: list[Assignment]) -> dict[str, object]:
    if not assignments:
        return {
            "best_locus": "",
            "mean_conf": 0.0,
            "dominant_status": "UNASSIGNED",
            "evidence_summary": "na",
            "locus_key": "",
            "gene_conv_flag": 0,
        }
    status_counts = Counter(a.status for a in assignments)
    dominant_status = status_counts.most_common(1)[0][0]
    mean_conf = sum(a.confidence for a in assignments) / len(assignments)
    best_locus_counts = Counter(a.assigned_locus for a in assignments if a.assigned_locus)
    best_locus = best_locus_counts.most_common(1)[0][0] if best_locus_counts else ""
    locus_key = assignments[0].locus_key
    evidence_summary = _format_evidence_summary(assignments)
    gene_conv = any("gene_conversion_suspected" in a.flags for a in assignments)
    return {
        "best_locus": best_locus,
        "mean_conf": round(mean_conf, 4),
        "dominant_status": dominant_status,
        "evidence_summary": evidence_summary,
        "locus_key": locus_key,
        "gene_conv_flag": 1 if gene_conv else 0,
    }


def _format_evidence_summary(assignments: list[Assignment]) -> str:
    """Aggregate evidence availability / match rate across reads, summarize."""
    if not assignments:
        return "na"
    sources: dict[str, list[float]] = {}
    for a in assignments:
        for ev in a.evidence_scores:
            if ev.available:
                sources.setdefault(ev.source, []).append(ev.normalized)
    if not sources:
        return "na"
    parts = []
    for src, values in sources.items():
        mean = sum(values) / len(values)
        parts.append(f"{src}:{mean:.2f}(n={len(values)})")
    return "|".join(parts)
