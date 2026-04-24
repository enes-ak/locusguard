"""Library-facing Python API for LocusGuard.

Users who embed LocusGuard in custom scripts use this module:

    from locusguard.api import Annotator
    from locusguard.config import load_config

    annotator = Annotator(
        configs=[load_config(p) for p in config_paths],
        reference_fasta=Path("/refs/grch38.fa"),
        data_type="wgs",
    )
    annotator.annotate_vcf(bam=..., vcf_in=..., vcf_out=...)
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path

from locusguard import __version__
from locusguard.assigner import LocusAssigner
from locusguard.capture_bed import PsvCoverage, compute_psv_coverage
from locusguard.config.schema import LocusConfig
from locusguard.io.bam import BamReader
from locusguard.io.fasta import FastaReader
from locusguard.io.vcf import VcfReader
from locusguard.preflight import run_preflight
from locusguard.projection.vcf import LocusRegion, VcfProjector
from locusguard.reporting.assignments_tsv import write_assignments_tsv
from locusguard.reporting.haplotypes_tsv import write_haplotypes_tsv
from locusguard.reporting.html_report import write_html_report
from locusguard.reporting.manifest import write_manifest
from locusguard.reporting.summary import write_summary
from locusguard.types import Assignment, HaplotypeCluster

_DATATYPE_TO_PROFILE = {
    "wgs": "ont_wgs",
    "wes": "ont_wes",
    "panel": "ont_panel",
}

_SCOPE_WARNING = (
    "Phase 2.8 release: ONT long-read only. PSV match, haplotype consistency, "
    "MAPQ pattern, soft-clip, and coverage-ratio evidence adapters are active, "
    "plus homozygous-deletion detection. --data-type wgs|wes|panel accepted; "
    "WES and Panel modes silence coverage-ratio adapter and support optional "
    "--capture-bed for PSV coverage validation. WES/Panel weight calibration "
    "pending real-data validation."
)


@dataclass(slots=True)
class AnnotationResult:
    variants_total: int
    variants_annotated: int
    assignments_by_locus: dict[str, list[Assignment]]
    haplotype_clusters_by_locus: dict[str, list[HaplotypeCluster]]


class Annotator:
    """High-level entry point for annotating a VCF from a BAM."""

    def __init__(
        self,
        configs: list[LocusConfig],
        reference_fasta: Path,
        data_type: str,
        capture_bed: Path | None = None,
    ) -> None:
        self._configs = configs
        self._reference_fasta = reference_fasta
        self._data_type = data_type
        self._capture_bed = capture_bed
        self._profile_name = _DATATYPE_TO_PROFILE.get(data_type)

    def annotate_vcf(
        self,
        bam: Path,
        vcf_in: Path,
        vcf_out: Path,
        summary_path: Path | None = None,
        manifest_path: Path | None = None,
        sample_name: str | None = None,
        assignments_tsv_path: Path | None = None,
        haplotypes_tsv_path: Path | None = None,
        html_report_path: Path | None = None,
    ) -> AnnotationResult:
        capture_regions = run_preflight(
            bam=bam, vcf=vcf_in, fasta=self._reference_fasta,
            capture_bed=self._capture_bed,
        )
        start = time.perf_counter()

        assignments_by_locus: dict[str, list[Assignment]] = {}
        clusters_by_locus: dict[str, list[HaplotypeCluster]] = {}
        warnings: list[str] = [_SCOPE_WARNING]
        psv_coverage_by_locus: dict[str, PsvCoverage] | None = None
        if capture_regions is not None:
            psv_coverage_by_locus = {
                cfg.locus.id: compute_psv_coverage(cfg, capture_regions)
                for cfg in self._configs
            }
            for cfg in self._configs:
                cov = psv_coverage_by_locus[cfg.locus.id]
                for missing_name in cov.missing:
                    psv = next(p for p in cfg.psvs if p.name == missing_name)
                    warnings.append(
                        f"{cfg.locus.id}: PSV {missing_name} ({psv.chrom}:{psv.pos}) "
                        f"not in capture bed — confidence for this locus may be reduced"
                    )
        elif self._data_type in ("wes", "panel"):
            warnings.append(
                f"--data-type {self._data_type} without --capture-bed — "
                f"PSV coverage cannot be verified"
            )
        with (
            BamReader(bam) as bam_reader,
            FastaReader(self._reference_fasta) as fasta_reader,
        ):
            for cfg in self._configs:
                assigner = LocusAssigner(cfg, profile_name=self._profile_name)
                assignments_by_locus[cfg.locus.id] = assigner.assign(
                    bam_reader, fasta_reader,
                )
                clusters_by_locus[cfg.locus.id] = assigner.haplotype_clusters
                warnings.extend(assigner.warnings)

        locus_regions: list[LocusRegion] = [
            (cfg.locus.id, cfg.coordinates.primary.chrom,
             cfg.coordinates.primary.start, cfg.coordinates.primary.end)
            for cfg in self._configs
        ]

        projector = VcfProjector(
            input_vcf=vcf_in,
            output_vcf=vcf_out,
            locus_regions=locus_regions,
        )
        projector.run(assignments_by_locus)

        variants_total, variants_annotated, counts_by_locus = _count_variants(
            vcf_in, locus_regions,
        )

        runtime = time.perf_counter() - start

        if summary_path is not None:
            gene_conv_flags, gene_conv_evidence = self._derive_gene_conv_maps(
                clusters_by_locus,
            )
            write_summary(
                output_path=summary_path,
                sample_name=sample_name or _infer_sample_name(vcf_in),
                reference="grch38",
                data_type=self._data_type,
                runtime_seconds=round(runtime, 3),
                assignments_by_locus=assignments_by_locus,
                variant_counts_by_locus=counts_by_locus,
                gene_conv_flags_by_locus=gene_conv_flags,
                gene_conv_evidence_by_locus=gene_conv_evidence,
                psv_coverage_by_locus=psv_coverage_by_locus,
            )

        if manifest_path is not None:
            degradations = self._collect_degradations(assignments_by_locus)
            write_manifest(
                output_path=manifest_path,
                locusguard_version=__version__,
                command_line=f"annotate bam={bam} vcf={vcf_in} out={vcf_out}",
                reference_fasta_path=str(self._reference_fasta),
                reference_fasta_md5=_md5(self._reference_fasta),
                config_hashes={cfg.locus.id: _config_hash(cfg) for cfg in self._configs},
                data_type=self._data_type,
                profile_used=self._profile_name,
                runtime_seconds=round(runtime, 3),
                warnings=warnings,
                degradations=degradations,
            )

        if assignments_tsv_path is not None:
            write_assignments_tsv(assignments_tsv_path, assignments_by_locus)

        if haplotypes_tsv_path is not None:
            write_haplotypes_tsv(haplotypes_tsv_path, clusters_by_locus)

        if html_report_path is not None:
            gene_conv_flags, _ = self._derive_gene_conv_maps(clusters_by_locus)
            write_html_report(
                output_path=html_report_path,
                sample_name=sample_name or _infer_sample_name(vcf_in),
                reference="grch38",
                data_type=self._data_type,
                runtime_seconds=round(runtime, 3),
                locusguard_version=__version__,
                assignments_by_locus=assignments_by_locus,
                clusters_by_locus=clusters_by_locus,
                variant_counts_by_locus=counts_by_locus,
                gene_conv_flags_by_locus=gene_conv_flags,
                warnings=warnings,
                degradations=self._collect_degradations(assignments_by_locus),
                psv_coverage_by_locus=psv_coverage_by_locus,
            )

        return AnnotationResult(
            variants_total=variants_total,
            variants_annotated=variants_annotated,
            assignments_by_locus=assignments_by_locus,
            haplotype_clusters_by_locus=clusters_by_locus,
        )

    def _derive_gene_conv_maps(
        self,
        clusters_by_locus: dict[str, list[HaplotypeCluster]],
    ) -> tuple[dict[str, bool], dict[str, str]]:
        flags: dict[str, bool] = {}
        evidence: dict[str, str] = {}
        for locus_id, clusters in clusters_by_locus.items():
            hits = [c for c in clusters if "gene_conversion_suspected" in c.notes]
            flags[locus_id] = bool(hits)
            if hits:
                parts = []
                for c in hits:
                    hotspots = [n for n in c.notes if n.startswith("hotspot_match:")]
                    suffix = f" at {','.join(hotspots)}" if hotspots else ""
                    parts.append(
                        f"cluster {c.hap_id} ({len(c.supporting_reads)} reads){suffix}"
                    )
                evidence[locus_id] = "; ".join(parts)
        return flags, evidence

    def _collect_degradations(
        self,
        assignments_by_locus: dict[str, list[Assignment]],
    ) -> list[dict[str, str]]:
        """Summarize evidence adapters that were disabled or skipped."""
        seen: dict[str, str] = {}
        for assignments in assignments_by_locus.values():
            for a in assignments:
                for ev in a.evidence_scores:
                    if ev.available:
                        continue
                    reason = ev.raw.get("reason", "unavailable")
                    seen.setdefault(ev.source, str(reason))
        return [{"evidence": k, "reason": v} for k, v in sorted(seen.items())]


def _count_variants(
    vcf_in: Path,
    regions: list[LocusRegion],
) -> tuple[int, int, dict[str, int]]:
    total = 0
    annotated = 0
    counts: dict[str, int] = {locus_id: 0 for locus_id, *_ in regions}
    reader = VcfReader(vcf_in)
    for v in reader.iter_variants():
        total += 1
        for locus_id, chrom, start, end in regions:
            if chrom == v.CHROM and start <= v.POS <= end:
                annotated += 1
                counts[locus_id] += 1
                break
    return total, annotated, counts


def _md5(path: Path) -> str:
    h = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _config_hash(config: LocusConfig) -> str:
    payload = json.dumps(
        config.model_dump(mode="json"),
        sort_keys=True,
    ).encode()
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def _infer_sample_name(vcf_in: Path) -> str:
    samples = VcfReader(vcf_in).samples
    return samples[0] if samples else "unknown"
