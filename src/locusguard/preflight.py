"""Pre-flight input validation - fail fast with actionable messages."""
from __future__ import annotations

from pathlib import Path

import cyvcf2
import pysam

from locusguard.capture_bed import CaptureBedError, CaptureRegion, load_capture_bed


class PreflightError(RuntimeError):
    """A pre-flight check failed; the run must not proceed."""


def run_preflight(
    bam: Path,
    vcf: Path | None,
    fasta: Path,
    capture_bed: Path | None = None,
) -> list[CaptureRegion] | None:
    """Run pre-flight validation. Returns parsed capture regions when
    ``capture_bed`` is provided, otherwise ``None``.

    Raises ``PreflightError`` on any failure. BAM/FASTA/VCF checks unchanged
    from prior behavior; ``capture_bed`` is validated if provided.
    """
    _check_bam(bam)
    _check_fasta(fasta)
    _check_bam_fasta_compat(bam, fasta)
    if vcf is not None:
        _check_vcf(vcf)

    if capture_bed is None:
        return None
    if not capture_bed.exists():
        raise PreflightError(f"Capture bed file not found: {capture_bed}")
    try:
        return load_capture_bed(capture_bed)
    except CaptureBedError as e:
        raise PreflightError(
            f"Failed to parse capture bed {capture_bed}: {e}"
        ) from e


def _check_bam(bam: Path) -> None:
    if not bam.exists():
        raise PreflightError(f"BAM file not found: {bam}")
    bai = bam.with_suffix(bam.suffix + ".bai")
    csi = bam.with_suffix(bam.suffix + ".csi")
    if not bai.exists() and not csi.exists():
        raise PreflightError(
            f"BAM index missing for {bam.name}. Run: samtools index {bam}"
        )
    with pysam.AlignmentFile(str(bam), "rb") as b:
        hd = b.header.get("HD", {})  # type: ignore[attr-defined]  # pysam AlignmentHeader supports dict-like .get at runtime
        if hd.get("SO") != "coordinate":
            raise PreflightError(
                f"BAM {bam.name} is not coordinate-sorted (SO={hd.get('SO', 'unset')}). "
                f"Run: samtools sort -o sorted.bam {bam}"
            )


def _check_fasta(fasta: Path) -> None:
    if not fasta.exists():
        raise PreflightError(f"Reference FASTA not found: {fasta}")
    fai = fasta.with_suffix(fasta.suffix + ".fai")
    if not fai.exists():
        pysam.faidx(str(fasta))


def _check_bam_fasta_compat(bam: Path, fasta: Path) -> None:
    with pysam.AlignmentFile(str(bam), "rb") as b, pysam.FastaFile(str(fasta)) as f:
        bam_chroms = set(b.references)
        fasta_chroms = set(f.references)
        missing = bam_chroms - fasta_chroms
        if missing:
            raise PreflightError(
                f"BAM contains chromosome(s) not found in FASTA: {sorted(missing)[:3]}... "
                f"This indicates a reference mismatch."
            )


def _check_vcf(vcf: Path) -> None:
    if not vcf.exists():
        raise PreflightError(f"VCF file not found: {vcf}")
    tbi = Path(str(vcf) + ".tbi")
    csi = Path(str(vcf) + ".csi")
    if not tbi.exists() and not csi.exists():
        raise PreflightError(
            f"VCF index missing for {vcf.name}. Run: tabix -p vcf {vcf}"
        )
    v = cyvcf2.VCF(str(vcf))
    if not list(v.samples):
        raise PreflightError(f"VCF {vcf.name} has no samples")
