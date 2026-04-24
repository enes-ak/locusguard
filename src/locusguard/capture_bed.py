"""Capture bed parsing + PSV intersection.

This module supports the Phase 2.8 capture-aware flow for ONT WES and Panel
inputs. Given a BED3+ file describing capture regions and a LocusConfig's
PSV positions, it reports which PSVs fall within the capture.

Coordinate conventions:
  - BED intervals are 0-based half-open: [start, end).
  - PSV positions in LocusConfig are 1-based.
  - A 1-based position ``pos`` is within a BED region when
    ``start < pos <= end`` (equivalent to the 0-based position ``pos - 1``
    being in ``[start, end)``).
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from locusguard.config.schema import LocusConfig

_MIN_BED_COLUMNS = 3


class CaptureBedError(RuntimeError):
    """A capture bed file failed to parse or validate."""


@dataclass(frozen=True, slots=True)
class CaptureRegion:
    """A single capture region (BED row)."""

    chrom: str
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive


@dataclass(frozen=True, slots=True)
class PsvCoverage:
    """Per-locus PSV coverage result.

    ``covered`` and ``missing`` are lists by convention; callers must not
    mutate them after construction. ``fraction_covered`` is
    ``len(covered) / total``, or ``0.0`` when the config defines no PSVs.
    """

    covered: list[str]
    missing: list[str]
    fraction_covered: float


def load_capture_bed(path: Path) -> list[CaptureRegion]:
    """Parse a BED3+ file into a list of ``CaptureRegion``.

    Skips track lines, browser lines, comment lines (``#``), and blank lines.
    Reads only the first three columns; additional columns are ignored.
    Raises ``CaptureBedError`` on malformed rows (non-int coordinates or
    fewer than 3 columns). Raises ``FileNotFoundError`` if the path does
    not exist.
    """
    regions: list[CaptureRegion] = []
    with path.open("r") as fh:
        for line_num, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                continue
            if stripped.startswith("track") or stripped.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < _MIN_BED_COLUMNS:
                raise CaptureBedError(
                    f"Line {line_num}: expected >={_MIN_BED_COLUMNS} tab-separated columns,"
                    f" got {len(parts)}"
                )
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError as e:
                raise CaptureBedError(
                    f"Line {line_num}: start/end must be integers,"
                    f" got '{parts[1]}' / '{parts[2]}'"
                ) from e
            regions.append(CaptureRegion(chrom=chrom, start=start, end=end))
    return regions


def position_in_capture(
    chrom: str, pos_1based: int, regions: Iterable[CaptureRegion],
) -> bool:
    """Return True iff the 1-based position ``pos_1based`` on ``chrom`` falls
    within any region.

    A 1-based position is in a 0-based half-open BED region ``[start, end)``
    iff ``start < pos_1based <= end``.
    """
    return any(r.chrom == chrom and r.start < pos_1based <= r.end for r in regions)


def compute_psv_coverage(
    config: LocusConfig, regions: list[CaptureRegion],
) -> PsvCoverage:
    """Compute per-locus PSV coverage against the capture regions.

    Returns a ``PsvCoverage`` listing covered and missing PSV names and the
    fraction ``len(covered) / len(config.psvs)`` (0.0 when the config
    defines no PSVs).
    """
    covered: list[str] = []
    missing: list[str] = []
    for psv in config.psvs:
        if position_in_capture(psv.chrom, psv.pos, regions):
            covered.append(psv.name)
        else:
            missing.append(psv.name)
    total = len(config.psvs)
    fraction = (len(covered) / total) if total > 0 else 0.0
    return PsvCoverage(
        covered=covered, missing=missing, fraction_covered=fraction,
    )
