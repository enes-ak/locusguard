"""Per-region depth measurement using pysam.count_coverage."""
from __future__ import annotations

import statistics
from dataclasses import dataclass

from locusguard.io.bam import BamReader


@dataclass(frozen=True, slots=True)
class DepthStats:
    """Summary statistics of read depth over a region."""
    region_name: str
    chrom: str
    start: int                # 1-based
    end: int                  # 1-based inclusive
    mean_depth: float
    median_depth: float
    length_bp: int
    reads_counted: int


def compute_region_depth(
    bam: BamReader,
    chrom: str,
    start: int,              # 1-based inclusive
    end: int,                # 1-based inclusive
    region_name: str = "",
) -> DepthStats:
    """Compute mean + median depth over a 1-based inclusive region.

    Uses pysam.count_coverage which returns per-base A/C/G/T counts; summed
    across alleles gives per-position total depth. Pure-Python, no subprocess.
    """
    if start > end:
        raise ValueError(f"start ({start}) must be <= end ({end})")

    # pysam uses 0-based half-open; convert from 1-based inclusive.
    py_start = start - 1
    py_end = end  # 1-based inclusive end == 0-based exclusive end

    counts_per_base = bam._bam.count_coverage(
        contig=chrom,
        start=py_start,
        stop=py_end,
    )
    # counts_per_base is a tuple of 4 arrays (A, C, G, T). Sum them per-position.
    per_position = [
        counts_per_base[0][i] + counts_per_base[1][i]
        + counts_per_base[2][i] + counts_per_base[3][i]
        for i in range(py_end - py_start)
    ]
    length_bp = end - start + 1
    if not per_position:
        return DepthStats(
            region_name=region_name,
            chrom=chrom,
            start=start,
            end=end,
            mean_depth=0.0,
            median_depth=0.0,
            length_bp=length_bp,
            reads_counted=0,
        )

    mean_d = sum(per_position) / len(per_position)
    median_d = statistics.median(per_position)

    # Count distinct reads overlapping the region (simpler than tracking per-base)
    reads_set = set()
    for read in bam.fetch(chrom, py_start, py_end):
        if read.query_name is not None:
            reads_set.add(read.query_name)

    return DepthStats(
        region_name=region_name,
        chrom=chrom,
        start=start,
        end=end,
        mean_depth=float(mean_d),
        median_depth=float(median_d),
        length_bp=length_bp,
        reads_counted=len(reads_set),
    )
