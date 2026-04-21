"""Copy number estimation using control-region-normalized depth.

Primary method: per_locus_cn = 2 * depth(locus) / depth(controls)
Secondary method: paralog_ratio = depth(primary) / depth(first_paralog)
"""
from __future__ import annotations

import math

from locusguard.config.schema import LocusConfig
from locusguard.depth.region import DepthStats
from locusguard.types import CnEstimate

_MIN_DEPTH_FOR_CN = 10.0
_CONTROL_COPIES = 2  # diploid assumption for control regions


def estimate_cn(
    config: LocusConfig,
    depths_by_region: dict[str, DepthStats],
    tech: str,
) -> CnEstimate:
    """Estimate copy number for `config.locus.id` given pre-measured depths.

    `depths_by_region` keys:
      - config.locus.id → primary locus DepthStats
      - each paralog name in config.coordinates.paralogs → paralog DepthStats
      - each control_region.name → control DepthStats
    """
    locus_id = config.locus.id

    if tech != "ont":
        return CnEstimate(
            locus_id=locus_id,
            absolute_cn=None,
            absolute_cn_rounded=None,
            paralog_ratio=None,
            cn_total_family=None,
            method="not_supported_for_tech",
            confidence=0.0,
            status="not_supported_for_tech",
        )

    if not config.control_regions:
        return CnEstimate(
            locus_id=locus_id,
            absolute_cn=None,
            absolute_cn_rounded=None,
            paralog_ratio=None,
            cn_total_family=None,
            method="no_control_regions",
            confidence=0.0,
            status="no_control_regions",
        )

    primary = depths_by_region.get(locus_id)
    # Deletion case: primary may be near zero yet controls fine — allow the
    # graceful "ok" path below. Only flag insufficient_depth for primary
    # missing entirely, or depth in the ambiguous band [1.0, _MIN_DEPTH_FOR_CN).
    if primary is None or (
        primary.mean_depth < _MIN_DEPTH_FOR_CN and primary.mean_depth >= 1.0
    ):
        return CnEstimate(
            locus_id=locus_id,
            absolute_cn=None,
            absolute_cn_rounded=None,
            paralog_ratio=None,
            cn_total_family=None,
            method="insufficient_depth",
            confidence=0.0,
            status="insufficient_depth",
        )

    control_depths = [
        depths_by_region[cr.name].mean_depth
        for cr in config.control_regions
        if cr.name in depths_by_region
    ]
    if not control_depths:
        return CnEstimate(
            locus_id=locus_id,
            absolute_cn=None,
            absolute_cn_rounded=None,
            paralog_ratio=None,
            cn_total_family=None,
            method="no_control_regions",
            confidence=0.0,
            status="no_control_regions",
        )

    control_mean = sum(control_depths) / len(control_depths)
    if control_mean < _MIN_DEPTH_FOR_CN:
        return CnEstimate(
            locus_id=locus_id,
            absolute_cn=None,
            absolute_cn_rounded=None,
            paralog_ratio=None,
            cn_total_family=None,
            method="insufficient_depth",
            confidence=0.0,
            status="insufficient_depth",
        )

    # Handle SMN1-deleted case: primary.mean_depth may be near zero but control is fine.
    if primary.mean_depth < 1.0:
        absolute_cn = 0.0
    else:
        absolute_cn = _CONTROL_COPIES * primary.mean_depth / control_mean

    # Paralog ratio + total family CN
    paralog_name = (
        config.locus.paralogs[0] if config.locus.paralogs else None
    )
    paralog_depth = (
        depths_by_region[paralog_name].mean_depth
        if paralog_name and paralog_name in depths_by_region
        else None
    )
    if paralog_depth is not None and paralog_depth > 0:
        paralog_ratio: float | None = primary.mean_depth / paralog_depth
    else:
        paralog_ratio = None

    if paralog_depth is not None:
        paralog_cn = _CONTROL_COPIES * paralog_depth / control_mean
        cn_total_family: float | None = absolute_cn + paralog_cn
    else:
        cn_total_family = absolute_cn

    confidence = _confidence_from_depth(primary)

    notes = [f"control:{cr.name}" for cr in config.control_regions if cr.name in depths_by_region]

    return CnEstimate(
        locus_id=locus_id,
        absolute_cn=round(absolute_cn, 3),
        absolute_cn_rounded=round(absolute_cn),
        paralog_ratio=round(paralog_ratio, 3) if paralog_ratio is not None else None,
        cn_total_family=round(cn_total_family, 3) if cn_total_family is not None else None,
        method="control_region_normalized",
        confidence=round(confidence, 3),
        status="ok",
        notes=notes,
    )


def _confidence_from_depth(depth_stats: DepthStats) -> float:
    if depth_stats.mean_depth <= 0 or depth_stats.length_bp <= 0:
        return 0.0
    # Poisson-based: relative error = σ/μ  * 1/sqrt(N_bp)
    relative_error = math.sqrt(depth_stats.mean_depth) / math.sqrt(
        depth_stats.length_bp * depth_stats.mean_depth
    )
    return max(0.0, min(1.0, 1.0 - 2.0 * relative_error))
