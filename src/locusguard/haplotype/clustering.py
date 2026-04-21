"""Exact PSV-pattern hash clustering.

Reads sharing an identical PSV-pattern (with `?` for reach=False positions)
form a single haplotype cluster. Deterministic, debuggable, no parameters.
"""
from __future__ import annotations

from collections.abc import Iterable

from locusguard.types import AnalyzedRead, HaplotypeCluster


def cluster_reads(
    reads: Iterable[AnalyzedRead],
    psv_names: list[str],
) -> list[HaplotypeCluster]:
    """Group reads by their PSV-pattern signature.

    Returns clusters with stable H1, H2... IDs assigned in pattern-sorted order
    so the same input always produces the same cluster IDs.
    """
    pattern_to_reads: dict[str, list[str]] = {}
    pattern_to_bases: dict[str, dict[str, str]] = {}

    for read in reads:
        pattern, bases = _encode_pattern(read, psv_names)
        pattern_to_reads.setdefault(pattern, []).append(read.read_id)
        pattern_to_bases[pattern] = bases

    clusters: list[HaplotypeCluster] = []
    for idx, pattern in enumerate(sorted(pattern_to_reads.keys()), start=1):
        clusters.append(
            HaplotypeCluster(
                hap_id=f"H{idx}",
                supporting_reads=pattern_to_reads[pattern],
                psv_pattern=pattern_to_bases[pattern],
                assigned_locus=None,
                confidence=0.0,
                notes=[],
            )
        )
    return clusters


def _encode_pattern(
    read: AnalyzedRead,
    psv_names: list[str],
) -> tuple[str, dict[str, str]]:
    chars: list[str] = []
    bases: dict[str, str] = {}
    for name in psv_names:
        obs = read.psv_observations.get(name)
        if obs is None or not obs.reach:
            chars.append("?")
            bases[name] = "?"
        else:
            chars.append(obs.base)
            bases[name] = obs.base
    return "".join(chars), bases
