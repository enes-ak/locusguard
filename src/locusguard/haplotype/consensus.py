"""Per-cluster locus assignment and gene-conversion detection."""
from __future__ import annotations

from locusguard.config.schema import LocusConfig
from locusguard.types import HaplotypeCluster

# Gene conversion requires at least two candidate paralogs to compare against.
_MIN_PARALOGS_FOR_GENE_CONV = 2


def assign_cluster_locus(cluster: HaplotypeCluster, config: LocusConfig) -> None:
    """Assign a locus to a cluster based on PSV pattern agreement.

    Score per candidate locus = (matches among observed PSVs) / (observed PSVs).
    Best-scoring locus wins. Confidence reflects both agreement rate and
    coverage (fraction of PSVs actually observed).
    """
    candidates = _candidate_loci(config)
    observed = {k: v for k, v in cluster.psv_pattern.items() if v != "?"}

    if not observed:
        cluster.assigned_locus = None
        cluster.confidence = 0.0
        return

    best_locus: str | None = None
    best_score = -1.0
    for locus_id in candidates:
        matches = 0
        total = 0
        for psv in config.psvs:
            expected = psv.alleles.get(locus_id)
            if expected is None or psv.name not in observed:
                continue
            total += 1
            if observed[psv.name] == expected:
                matches += 1
        if total == 0:
            continue
        score = matches / total
        if score > best_score:
            best_score = score
            best_locus = locus_id

    cluster.assigned_locus = best_locus
    coverage = len(observed) / max(1, len(config.psvs))
    cluster.confidence = best_score * coverage


def detect_gene_conversion(cluster: HaplotypeCluster, config: LocusConfig) -> None:
    """Add notes to cluster when PSV pattern is internally inconsistent.

    - 'gene_conversion_suspected' when some PSVs match one paralog and
      others match a different paralog.
    - 'hotspot_match:<name>' when the mixed positions fall inside a
      configured known hotspot.
    """
    candidates = _candidate_loci(config)
    if cluster.assigned_locus is None or len(candidates) < _MIN_PARALOGS_FOR_GENE_CONV:
        return

    mismatched_psv_names = _find_mismatched_psvs(cluster, config, candidates)
    if not mismatched_psv_names:
        return

    if "gene_conversion_suspected" not in cluster.notes:
        cluster.notes.append("gene_conversion_suspected")

    if config.gene_conversion is None:
        return

    name_to_pos = {p.name: p.pos for p in config.psvs}
    for hotspot in config.gene_conversion.known_hotspots:
        hit = any(
            hotspot.start <= name_to_pos[n] <= hotspot.end
            for n in mismatched_psv_names
            if n in name_to_pos
        )
        if hit:
            note = f"hotspot_match:{hotspot.name}"
            if note not in cluster.notes:
                cluster.notes.append(note)


def _find_mismatched_psvs(
    cluster: HaplotypeCluster,
    config: LocusConfig,
    candidates: list[str],
) -> list[str]:
    """Return PSV names whose observed base matches a paralog other than the assigned locus."""
    mismatched: list[str] = []
    assigned = cluster.assigned_locus
    for psv in config.psvs:
        observed = cluster.psv_pattern.get(psv.name)
        if observed is None or observed == "?":
            continue
        expected_for_assigned = psv.alleles.get(assigned) if assigned is not None else None
        if expected_for_assigned is None or observed == expected_for_assigned:
            continue
        # Does this observed base match any OTHER paralog's allele?
        for other_locus in candidates:
            if other_locus == assigned:
                continue
            if psv.alleles.get(other_locus) == observed:
                mismatched.append(psv.name)
                break
    return mismatched


def _candidate_loci(config: LocusConfig) -> list[str]:
    candidates = {config.locus.id}
    for psv in config.psvs:
        candidates.update(psv.alleles.keys())
    return sorted(candidates)
