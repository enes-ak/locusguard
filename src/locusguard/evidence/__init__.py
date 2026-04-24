"""Evidence adapters — each computes one evidence score from reads + config."""
from locusguard.evidence.base import EvidenceSource
from locusguard.evidence.coverage_ratio import CoverageRatioEvidence
from locusguard.evidence.haplotype_consistency import HaplotypeConsistencyEvidence
from locusguard.evidence.mapq_pattern import MapqPatternEvidence
from locusguard.evidence.psv import PSVEvidence
from locusguard.evidence.softclip import SoftclipEvidence

__all__ = [
    "CoverageRatioEvidence",
    "EvidenceSource",
    "HaplotypeConsistencyEvidence",
    "MapqPatternEvidence",
    "PSVEvidence",
    "SoftclipEvidence",
]
