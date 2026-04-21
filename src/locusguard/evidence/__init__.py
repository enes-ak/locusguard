"""Evidence adapters — each computes one evidence score from reads + config."""
from locusguard.evidence.base import EvidenceSource, ReadTech
from locusguard.evidence.mapq_pattern import MapqPatternEvidence
from locusguard.evidence.psv import PSVEvidence
from locusguard.evidence.softclip import SoftclipEvidence

__all__ = ["EvidenceSource", "MapqPatternEvidence", "PSVEvidence", "ReadTech", "SoftclipEvidence"]
