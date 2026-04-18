"""Evidence adapters — each computes one evidence score from reads + config."""
from locusguard.evidence.base import EvidenceSource, ReadTech

__all__ = ["EvidenceSource", "ReadTech"]
