"""Thin pysam.AlignmentFile wrapper for region iteration."""
from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pysam


class BamReader:
    """Reads alignments from an indexed BAM/CRAM file.

    Requires a `.bai` / `.csi` index. Raises on missing index.
    """

    _LONG_READ_MEDIAN_THRESHOLD = 500  # bp; above this median we call it long-read

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)
        self._bam = pysam.AlignmentFile(str(self._path), "rb", require_index=True)

    def fetch(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> Iterator[pysam.AlignedSegment]:
        yield from self._bam.fetch(chrom, start, end)

    def chromosomes(self) -> list[str]:
        return list(self._bam.references)

    @property
    def is_sorted_by_coordinate(self) -> bool:
        hd = self._bam.header.get("HD", {})  # type: ignore[attr-defined]  # pysam AlignmentHeader supports dict-like .get at runtime
        return bool(hd.get("SO") == "coordinate")

    def estimated_is_long_read(self, sample_size: int = 100) -> bool:
        """Peek at up to `sample_size` reads; return True if median length > threshold."""
        lengths: list[int] = []
        for read in self._bam.head(sample_size):
            if read.query_length:
                lengths.append(read.query_length)
        if not lengths:
            return False
        lengths.sort()
        median = lengths[len(lengths) // 2]
        return median > self._LONG_READ_MEDIAN_THRESHOLD

    def close(self) -> None:
        self._bam.close()

    def __enter__(self) -> BamReader:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()
