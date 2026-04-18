"""Thin pysam.FastaFile wrapper for position-level reads."""
from __future__ import annotations

from pathlib import Path

import pysam


class FastaReader:
    """Reads bases from an indexed FASTA.

    Expects a `.fai` index alongside the FASTA; if missing, builds one via pysam.
    """

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)
        fai = self._path.with_suffix(self._path.suffix + ".fai")
        if not fai.exists():
            pysam.faidx(str(self._path))
        self._fa = pysam.FastaFile(str(self._path))

    def base_at(self, chrom: str, pos: int) -> str:
        """Return the single base at 0-based position `pos` on `chrom`."""
        if chrom not in self.chromosomes():
            raise KeyError(f"Chromosome '{chrom}' not found in FASTA")
        return self._fa.fetch(chrom, pos, pos + 1).upper()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Return the sequence on `chrom` between 0-based `[start, end)`."""
        return self._fa.fetch(chrom, start, end).upper()

    def chromosomes(self) -> list[str]:
        return list(self._fa.references)

    def close(self) -> None:
        self._fa.close()

    def __enter__(self) -> FastaReader:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()
