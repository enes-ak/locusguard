"""cyvcf2 wrappers: streaming VCF read and annotated write."""
from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import cyvcf2


class VcfReader:
    """Thin wrapper over cyvcf2.VCF with iter/fetch."""

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)
        self._vcf = cyvcf2.VCF(str(self._path))

    def iter_variants(self) -> Iterator[cyvcf2.Variant]:
        for v in self._vcf:
            yield v

    def fetch(self, chrom: str, start: int, end: int) -> Iterator[cyvcf2.Variant]:
        region = f"{chrom}:{start}-{end}"
        for v in self._vcf(region):
            yield v

    @property
    def samples(self) -> list[str]:
        return list(self._vcf.samples)

    @property
    def raw_header(self) -> str:
        return self._vcf.raw_header

    @property
    def cyvcf2_handle(self) -> cyvcf2.VCF:
        return self._vcf


InfoFieldDef = tuple[str, str, str, str]  # (id, number, type, description)


class VcfWriter:
    """Write annotated VCF by adding INFO fields to variants from a template."""

    def __init__(
        self,
        path: Path | str,
        template_reader: VcfReader,
        extra_info_fields: list[InfoFieldDef],
    ) -> None:
        handle = template_reader.cyvcf2_handle
        for field_id, number, type_, description in extra_info_fields:
            handle.add_info_to_header({
                "ID": field_id,
                "Number": number,
                "Type": type_,
                "Description": description,
            })
        self._writer = cyvcf2.Writer(str(path), handle)
        self._path = Path(path)

    def write_annotated(
        self,
        variant: cyvcf2.Variant,
        info_updates: dict[str, object],
    ) -> None:
        for key, value in info_updates.items():
            variant.INFO[key] = value
        self._writer.write_record(variant)

    def close(self) -> None:
        self._writer.close()
