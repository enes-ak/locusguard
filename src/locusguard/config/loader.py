"""Load LocusConfig from YAML files."""
from __future__ import annotations

from pathlib import Path

import yaml

from locusguard.config.schema import LocusConfig


def load_config(path: Path | str) -> LocusConfig:
    """Load and validate a locus YAML config.

    Raises:
        FileNotFoundError: the path does not exist
        yaml.YAMLError: the YAML is malformed
        pydantic.ValidationError: the document does not satisfy the schema
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Config file not found: {p}")
    with p.open("r", encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    return LocusConfig.model_validate(doc)
