"""Write run manifest JSON (versions, command, hashes, warnings)."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


def write_manifest(
    output_path: Path,
    locusguard_version: str,
    command_line: str,
    reference_fasta_path: str,
    reference_fasta_md5: str,
    config_hashes: dict[str, str],
    tech: str,
    data_type: str,
    profile_used: str | None,
    runtime_seconds: float,
    warnings: list[str],
    degradations: list[dict[str, str]] | None = None,
    cn_method: str | None = None,
    cn_controls_used: dict[str, list[str]] | None = None,
) -> None:
    doc = {
        "locusguard_version": locusguard_version,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command_line": command_line,
        "reference_fasta": {
            "path": reference_fasta_path,
            "md5": reference_fasta_md5,
        },
        "config_hashes": config_hashes,
        "tech": tech,
        "data_type": data_type,
        "profile_used": profile_used,
        "runtime_sec": runtime_seconds,
        "warnings": warnings,
        "degradations": degradations or [],
        "cn_method": cn_method,
        "cn_controls_used": cn_controls_used or {},
    }
    output_path.write_text(json.dumps(doc, indent=2, sort_keys=True))
