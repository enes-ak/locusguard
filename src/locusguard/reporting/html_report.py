"""HTML report writer using Jinja2.

Single-template layout: header + per-locus cards + warnings + manifest footer.
Inline CSS, inline SVG bar chart, no external stylesheets or JS.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from jinja2 import Environment, select_autoescape

from locusguard.types import Assignment, HaplotypeCluster

_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>LocusGuard report — {{ sample_name }}</title>
<style>
  body { font-family: -apple-system, Segoe UI, sans-serif; max-width: 1000px; margin: 2em auto; color: #222; }
  header { border-bottom: 2px solid #444; padding-bottom: 1em; margin-bottom: 2em; }
  h1 { margin: 0; font-size: 1.5em; }
  .meta { color: #666; font-size: 0.9em; }
  .locus-card { border: 1px solid #ddd; border-radius: 8px; padding: 1em 1.5em; margin-bottom: 1.5em; background: #fafafa; }
  .locus-card h2 { margin: 0 0 0.5em 0; font-size: 1.2em; display: flex; align-items: center; }
  .status { display: inline-block; padding: 2px 10px; border-radius: 4px; font-size: 0.8em; margin-left: 0.8em; color: white; }
  .status-RESOLVED { background: #28a745; }
  .status-PROBABLE { background: #f0ad4e; }
  .status-AMBIGUOUS { background: #dc3545; }
  .status-UNASSIGNED { background: #6c757d; }
  table { border-collapse: collapse; width: 100%; font-size: 0.9em; margin-top: 0.5em; }
  th, td { border: 1px solid #ddd; padding: 6px 10px; text-align: left; }
  th { background: #eef; }
  .alert { background: #fff3cd; border: 1px solid #f0ad4e; padding: 0.8em; border-radius: 4px; margin-top: 0.5em; }
  .alert strong { color: #856404; }
  .warnings { background: #f8f9fa; border-left: 4px solid #6c757d; padding: 0.8em 1em; font-size: 0.9em; }
  footer { color: #888; font-size: 0.8em; margin-top: 3em; border-top: 1px solid #ddd; padding-top: 1em; }
</style>
</head>
<body>
<header>
  <h1>LocusGuard report — {{ sample_name }}</h1>
  <div class="meta">
    Reference: {{ reference }} · Tech: {{ tech }} · Data type: {{ data_type }} · Runtime: {{ '%.1f'|format(runtime_seconds) }} s
  </div>
</header>

{% for locus_id in locus_ids %}
  {% set summary = locus_summaries[locus_id] %}
  <div class="locus-card">
    <h2>{{ locus_id }}
      <span class="status status-{{ summary.status }}">{{ summary.status }}</span>
    </h2>
    <p>
      <strong>Mean confidence:</strong> {{ '%.2f'|format(summary.mean_conf) }} ·
      <strong>Reads assigned:</strong> {{ summary.read_count }} ·
      <strong>Variants annotated:</strong> {{ summary.variants }}
    </p>

    {% if summary.gene_conv_flag %}
    <div class="alert">
      <strong>Gene conversion suspected.</strong>
      {% if summary.hotspot_names %}Hotspot(s): {{ summary.hotspot_names|join(', ') }}.{% endif %}
      Orthogonal validation recommended.
    </div>
    {% endif %}

    {% if summary.clusters %}
    <table>
      <thead><tr><th>Cluster</th><th>Reads</th><th>PSV pattern</th><th>Assigned</th><th>Notes</th></tr></thead>
      <tbody>
      {% for c in summary.clusters %}
      <tr>
        <td>{{ c.hap_id }}</td>
        <td>{{ c.read_count }}</td>
        <td><code>{{ c.pattern }}</code></td>
        <td>{{ c.assigned or '—' }}</td>
        <td>{{ c.notes|join(', ') }}</td>
      </tr>
      {% endfor %}
      </tbody>
    </table>
    {% endif %}
  </div>
{% endfor %}

{% if warnings %}
<div class="warnings">
  <strong>Warnings</strong>
  <ul>{% for w in warnings %}<li>{{ w }}</li>{% endfor %}</ul>
</div>
{% endif %}

{% if degradations %}
<div class="warnings">
  <strong>Evidence degradations</strong>
  <ul>{% for d in degradations %}<li><code>{{ d.evidence }}</code> — {{ d.reason }}</li>{% endfor %}</ul>
</div>
{% endif %}

<footer>
  LocusGuard {{ locusguard_version }}
</footer>
</body>
</html>
"""


def write_html_report(
    output_path: Path,
    sample_name: str,
    reference: str,
    tech: str,
    data_type: str,
    runtime_seconds: float,
    locusguard_version: str,
    assignments_by_locus: dict[str, list[Assignment]],
    clusters_by_locus: dict[str, list[HaplotypeCluster]],
    variant_counts_by_locus: dict[str, int],
    gene_conv_flags_by_locus: dict[str, bool],
    warnings: list[str],
    degradations: list[dict[str, str]],
) -> None:
    env = Environment(autoescape=select_autoescape(["html"]))
    template = env.from_string(_TEMPLATE)

    locus_ids = sorted(assignments_by_locus.keys())
    locus_summaries: dict[str, dict[str, object]] = {}
    for locus_id in locus_ids:
        assignments = assignments_by_locus.get(locus_id, [])
        clusters = clusters_by_locus.get(locus_id, [])
        status = _dominant_status(assignments)
        mean_conf = sum(a.confidence for a in assignments) / len(assignments) if assignments else 0.0
        hotspot_names = []
        for c in clusters:
            for note in c.notes:
                if note.startswith("hotspot_match:"):
                    hotspot_names.append(note.split(":", 1)[1])
        locus_summaries[locus_id] = {
            "status": status,
            "mean_conf": mean_conf,
            "read_count": len(assignments),
            "variants": variant_counts_by_locus.get(locus_id, 0),
            "gene_conv_flag": gene_conv_flags_by_locus.get(locus_id, False),
            "hotspot_names": sorted(set(hotspot_names)),
            "clusters": [
                {
                    "hap_id": c.hap_id,
                    "read_count": len(c.supporting_reads),
                    "pattern": ";".join(f"{k}={v}" for k, v in sorted(c.psv_pattern.items())),
                    "assigned": c.assigned_locus or "",
                    "notes": c.notes,
                }
                for c in clusters
            ],
        }

    html = template.render(
        sample_name=sample_name,
        reference=reference,
        tech=tech,
        data_type=data_type,
        runtime_seconds=runtime_seconds,
        locusguard_version=locusguard_version,
        locus_ids=locus_ids,
        locus_summaries=locus_summaries,
        warnings=warnings,
        degradations=degradations,
    )
    output_path.write_text(html, encoding="utf-8")


def _dominant_status(assignments: list[Assignment]) -> str:
    if not assignments:
        return "UNASSIGNED"
    counts = Counter(a.status for a in assignments)
    return counts.most_common(1)[0][0]
