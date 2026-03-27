#!/usr/bin/env python3
"""
Create comprehensive Cordlife Quality Assessment Report with embedded figures (base64).

This is intentionally aligned with the style/structure of:
  - `bin/scripts/utils/reports/generate_ucl_report_final.py`

Inputs (expected to already exist from batch analysis):
  - analysis/batch_analysis_cordlife_complete/*.png
  - analysis/batch_analysis_cordlife_complete/*.txt
  - analysis/GC_Chr_Bias_Relationship.png

Output:
  - Cordlife_Quality_Assessment_Report_Final.html
"""

from __future__ import annotations

import base64
import os
import re


def image_to_base64(image_path: str) -> str | None:
    """Convert image to base64 for embedding in HTML."""
    if not os.path.exists(image_path):
        return None
    with open(image_path, "rb") as f:
        return base64.b64encode(f.read()).decode()


def read_text(path: str) -> str:
    if not os.path.exists(path):
        return ""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def parse_total_samples_from_summary(summary_txt: str) -> int | None:
    # Looks like: "Total samples: 901"
    m = re.search(r"Total\s+samples:\s*(\d+)", summary_txt)
    return int(m.group(1)) if m else None


def parse_months_from_summary(summary_txt: str) -> str | None:
    # Looks like: "Months covered: 2507, 2508, ..."
    m = re.search(r"Months\s+covered:\s*([0-9,\s]+)", summary_txt)
    if not m:
        return None
    return m.group(1).strip()


def parse_unique_batches_from_summary(summary_txt: str) -> int | None:
    m = re.search(r"Unique\s+batches:\s*(\d+)", summary_txt)
    return int(m.group(1)) if m else None


def parse_cluster_sizes(pca_cluster_txt: str) -> dict:
    """
    Parse cluster sizes from `pca_clustering_report.txt`:
      Group 2 (Right): 557 samples (61.8%)
      Group 1 (Left): 344 samples (38.2%)
    """
    out = {}
    for line in pca_cluster_txt.splitlines():
        line = line.strip()
        m = re.match(r"Group\s+(\d+)\s+\(([^)]+)\):\s+(\d+)\s+samples\s+\(([\d.]+)%\)", line)
        if m:
            gid = int(m.group(1))
            label = m.group(2)
            n = int(m.group(3))
            pct = float(m.group(4))
            out[gid] = {"label": label, "n": n, "pct": pct}
    return out


def parse_group_metric_means(pca_cluster_txt: str) -> dict:
    """
    Parse means from the "STATISTICAL COMPARISON" table (Group1 vs Group2).
    Example lines:
      duplication_rate(%)                   2.487        3.826    9.88e-140    ***
      mean_coverageData(X)                  0.169        0.255    1.62e-139    ***
    """
    metrics = {}
    for line in pca_cluster_txt.splitlines():
        line = line.rstrip()
        if not line or line.startswith("Metric") or line.startswith("-"):
            continue
        # metric name can have spaces? here it's padded.
        m = re.match(r"^([A-Za-z0-9_\-\(\)%\.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9eE\.\-]+)\s+(\S+)\s*$", line)
        if not m:
            continue
        metric = m.group(1)
        g1 = float(m.group(2))
        g2 = float(m.group(3))
        pval = m.group(4)
        sig = m.group(5)
        metrics[metric] = {"group1": g1, "group2": g2, "pval": pval, "sig": sig}
    return metrics


def create_report_with_figures() -> str:
    # Read existing analysis text outputs for grounded numbers
    summary_txt = read_text("analysis/batch_analysis_cordlife_complete/batch_effect_summary.txt")
    pca_cluster_txt = read_text("analysis/batch_analysis_cordlife_complete/pca_clustering_report.txt")

    total_samples = parse_total_samples_from_summary(summary_txt) or 0
    months_covered = parse_months_from_summary(summary_txt) or "2507, 2508, 2509, 2510, 2511, 2512, 2601"
    unique_batches = parse_unique_batches_from_summary(summary_txt) or 0

    clusters = parse_cluster_sizes(pca_cluster_txt)
    metric_means = parse_group_metric_means(pca_cluster_txt)

    # Key metrics for narrative (fallbacks if parsing fails)
    g1 = clusters.get(1, {"n": None, "pct": None})
    g2 = clusters.get(2, {"n": None, "pct": None})

    dup = metric_means.get("duplication_rate(%)", {})
    cov = metric_means.get("mean_coverageData(X)", {})
    mapr = metric_means.get("mapping_rate(%)", {})
    gc = metric_means.get("GC_content(%)", {})

    images = {
        "cord_pca": image_to_base64("analysis/batch_analysis_cordlife_complete/pca_clustering_groups.png"),
        "cord_dup": image_to_base64("analysis/batch_analysis_cordlife_complete/duplication_rate_timeline.png"),
        "cord_cov": image_to_base64("analysis/batch_analysis_cordlife_complete/coverage_timeline.png"),
        "cord_chr_dist": image_to_base64("analysis/batch_analysis_cordlife_complete/05_chromosome_distribution.png"),
        "cord_batch_over_time": image_to_base64("analysis/batch_analysis_cordlife_complete/04_batch_effect_over_time.png"),
        "cord_readcount_vs_chrprop": image_to_base64("analysis/batch_analysis_cordlife_complete/readcount_vs_chrprop.png"),
        "gc_chr_bias": image_to_base64("analysis/GC_Chr_Bias_Relationship.png"),
    }

    # Basic robustness: if any image is missing, keep placeholder text.
    def img_tag(key: str, alt: str) -> str:
        b64 = images.get(key)
        if not b64:
            return f"<div style='padding:12px;background:#fff3cd;border:1px solid #ffeeba;border-radius:6px;'>Missing image: <code>{key}</code></div>"
        return f"<img src=\"data:image/png;base64,{b64}\" alt=\"{alt}\">"

    report_date = "February 5, 2026"
    analysis_period = months_covered.replace(" ", "")

    # Platform: Cordlife is used as the internal reference lab in other docs; keep wording neutral.
    platform_label = "NextSeq 550 (reference lab dataset)"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NIPT Sample Quality Assessment Report - Cordlife Laboratory</title>
  <style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{
      font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
      line-height: 1.6;
      color: #333;
      background: #f5f5f5;
      padding: 20px;
    }}
    .container {{
      max-width: 1200px;
      margin: 0 auto;
      background: white;
      box-shadow: 0 0 20px rgba(0,0,0,0.1);
    }}
    .header {{
      background: linear-gradient(135deg, #2C3E50 0%, #27AE60 100%);
      color: white;
      padding: 40px;
      text-align: center;
    }}
    .header h1 {{ font-size: 2.5em; margin-bottom: 10px; }}
    .header .subtitle {{ font-size: 1.2em; opacity: 0.9; }}
    .meta-info {{
      background: #34495E;
      color: white;
      padding: 20px 40px;
      display: flex;
      justify-content: space-between;
      flex-wrap: wrap;
    }}
    .meta-item {{ margin: 5px 0; }}
    .content {{ padding: 40px; }}
    .section {{ margin-bottom: 50px; }}
    .section h2 {{
      color: #2C3E50;
      font-size: 2em;
      margin-bottom: 20px;
      padding-bottom: 10px;
      border-bottom: 3px solid #27AE60;
    }}
    .section h3 {{
      color: #34495E;
      font-size: 1.5em;
      margin: 25px 0 15px 0;
    }}
    .executive-summary {{
      background: #E8F5E9;
      border-left: 5px solid #27AE60;
      padding: 25px;
      margin: 20px 0;
      font-size: 1.1em;
    }}
    .observation-box {{
      background: #F8F9FA;
      border-left: 5px solid #6C757D;
      padding: 25px;
      margin: 20px 0;
    }}
    .observation-box h4 {{
      color: #495057;
      font-size: 1.3em;
      margin-bottom: 15px;
    }}
    .success-box {{
      background: #D4EDDA;
      border-left: 5px solid #28A745;
      padding: 25px;
      margin: 20px 0;
    }}
    .success-box h4 {{
      color: #28A745;
      font-size: 1.3em;
      margin-bottom: 15px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin: 20px 0;
      box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }}
    table th {{
      background: #2C3E50;
      color: white;
      padding: 15px;
      text-align: left;
      font-weight: 600;
    }}
    table td {{
      padding: 12px 15px;
      border-bottom: 1px solid #ddd;
    }}
    table tr:hover {{ background: #f5f5f5; }}
    .metric-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
      gap: 20px;
      margin: 25px 0;
    }}
    .metric-card {{
      background: white;
      border: 2px solid #ddd;
      border-radius: 8px;
      padding: 20px;
      text-align: center;
      box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }}
    .metric-card.good {{ border-color: #28A745; background: #F5FFF5; }}
    .metric-card h4 {{
      font-size: 0.9em;
      color: #666;
      margin-bottom: 10px;
      text-transform: uppercase;
    }}
    .metric-card .value {{ font-size: 2.5em; font-weight: bold; margin: 10px 0; color: #28A745; }}
    .metric-card .label {{ font-size: 0.9em; color: #666; }}
    .figure-container {{
      margin: 30px 0;
      padding: 20px;
      background: #F8F9FA;
      border-radius: 8px;
    }}
    .figure-container img {{
      max-width: 100%;
      height: auto;
      border: 1px solid #ddd;
      border-radius: 8px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.1);
      display: block;
      margin: 20px auto;
    }}
    .figure-caption {{
      text-align: center;
      font-style: italic;
      color: #666;
      margin-top: 10px;
      font-size: 0.95em;
    }}
    .figure-title {{
      font-weight: bold;
      color: #2C3E50;
      margin-bottom: 15px;
      font-size: 1.1em;
    }}
    .recommendations {{
      background: #E8F5E9;
      border-left: 5px solid #27AE60;
      padding: 25px;
      margin: 20px 0;
    }}
    .recommendations h4 {{
      color: #1E8449;
      font-size: 1.3em;
      margin-bottom: 15px;
    }}
    .recommendations ul {{ margin-left: 20px; }}
    .recommendations li {{ margin: 10px 0; line-height: 1.8; }}
    .footer {{
      background: #2C3E50;
      color: white;
      padding: 30px 40px;
      text-align: center;
    }}
    .footer p {{ margin: 5px 0; opacity: 0.85; }}
    .page-break {{ page-break-after: always; }}
    @media print {{
      body {{ background: white; padding: 0; }}
      .container {{ box-shadow: none; }}
      .page-break {{ page-break-after: always; }}
    }}
    .status-good {{ color: #28A745; font-weight: bold; }}
    code {{ background: #f1f1f1; padding: 2px 6px; border-radius: 4px; }}
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <h1>NIPT Sample Quality Assessment Report</h1>
      <div class="subtitle">Cordlife dataset (batch-level and sample-level QC patterns)</div>
    </div>

    <div class="meta-info">
      <div class="meta-item"><strong>Laboratory:</strong> Cordlife</div>
      <div class="meta-item"><strong>Analysis Period:</strong> {analysis_period}</div>
      <div class="meta-item"><strong>Total Samples:</strong> {total_samples}</div>
      <div class="meta-item"><strong>Unique Batches:</strong> {unique_batches}</div>
      <div class="meta-item"><strong>Sequencing Platform:</strong> {platform_label}</div>
      <div class="meta-item"><strong>Report Date:</strong> {report_date}</div>
    </div>

    <div class="content">
      <div class="section">
        <h2>Executive Summary</h2>

        <div class="executive-summary">
          <p><strong>Purpose:</strong> This report summarizes Cordlife sample quality patterns and batch-level variation using the existing batch analysis outputs.</p>
          <p style="margin-top:10px;"><strong>Scope:</strong> {total_samples} samples across months <code>{months_covered}</code> (unique batches: {unique_batches}).</p>
        </div>

        <div class="success-box">
          <h4>Key Summary</h4>
          <ul style="margin-left: 20px; margin-top: 10px;">
            <li><strong>Two major clusters are observed</strong> in PCA clustering. Separation is driven mainly by <strong>duplication rate</strong> and <strong>coverage</strong> (both highly significant).</li>
            <li>GC content and mapping rate are broadly stable between clusters (not significant in cluster comparison).</li>
          </ul>
        </div>
      </div>

      <div class="section">
        <h2>1. PCA Clustering (QC-driven)</h2>

        <div class="metric-grid">
          <div class="metric-card good">
            <h4>Group 1</h4>
            <div class="value">{"" if g1.get("pct") is None else f"{g1.get('pct'):.1f}%"} </div>
            <div class="label">{"" if g1.get("n") is None else f"{g1.get('n')} samples"}</div>
          </div>
          <div class="metric-card good">
            <h4>Group 2</h4>
            <div class="value">{"" if g2.get("pct") is None else f"{g2.get('pct'):.1f}%"} </div>
            <div class="label">{"" if g2.get("n") is None else f"{g2.get('n')} samples"}</div>
          </div>
          <div class="metric-card good">
            <h4>Primary Drivers</h4>
            <div class="value">Dup + Cov</div>
            <div class="label">Significant separation</div>
          </div>
        </div>

        <div class="figure-container">
          <div class="figure-title">Figure 1: Cordlife PCA clustering (QC metrics)</div>
          {img_tag("cord_pca", "Cordlife PCA Clustering")}
          <div class="figure-caption">
            Two clusters are visible. Based on the clustering report, <strong>duplication rate</strong> and <strong>mean coverage</strong> differ substantially between clusters.
          </div>
        </div>

        <div class="observation-box">
          <h4>Cluster-level comparison (from report)</h4>
          <table>
            <thead>
              <tr>
                <th>Metric</th>
                <th>Group 1 mean</th>
                <th>Group 2 mean</th>
                <th>Significance</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td><strong>duplication_rate(%)</strong></td>
                <td>{dup.get("group1","")}</td>
                <td>{dup.get("group2","")}</td>
                <td>{dup.get("sig","")}</td>
              </tr>
              <tr>
                <td><strong>mean_coverageData(X)</strong></td>
                <td>{cov.get("group1","")}</td>
                <td>{cov.get("group2","")}</td>
                <td>{cov.get("sig","")}</td>
              </tr>
              <tr>
                <td><strong>mapping_rate(%)</strong></td>
                <td>{mapr.get("group1","")}</td>
                <td>{mapr.get("group2","")}</td>
                <td>{mapr.get("sig","")}</td>
              </tr>
              <tr>
                <td><strong>GC_content(%)</strong></td>
                <td>{gc.get("group1","")}</td>
                <td>{gc.get("group2","")}</td>
                <td>{gc.get("sig","")}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>

      <div class="page-break"></div>

      <div class="section">
        <h2>2. Temporal Patterns</h2>
        <p>Time-series plots provide context on stability over batches and months.</p>

        <div class="figure-container">
          <div class="figure-title">Figure 2: Duplication rate timeline</div>
          {img_tag("cord_dup", "Cordlife Duplication Timeline")}
          <div class="figure-caption">Duplication rate across time/batches. Use this to spot batch-specific deviations.</div>
        </div>

        <div class="figure-container">
          <div class="figure-title">Figure 3: Coverage timeline</div>
          {img_tag("cord_cov", "Cordlife Coverage Timeline")}
          <div class="figure-caption">Coverage trends and outliers across time/batches.</div>
        </div>

        <div class="figure-container">
          <div class="figure-title">Figure 4: Batch effect overview (multiple QC/chr metrics)</div>
          {img_tag("cord_batch_over_time", "Cordlife Batch Effect Over Time")}
          <div class="figure-caption">Summary view for batch-level variation patterns across the analysis period.</div>
        </div>
      </div>

      <div class="section">
        <h2>3. Chromosome Proportion & Bias Context</h2>
        <p>This section summarizes chromosome distribution patterns and provides context for GC/duplication-driven bias.</p>

        <div class="figure-container">
          <div class="figure-title">Figure 5: Chromosome proportion distribution</div>
          {img_tag("cord_chr_dist", "Cordlife Chromosome Distribution")}
          <div class="figure-caption">Chromosome proportion distribution across samples.</div>
        </div>

        <div class="figure-container">
          <div class="figure-title">Figure 6: Readcount vs chromosome proportion</div>
          {img_tag("cord_readcount_vs_chrprop", "Readcount vs Chromosome Proportion")}
          <div class="figure-caption">Relationship between total readcount and chromosome proportion metrics.</div>
        </div>

        <div class="figure-container">
          <div class="figure-title">Figure 7: GC / duplication vs Chr bias relationship (context figure)</div>
          {img_tag("gc_chr_bias", "GC Chr Bias Relationship")}
          <div class="figure-caption">Context figure showing how GC and duplication can relate to chromosome proportion shifts.</div>
        </div>
      </div>

      <div class="section">
        <h2>4. Recommendations (Reference-building focused)</h2>
        <div class="recommendations">
          <h4>Practical selection principles for new references</h4>
          <ul>
            <li><strong>Exclude High Risk / No call</strong> samples from Normal reference pools (keep them only for positive controls where relevant).</li>
            <li><strong>Prefer cluster-stable region</strong>: prioritize samples close to the dominant cluster centroid to minimize reference variance.</li>
            <li><strong>Control duplication & coverage</strong>: because these drive cluster separation, avoid extreme duplication/coverage tails when building references.</li>
            <li><strong>Balance by time/batch</strong>: include multiple batches to reduce batch-specific bias (but avoid mixing incompatible process regimes if they exist).</li>
          </ul>
          <p style="margin-top:12px;">
            Implementation-wise, the repo already supports reference sample selection via <code>bin/scripts/utils/reference/create_reference.py</code> using
            a <code>reference_sample_list*.tsv</code> plus filters (mapping rate, SeqFF, duplication threshold, etc.).
          </p>
        </div>
      </div>

      <div class="section">
        <h2>Appendix: Data sources used</h2>
        <ul style="margin-left: 20px;">
          <li><code>analysis/batch_analysis_cordlife_complete/batch_effect_summary.txt</code></li>
          <li><code>analysis/batch_analysis_cordlife_complete/pca_clustering_report.txt</code></li>
          <li><code>analysis/batch_analysis_cordlife_complete/*.png</code> figures embedded as base64</li>
        </ul>
      </div>
    </div>

    <div class="footer">
      <p><strong>NIPT Sample Quality Assessment Report</strong></p>
      <p>Cordlife dataset | Analysis Period: {analysis_period}</p>
      <p>Report Generated: {report_date}</p>
      <p style="margin-top: 18px; font-size: 0.9em;">This report is generated from existing analysis outputs for review and reference-building support.</p>
    </div>
  </div>
</body>
</html>
"""
    return html_content


if __name__ == "__main__":
    print("Generating Cordlife Quality Assessment Report with figures...")
    html = create_report_with_figures()
    out_path = "Cordlife_Quality_Assessment_Report_Final.html"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\n✅ Report generated: {out_path}")
    print("\nOpen with: firefox Cordlife_Quality_Assessment_Report_Final.html")


