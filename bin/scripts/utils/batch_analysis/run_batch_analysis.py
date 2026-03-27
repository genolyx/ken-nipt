#!/usr/bin/env python3
"""
NIPT Batch Effect - 통합 분석 파이프라인
모든 batch effect 분석을 순차적으로 실행하고 종합 리포트 생성
"""

import sys
import os
import subprocess
import time
from datetime import datetime

def run_command(cmd, description):
    """
    Run a shell command and handle errors
    """
    print(f"\n{'='*80}")
    print(f"STEP: {description}")
    print(f"{'='*80}")
    print(f"Command: {cmd}\n")
    
    start_time = time.time()
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, 
                              capture_output=True, text=True)
        elapsed = time.time() - start_time
        
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        
        print(f"\n✓ Completed in {elapsed:.1f} seconds")
        return True
        
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"\n❌ Failed after {elapsed:.1f} seconds")
        print(f"Exit code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        return False

def create_html_report(output_dir, lab_name):
    """
    Create an HTML report with all figures
    """
    html_file = os.path.join(output_dir, 'batch_analysis_report.html')
    
    # List of expected figures
    figures = [
        ('01_batch_sample_counts.png', 'Sample Counts by Batch'),
        ('02_qc_by_batch.png', 'QC Metrics Distribution by Batch'),
        ('03_qc_trends.png', 'QC Metrics Trends Over Time'),
        ('04_ff_by_batch.png', 'Fetal Fraction by Batch'),
        ('05_chromosome_distributions.png', 'Chromosome Distribution by Batch'),
        ('06_pca_batch_and_lab.png', 'PCA Analysis - All Metrics'),
        ('07_pca_chromosome_only.png', 'PCA Analysis - Chromosome Proportions Only'),
        ('08_pca_qc_only.png', 'PCA Analysis - QC Metrics Only'),
        ('pca_clustering_groups.png', 'K-means Clustering (k=2)'),
        ('pca_groups_comparison.png', 'Statistical Comparison Between Groups'),
        ('duplication_rate_timeline.png', 'Duplication Rate Timeline'),
        ('duplication_rate_timeline_with_means.png', 'Duplication Rate Timeline with Batch Means'),
        ('coverage_timeline.png', 'Coverage Timeline'),
        ('coverage_timeline_with_means.png', 'Coverage Timeline with Batch Means'),
    ]
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Batch Effect Analysis - {lab_name}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 40px;
            border-left: 5px solid #3498db;
            padding-left: 15px;
        }}
        .figure-container {{
            background: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .figure-title {{
            font-size: 18px;
            font-weight: bold;
            color: #2c3e50;
            margin-bottom: 15px;
        }}
        img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}
        .info-box {{
            background: #e8f4f8;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .timestamp {{
            color: #7f8c8d;
            text-align: center;
            margin-top: 40px;
            font-size: 14px;
        }}
        .section {{
            margin-bottom: 60px;
        }}
        .report-link {{
            display: inline-block;
            background: #3498db;
            color: white;
            padding: 10px 20px;
            text-decoration: none;
            border-radius: 5px;
            margin: 10px 5px;
        }}
        .report-link:hover {{
            background: #2980b9;
        }}
    </style>
</head>
<body>
    <h1>Batch Effect Analysis Report - {lab_name}</h1>
    
    <div class="info-box">
        <p><strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Laboratory:</strong> {lab_name}</p>
        <p><strong>Output Directory:</strong> {output_dir}</p>
    </div>
    
    <div class="info-box">
        <strong>Text Reports:</strong><br>
        <a href="batch_effect_summary.txt" class="report-link">Batch Effect Summary</a>
        <a href="pca_clustering_report.txt" class="report-link">PCA Clustering Report</a>
        <a href="pca_group_analysis.txt" class="report-link">PCA Group Analysis</a>
        <a href="sample_list_with_clusters.tsv" class="report-link">Sample Data with Clusters (TSV)</a>
    </div>
"""
    
    # Section 1: Basic Batch Analysis
    html_content += """
    <div class="section">
        <h2>1. Basic Batch Analysis</h2>
"""
    
    for fig_file, fig_title in figures[:5]:
        fig_path = os.path.join(output_dir, fig_file)
        if os.path.exists(fig_path):
            html_content += f"""
        <div class="figure-container">
            <div class="figure-title">{fig_title}</div>
            <img src="{fig_file}" alt="{fig_title}">
        </div>
"""
    
    html_content += "    </div>\n"
    
    # Section 2: PCA Analysis
    html_content += """
    <div class="section">
        <h2>2. Principal Component Analysis (PCA)</h2>
        <p>PCA를 통해 샘플들의 전체적인 패턴과 배치 효과를 시각화합니다. 
        세 가지 PCA 분석을 수행했습니다: (1) 모든 메트릭, (2) 염색체 비율만, (3) QC 메트릭만</p>
"""
    
    for fig_file, fig_title in figures[5:8]:
        fig_path = os.path.join(output_dir, fig_file)
        if os.path.exists(fig_path):
            html_content += f"""
        <div class="figure-container">
            <div class="figure-title">{fig_title}</div>
            <img src="{fig_file}" alt="{fig_title}">
        </div>
"""
    
    html_content += "    </div>\n"
    
    # Section 3: Clustering Analysis
    html_content += """
    <div class="section">
        <h2>3. K-means Clustering Analysis</h2>
        <p>K-means 클러스터링(k=2)을 통해 샘플들을 두 그룹으로 나누고, 
        각 그룹의 특성을 통계적으로 비교했습니다.</p>
"""
    
    for fig_file, fig_title in figures[8:10]:
        fig_path = os.path.join(output_dir, fig_file)
        if os.path.exists(fig_path):
            html_content += f"""
        <div class="figure-container">
            <div class="figure-title">{fig_title}</div>
            <img src="{fig_file}" alt="{fig_title}">
        </div>
"""
    
    html_content += "    </div>\n"
    
    # Section 4: Timeline Analysis
    html_content += """
    <div class="section">
        <h2>4. Timeline Analysis</h2>
        <p>시간 흐름에 따른 Duplication Rate와 Coverage의 변화를 배치별로 시각화했습니다. 
        이를 통해 프로토콜 변경이나 시스템 변화를 감지할 수 있습니다.</p>
"""
    
    for fig_file, fig_title in figures[10:]:
        fig_path = os.path.join(output_dir, fig_file)
        if os.path.exists(fig_path):
            html_content += f"""
        <div class="figure-container">
            <div class="figure-title">{fig_title}</div>
            <img src="{fig_file}" alt="{fig_title}">
        </div>
"""
    
    html_content += "    </div>\n"
    
    # Footer
    html_content += f"""
    <div class="timestamp">
        Report generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    </div>
</body>
</html>
"""
    
    with open(html_file, 'w') as f:
        f.write(html_content)
    
    print(f"\n✓ HTML report generated: {html_file}")
    return html_file

def run_batch_analysis_pipeline(sample_list, analysis_dir, output_dir, lab_name):
    """
    Run complete batch effect analysis pipeline
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    print(f"\n{'#'*80}")
    print(f"# NIPT BATCH EFFECT ANALYSIS PIPELINE")
    print(f"#")
    print(f"# Sample List: {sample_list}")
    print(f"# Analysis Dir: {analysis_dir}")
    print(f"# Output Dir: {output_dir}")
    print(f"# Lab: {lab_name}")
    print(f"{'#'*80}\n")
    
    start_time = time.time()
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Basic batch effect analysis
    step1_cmd = (f"python3 {script_dir}/analyze_batch_effect_v2.py "
                f"--sample-list {sample_list} "
                f"--analysis-dir {analysis_dir} "
                f"--output-dir {output_dir} "
                f"--lab {lab_name}")
    
    if not run_command(step1_cmd, "Step 1: Basic Batch Effect Analysis"):
        print("\n❌ Pipeline failed at Step 1")
        return False
    
    # Step 2: PCA group analysis
    step2_cmd = (f"python3 {script_dir}/analyze_pca_groups.py "
                f"{sample_list} {output_dir}")
    
    if not run_command(step2_cmd, "Step 2: PCA Group Analysis"):
        print("\n⚠️ Warning: Step 2 failed, continuing...")
    
    # Step 3: PCA clustering analysis
    step3_cmd = (f"python3 {script_dir}/analyze_pca_clustering.py "
                f"{sample_list} {output_dir} {lab_name}")
    
    if not run_command(step3_cmd, "Step 3: K-means Clustering Analysis"):
        print("\n❌ Pipeline failed at Step 3")
        return False
    
    # Step 4: Timeline analysis
    step4_cmd = (f"python3 {script_dir}/plot_duplication_timeline.py "
                f"{output_dir} {lab_name}")
    
    if not run_command(step4_cmd, "Step 4: Timeline Analysis"):
        print("\n⚠️ Warning: Step 4 failed, continuing...")
    
    # Step 5: Generate HTML report
    print(f"\n{'='*80}")
    print("STEP: Generating HTML Report")
    print(f"{'='*80}\n")
    
    html_file = create_html_report(output_dir, lab_name)
    
    # Summary
    total_time = time.time() - start_time
    
    print(f"\n{'#'*80}")
    print(f"# PIPELINE COMPLETED SUCCESSFULLY")
    print(f"#")
    print(f"# Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    print(f"# Output directory: {output_dir}")
    print(f"# HTML report: {html_file}")
    print(f"{'#'*80}\n")
    
    print(f"\n📊 Open the HTML report to view all results:")
    print(f"   file://{os.path.abspath(html_file)}\n")
    
    return True

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 run_batch_analysis.py <sample_list.tsv> <analysis_dir> <output_dir> <lab_name>")
        print("\nExample:")
        print("  python3 run_batch_analysis.py \\")
        print("    data/refs/ucl/reference_make/reference_sample_list_UCL_all.tsv \\")
        print("    ~/ken-nipt/analysis \\")
        print("    analysis/batch_analysis_ucl \\")
        print("    UCL")
        sys.exit(1)
    
    sample_list = sys.argv[1]
    analysis_dir = sys.argv[2]
    output_dir = sys.argv[3]
    lab_name = sys.argv[4]
    
    # Validate inputs
    if not os.path.exists(sample_list):
        print(f"❌ Error: Sample list not found: {sample_list}")
        sys.exit(1)
    
    if not os.path.exists(analysis_dir):
        print(f"❌ Error: Analysis directory not found: {analysis_dir}")
        sys.exit(1)
    
    # Run pipeline
    success = run_batch_analysis_pipeline(sample_list, analysis_dir, output_dir, lab_name)
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
