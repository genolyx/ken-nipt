# Microdeletion Detection Performance Dashboard

Interactive web-based dashboard for analyzing microdeletion detection performance across multiple methods, fetal fractions, and deletion lengths.

## Overview

This dashboard provides a two-phase workflow:

1. **Phase 1: Threshold Optimization** - Interactive ROC analysis to determine optimal z-score thresholds and minimum detection lengths
2. **Phase 2: Performance Visualization** - Heatmaps and line plots showing sensitivity across FF and deletion length combinations

## Prerequisites

Install required Python packages:

```bash
pip install pandas numpy plotly dash
```

## Workflow

### Step 1: Prepare Compact Analysis Data

First, copy only the essential files (JSON, report.txt, aberrations.bed) to a compact directory:

#### Option A: Single Disease

```bash
python copy_analysis_results.py \
    --source /data/md_validation/1p36 \
    --dest /data/md_validation/analysis_result/1p36 \
    --zscore-tsv ~/ken-nipt/analysis/md_validation/zscore/zscore_extraction_1p36.tsv
```

#### Option B: All 8 Diseases at Once

```bash
# Copy from /data/md_validation (4 diseases)
python copy_analysis_results.py \
    --source-dirs /data/md_validation/1p36,/data/md_validation/2q33,/data/md_validation/CDC,/data/md_validation/DGS \
    --dest /data/md_validation/analysis_result \
    --zscore-dir ~/ken-nipt/analysis/md_validation/zscore

# Copy from ~/ken-nipt/analysis/md_validation (4 diseases)
python copy_analysis_results.py \
    --source-dirs ~/ken-nipt/analysis/md_validation/Jacobsen,~/ken-nipt/analysis/md_validation/PWS,~/ken-nipt/analysis/md_validation/WBS,~/ken-nipt/analysis/md_validation/WHS \
    --dest /data/md_validation/analysis_result \
    --zscore-dir ~/ken-nipt/analysis/md_validation/zscore
```

This will create:

```
/data/md_validation/analysis_result/
├── zscore_data/
│   ├── zscore_extraction_1p36.tsv
│   ├── zscore_extraction_2q33.tsv
│   ├── zscore_extraction_CDC.tsv
│   ├── zscore_extraction_DGS.tsv
│   ├── zscore_extraction_Jacobsen.tsv
│   ├── zscore_extraction_PWS.tsv
│   ├── zscore_extraction_WBS.tsv
│   └── zscore_extraction_WHS.tsv
├── 1p36/
│   ├── {sample1}/
│   │   ├── {sample1}.json
│   │   └── results/
│   │       ├── wc_orig_report.txt
│   │       ├── wc_fetus_report.txt
│   │       ├── wcx_orig_aberrations.bed
│   │       └── wcx_fetus_aberrations.bed
│   └── ... (2100 samples)
├── 2q33/
├── CDC/
└── ... (8 diseases total)
```

#### Option C: With Cleanup (Free Disk Space)

⚠️ **Warning**: This will DELETE BAM, NPZ, and plot files from the source directories!

```bash
# Dry-run first to see what would be deleted
python copy_analysis_results.py \
    --source /data/md_validation/1p36 \
    --dest /data/md_validation/analysis_result/1p36 \
    --zscore-tsv ~/ken-nipt/analysis/md_validation/zscore/zscore_extraction_1p36.tsv \
    --cleanup \
    --dry-run

# If satisfied, run without --dry-run
python copy_analysis_results.py \
    --source /data/md_validation/1p36 \
    --dest /data/md_validation/analysis_result/1p36 \
    --zscore-tsv ~/ken-nipt/analysis/md_validation/zscore/zscore_extraction_1p36.tsv \
    --cleanup
```

Expected space savings:
- **Before**: ~500 GB (with BAM/NPZ files)
- **After**: ~50 MB (JSON + result files only)
- **Space freed**: ~499.95 GB per disease

### Step 2: Launch Dashboard

```bash
python md_dashboard.py --port 8050
```

Then open your browser to: `http://localhost:8050`

### Step 3: Use the Dashboard

#### Phase 1: Threshold Optimization

1. **Load Data**:
   - Z-score TSV: `/data/md_validation/analysis_result/zscore_data/zscore_extraction_1p36.tsv`
   - Sample Directory: `/data/md_validation/analysis_result/1p36`
   - Click "Load Data"

2. **Configure ROC Analysis**:
   - Method: Select one (WC_orig, WC_fetus, WCX_orig, WCX_fetus)
   - Fetal Fraction: Select one (5%, 10%, 15%)
   - Min Detect Length: Select threshold (0, 0.5, 1, 2 Mb)
   - Click "Analyze ROC"

3. **Interpret Results**:
   - **ROC Curve**: Shows trade-off between TPR and FPR
   - **Sensitivity vs Threshold**: How sensitivity changes with z-score cutoff
   - **PPV vs Threshold**: How PPV changes with z-score cutoff
   - **Sensitivity-PPV Trade-off**: Interactive plot to select optimal point

4. **Select Optimal Threshold**:
   - Review suggested thresholds for PPV ≥ 90%, 80%, 70%
   - Hover over points to see exact values
   - Record the z-score threshold for use in Phase 2

#### Phase 2: Performance Visualization

1. **Configure Analysis**:
   - Method: Select detection method or aggregated mode (ORIG, FETUS, ANY)
   - Z-score Threshold: Enter value from Phase 1 (e.g., -5.0)
   - Min Length: Enter threshold in Mb (e.g., 1.0)
   - Click "Calculate Performance"

2. **Interpret Results**:
   - **Heatmap**: Sensitivity matrix (FF vs Deletion Length)
     - Green = High sensitivity
     - Red = Low sensitivity
   - **Line Plot**: Sensitivity trends across deletion lengths
     - Separate lines for each FF level
     - Shows how performance varies with deletion size

### Step 4: Compare Across Diseases

Repeat Steps 2-3 for each disease:

```bash
# 1p36 deletion syndrome
# TSV: /data/md_validation/analysis_result/zscore_data/zscore_extraction_1p36.tsv
# Dir: /data/md_validation/analysis_result/1p36

# 2q33.1 deletion syndrome
# TSV: /data/md_validation/analysis_result/zscore_data/zscore_extraction_2q33.tsv
# Dir: /data/md_validation/analysis_result/2q33

# ... etc for all 8 diseases
```

## Detection Methods

### Individual Methods
- **WC_orig**: Wisecondor with original fetus estimation
- **WC_fetus**: Wisecondor with optimized fetus mode
- **WCX_orig**: WisecondorX with original fetus estimation
- **WCX_fetus**: WisecondorX with optimized fetus mode

### Aggregated Methods (OR logic)
- **ORIG**: WC_orig OR WCX_orig (any orig method detects)
- **FETUS**: WC_fetus OR WCX_fetus (any fetus method detects)
- **ANY**: Detection by any of the 4 methods

## Performance Metrics

- **Sensitivity (Recall)**: TP / (TP + FN) - How many true deletions are detected
- **Specificity**: TN / (TN + FP) - How many true negatives are correctly identified
- **PPV (Precision)**: TP / (TP + FP) - Of detected deletions, how many are real
- **NPV**: TN / (TN + FN) - Of non-detected samples, how many truly have no deletion

## File Structure

### Input Files

#### Z-score Extraction TSV
Columns: `sample_name`, `disease`, `ff`, `deletion_length_mb`, `expected_deletion_chr`, `expected_deletion_start`, `expected_deletion_end`, `WC_orig_zscore`, `WC_fetus_zscore`, `WCX_orig_zscore`, `WCX_fetus_zscore`

Example:
```
sample_name	disease	ff	deletion_length_mb	expected_deletion_chr	expected_deletion_start	expected_deletion_end	WC_orig_zscore	WC_fetus_zscore	WCX_orig_zscore	WCX_fetus_zscore
1_1_williamsbeurensyndrome_FF5_15M_10Mb_F	WBS	5.0	10	7	72744455	82744455	-8.23	-7.45	-9.12	-8.88
```

#### Sample Directory Structure
```
{disease}/
└── {sample_name}/
    ├── {sample_name}.json          # Metadata (FF, deletion coordinates, gender)
    └── results/
        ├── wc_orig_report.txt      # Wisecondor orig detection results
        ├── wc_fetus_report.txt     # Wisecondor fetus detection results
        ├── wcx_orig_aberrations.bed # WisecondorX orig detection results
        └── wcx_fetus_aberrations.bed # WisecondorX fetus detection results
```

## Tips

### Recommended Thresholds

Based on clinical requirements:
- **High PPV (≥90%)**: For confirmatory testing, minimize false positives
- **Balanced (≥80% PPV, ≥80% Sensitivity)**: For screening
- **High Sensitivity (≥95%)**: For critical deletions, willing to accept more false positives

### Min Detection Length

- **0 Mb**: Detect all aberrations (may include noise)
- **0.5 Mb**: Filter very small aberrations
- **1 Mb**: Standard clinical threshold
- **2 Mb**: Very conservative, larger deletions only

### Performance Expectations

- **FF 5%**: Lower sensitivity, especially for small deletions (<3 Mb)
- **FF 10%**: Moderate sensitivity, good balance
- **FF 15%**: High sensitivity, even for smaller deletions

## Troubleshooting

### Dashboard not loading data
- Check file paths are correct (use absolute paths)
- Verify TSV file format matches expected structure
- Check sample directories contain required result files

### ROC analysis shows no data
- Verify the selected FF (5%, 10%, 15%) has samples in the dataset
- Check that result files exist for the selected method
- Ensure min_detect_length is appropriate for your data

### Performance plots show unexpected results
- Confirm z-score threshold is negative (deletions have negative z-scores)
- Verify min_length is in Mb, not bp
- Check that the selected method has results for all samples

## Advanced Usage

### Remote Access

To access dashboard from another machine:

```bash
# On server
python md_dashboard.py --host 0.0.0.0 --port 8050

# On client
# Open browser to: http://{server-ip}:8050
```

### Debug Mode

```bash
python md_dashboard.py --debug
```

This enables:
- Hot reload on code changes
- Detailed error messages
- Interactive debugger

## Contact

For issues or questions, contact the bioinformatics team.


