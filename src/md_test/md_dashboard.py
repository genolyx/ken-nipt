#!/usr/bin/env python3
"""
Microdeletion Detection Performance Dashboard

Interactive web-based dashboard for analyzing microdeletion detection performance.

Features:
- Phase 1: ROC Analysis & Threshold Optimization
- Phase 2: Performance Visualization with Selected Thresholds

Usage:
    python md_dashboard.py --port 8050
    
Then open browser: http://localhost:8050
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np

import dash
from dash import dcc, html, Input, Output, State, callback_context
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


# ============================================================================
# Data Processing Functions
# ============================================================================

def parse_wc_report_all_regions(report_file: Path) -> list:
    """Parse WC report.txt and extract all detected regions"""
    regions = []
    if not report_file.exists():
        return regions
    
    try:
        with open(report_file, 'r') as f:
            lines = f.readlines()
        
        in_test_section = False
        for line in lines:
            if '# test results:' in line.lower():
                in_test_section = True
                continue
            
            if in_test_section:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if 'z-score' in line.lower() and 'effect' in line.lower():
                    continue
                
                fields = line.split()
                if len(fields) < 4:
                    continue
                
                try:
                    zscore = float(fields[0])
                    location = fields[3]
                    if ':' in location and '-' in location:
                        chr_part, coord_part = location.split(':', 1)
                        start_str, end_str = coord_part.split('-', 1)
                        chr_name = str(chr_part).replace('chr', '')
                        start = int(start_str)
                        end = int(end_str)
                        regions.append({
                            'chr': chr_name,
                            'start': start,
                            'end': end,
                            'zscore': zscore,
                            'length': end - start
                        })
                except (ValueError, IndexError):
                    continue
    except Exception:
        pass
    
    return regions


def parse_wcx_bed_all_regions(bed_file: Path) -> list:
    """Parse WCX aberrations.bed and extract all detected regions"""
    regions = []
    if not bed_file.exists():
        return regions
    
    try:
        with open(bed_file, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 5:
                fields = line.split()
                if len(fields) < 5:
                    continue
            
            try:
                chr_name = str(fields[0]).replace('chr', '')
                start = int(fields[1])
                end = int(fields[2])
                zscore = float(fields[4])
                regions.append({
                    'chr': chr_name,
                    'start': start,
                    'end': end,
                    'zscore': zscore,
                    'length': end - start
                })
            except (ValueError, IndexError):
                continue
    except Exception:
        pass
    
    return regions


def check_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Check if two regions overlap"""
    return max(start1, start2) < min(end1, end2)


def load_data_from_dirs(
    zscore_tsv_path: str,
    sample_dirs: Dict[str, str]
) -> Tuple[pd.DataFrame, Dict]:
    """Load zscore TSV and sample directory info"""
    
    # Load zscore extraction data
    df = pd.read_csv(zscore_tsv_path, sep='\t')
    
    # Convert to numeric
    df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
    df['deletion_length_mb'] = pd.to_numeric(df['deletion_length_mb'], errors='coerce')
    
    for col in ['WC_orig_zscore', 'WC_fetus_zscore', 'WCX_orig_zscore', 'WCX_fetus_zscore']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Add disease column if not present
    if 'disease' not in df.columns:
        # Try to infer from sample_name or directory
        df['disease'] = 'Unknown'
    
    return df, sample_dirs


def calculate_roc_metrics(
    df: pd.DataFrame,
    sample_dir: Path,
    method: str,
    output_type: str,
    ff_value: float,
    min_detect_length: int,
    zcut_candidates: np.ndarray
) -> pd.DataFrame:
    """Calculate TPR, FPR, Sensitivity, Specificity, PPV for ROC curve"""
    
    ff_df = df[df['ff'] == ff_value].copy()
    
    if len(ff_df) == 0:
        return None
    
    logger.info(f"Calculating ROC for {method}_{output_type}, FF={ff_value}%, {len(ff_df)} samples")
    
    results = []
    
    for idx, zcut in enumerate(zcut_candidates):
        if idx % 20 == 0:
            logger.info(f"  Processing threshold {idx+1}/{len(zcut_candidates)}: z={zcut:.1f}")
        
        tp = 0
        fp = 0
        fn = 0
        tn = 0
        
        for _, row in ff_df.iterrows():
            sample_name = row['sample_name']
            sample_dir_path = sample_dir / sample_name
            
            expected_deletion = {
                'chr': str(row['expected_deletion_chr']).replace('chr', ''),
                'start': row['expected_deletion_start'],
                'end': row['expected_deletion_end']
            }
            
            # Get detected regions
            if method == 'WCX':
                bed_file = sample_dir_path / "results" / f"wcx_{output_type}_aberrations.bed"
                detected_regions = parse_wcx_bed_all_regions(bed_file)
            else:  # WC
                report_file = sample_dir_path / "results" / f"wc_{output_type}_report.txt"
                detected_regions = parse_wc_report_all_regions(report_file)
            
            # Filter by threshold
            filtered_regions = [r for r in detected_regions if r['zscore'] <= zcut]
            
            # Check target detection
            target_detected = False
            other_regions_detected = False
            
            for region in filtered_regions:
                if region['chr'] == expected_deletion['chr']:
                    if check_overlap(
                        expected_deletion['start'], expected_deletion['end'],
                        region['start'], region['end']
                    ):
                        target_detected = True
                        continue
                
                # Other region - check length
                if region['length'] >= min_detect_length:
                    other_regions_detected = True
            
            # Count TP, FP, FN, TN
            if target_detected:
                tp += 1
            else:
                fn += 1
            
            if other_regions_detected:
                fp += 1
            else:
                tn += 1
        
        # Calculate metrics
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 1.0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        fpr = 1 - specificity
        
        results.append({
            'threshold': zcut,
            'tp': tp,
            'fp': fp,
            'fn': fn,
            'tn': tn,
            'sensitivity': sensitivity,
            'specificity': specificity,
            'ppv': ppv,
            'tpr': sensitivity,
            'fpr': fpr
        })
    
    logger.info(f"  ROC calculation complete: {len(results)} thresholds evaluated")
    
    return pd.DataFrame(results)


# ============================================================================
# Dash App
# ============================================================================

# Initialize app
app = dash.Dash(__name__, suppress_callback_exceptions=True)

# Pre-scan files for initial load
initial_zscore_files = []
initial_sample_dirs = []

try:
    initial_zscore_files = scan_zscore_files()
    initial_sample_dirs = scan_sample_directories()
    logger.info(f"Found {len(initial_zscore_files)} z-score files and {len(initial_sample_dirs)} sample directories")
except Exception as e:
    logger.warning(f"Could not scan initial files: {e}")


def scan_zscore_files():
    """Scan for available z-score TSV files"""
    search_dirs = [
        Path('/home/ken/ken-nipt/analysis/md_validation/zscore'),
        Path('/data/md_validation/zscore_data')
    ]
    
    options = []
    for search_dir in search_dirs:
        if search_dir.exists():
            tsv_files = sorted(search_dir.glob('*.tsv'))
            for tsv_file in tsv_files:
                options.append({
                    'label': f'{tsv_file.stem} ({tsv_file.parent.name})',
                    'value': str(tsv_file)
                })
    
    return options


def scan_sample_directories():
    """Scan for available sample directories"""
    search_dirs = [
        Path('/data/md_validation/analysis_result'),
        Path('/home/ken/ken-nipt/analysis/md_validation')
    ]
    
    options = []
    for search_dir in search_dirs:
        if search_dir.exists():
            subdirs = sorted([d for d in search_dir.iterdir() if d.is_dir() and not d.name.startswith('.')])
            for subdir in subdirs:
                # Check if it's a disease directory (has sample subdirectories)
                sample_count = len([d for d in subdir.iterdir() if d.is_dir()])
                if sample_count > 0:
                    options.append({
                        'label': f'{subdir.name} ({sample_count} samples) - {search_dir.name}',
                        'value': str(subdir)
                    })
    
    return options


app.layout = html.Div([
    html.H1("Microdeletion Detection Performance Dashboard",
            style={'textAlign': 'center', 'color': '#2c3e50', 'marginBottom': 30}),
    
    # Configuration Panel
    html.Div([
        html.H3("Configuration", style={'color': '#34495e'}),
        
        html.Div([
            html.Div([
                html.Label("Z-score Extraction TSV:", style={'fontWeight': 'bold', 'marginBottom': 5}),
                dcc.Dropdown(
                    id='zscore-tsv-path',
                    options=initial_zscore_files,
                    placeholder='Select a z-score TSV file...',
                    style={'marginBottom': 10}
                ),
                html.Button('🔄 Refresh', id='refresh-files-btn', n_clicks=0,
                           style={'fontSize': '12px', 'padding': '5px 10px', 'marginBottom': 10})
            ], style={'marginBottom': 15}),
            
            html.Div([
                html.Label("Sample Directory:", style={'fontWeight': 'bold', 'marginBottom': 5}),
                dcc.Dropdown(
                    id='sample-dir-path',
                    options=initial_sample_dirs,
                    placeholder='Select a sample directory...',
                    style={'marginBottom': 10}
                ),
            ], style={'marginBottom': 15}),
        ]),
        
        html.Button('Load Data', id='load-data-btn', n_clicks=0,
                   style={'backgroundColor': '#3498db', 'color': 'white', 'padding': '10px 20px',
                         'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer', 'marginRight': 10}),
        
        html.Div(id='load-status', style={'marginTop': 10, 'color': '#27ae60'})
    ], style={'padding': 20, 'backgroundColor': '#ecf0f1', 'borderRadius': 10, 'marginBottom': 20}),
    
    # Tabs for two phases
    dcc.Tabs(id='main-tabs', value='phase1', children=[
        # Phase 1: Threshold Optimization
        dcc.Tab(label='Phase 1: Threshold Optimization', value='phase1', children=[
            html.Div([
                html.H3("ROC Analysis & Threshold Selection", style={'marginTop': 20}),
                
                html.Div([
                    html.Div([
                        html.Label("Method:"),
                        dcc.Dropdown(
                            id='method-selector',
                            options=[
                                {'label': 'WC orig', 'value': 'WC_orig'},
                                {'label': 'WC fetus', 'value': 'WC_fetus'},
                                {'label': 'WCX orig', 'value': 'WCX_orig'},
                                {'label': 'WCX fetus', 'value': 'WCX_fetus'}
                            ],
                            value='WCX_fetus'
                        )
                    ], style={'width': '30%', 'display': 'inline-block', 'marginRight': 20}),
                    
                    html.Div([
                        html.Label("Fetal Fraction (%):"),
                        dcc.Dropdown(
                            id='ff-selector',
                            options=[
                                {'label': '5%', 'value': 5},
                                {'label': '10%', 'value': 10},
                                {'label': '15%', 'value': 15}
                            ],
                            value=10
                        )
                    ], style={'width': '30%', 'display': 'inline-block', 'marginRight': 20}),
                    
                    html.Div([
                        html.Label("Min Detect Length (Mb):"),
                        dcc.Dropdown(
                            id='min-length-selector',
                            options=[
                                {'label': '0 Mb', 'value': 0},
                                {'label': '0.5 Mb', 'value': 500000},
                                {'label': '1 Mb', 'value': 1000000},
                                {'label': '2 Mb', 'value': 2000000}
                            ],
                            value=1000000
                        )
                    ], style={'width': '30%', 'display': 'inline-block'}),
                ], style={'marginBottom': 20}),
                
                html.Button('Analyze ROC', id='analyze-roc-btn', n_clicks=0,
                           style={'backgroundColor': '#e74c3c', 'color': 'white', 'padding': '10px 20px',
                                 'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer'}),
                
                dcc.Loading(
                    id='roc-loading',
                    type='default',
                    children=html.Div(id='roc-graph-container', style={'marginTop': 20})
                ),
                
                html.Div([
                    html.H4("Selected Optimal Thresholds", style={'marginTop': 30}),
                    html.Div(id='optimal-threshold-display')
                ])
            ])
        ]),
        
        # Phase 2: Performance Visualization
        dcc.Tab(label='Phase 2: Performance Analysis', value='phase2', children=[
            html.Div([
                html.H3("Performance Visualization", style={'marginTop': 20}),
                
                html.Div([
                    html.Div([
                        html.Label("Method:"),
                        dcc.Dropdown(
                            id='perf-method-selector',
                            options=[
                                {'label': 'WC orig', 'value': 'wc_orig'},
                                {'label': 'WC fetus', 'value': 'wc_fetus'},
                                {'label': 'WCX orig', 'value': 'wcx_orig'},
                                {'label': 'WCX fetus', 'value': 'wcx_fetus'},
                                {'label': 'ORIG (OR)', 'value': 'orig'},
                                {'label': 'FETUS (OR)', 'value': 'fetus'},
                                {'label': 'ANY (OR)', 'value': 'any'}
                            ],
                            value='wcx_fetus'
                        )
                    ], style={'width': '45%', 'display': 'inline-block', 'marginRight': 20}),
                    
                    html.Div([
                        html.Label("Z-score Threshold:"),
                        dcc.Input(
                            id='perf-zcut-input',
                            type='number',
                            value=-5.0,
                            step=0.1,
                            style={'width': '100%'}
                        )
                    ], style={'width': '20%', 'display': 'inline-block', 'marginRight': 20}),
                    
                    html.Div([
                        html.Label("Min Length (Mb):"),
                        dcc.Input(
                            id='perf-min-length-input',
                            type='number',
                            value=1.0,
                            step=0.1,
                            style={'width': '100%'}
                        )
                    ], style={'width': '20%', 'display': 'inline-block'}),
                ], style={'marginBottom': 20}),
                
                html.Button('Calculate Performance', id='calc-perf-btn', n_clicks=0,
                           style={'backgroundColor': '#27ae60', 'color': 'white', 'padding': '10px 20px',
                                 'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer'}),
                
                dcc.Loading(
                    id='perf-loading',
                    type='default',
                    children=[
                        html.Div(id='heatmap-container', style={'marginTop': 20}),
                        html.Div(id='lineplot-container', style={'marginTop': 20})
                    ]
                )
            ])
        ])
    ]),
    
    # Hidden data stores
    dcc.Store(id='loaded-data-store'),
    dcc.Store(id='optimal-thresholds-store')
], style={'padding': 30, 'fontFamily': 'Arial, sans-serif'})


# ============================================================================
# Callbacks
# ============================================================================

@app.callback(
    [Output('zscore-tsv-path', 'options'),
     Output('sample-dir-path', 'options')],
    [Input('refresh-files-btn', 'n_clicks')]
)
def refresh_file_lists(n_clicks):
    """Refresh available files and directories"""
    return scan_zscore_files(), scan_sample_directories()


@app.callback(
    [Output('load-status', 'children'),
     Output('loaded-data-store', 'data')],
    [Input('load-data-btn', 'n_clicks')],
    [State('zscore-tsv-path', 'value'),
     State('sample-dir-path', 'value')]
)
def load_data(n_clicks, zscore_path, sample_dir):
    """Load data callback"""
    if n_clicks == 0:
        return "", None
    
    if not zscore_path or not sample_dir:
        return html.Div("Please provide both TSV and sample directory paths", 
                       style={'color': '#e74c3c'}), None
    
    try:
        # Load TSV
        df = pd.read_csv(zscore_path, sep='\t')
        
        # Convert to records
        data_dict = {
            'df': df.to_dict('records'),
            'columns': df.columns.tolist(),
            'sample_dir': sample_dir,
            'n_samples': len(df)
        }
        
        return html.Div(
            f"✓ Loaded {len(df)} samples successfully!",
            style={'color': '#27ae60', 'fontWeight': 'bold'}
        ), data_dict
        
    except Exception as e:
        return html.Div(
            f"✗ Error loading data: {str(e)}",
            style={'color': '#e74c3c'}
        ), None


@app.callback(
    [Output('roc-graph-container', 'children'),
     Output('optimal-threshold-display', 'children')],
    [Input('analyze-roc-btn', 'n_clicks')],
    [State('loaded-data-store', 'data'),
     State('method-selector', 'value'),
     State('ff-selector', 'value'),
     State('min-length-selector', 'value')]
)
def update_roc_analysis(n_clicks, data_store, method, ff, min_length):
    """Update ROC analysis"""
    if n_clicks == 0 or data_store is None:
        return html.Div("Load data and click 'Analyze ROC'"), ""
    
    try:
        # Reconstruct dataframe
        df = pd.DataFrame(data_store['df'])
        sample_dir = Path(data_store['sample_dir'])
        
        # Convert to numeric
        df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
        
        # Parse method
        method_type, output_type = method.split('_')
        
        # Filter by FF
        ff_df = df[df['ff'] == ff].copy()
        
        if len(ff_df) == 0:
            return html.Div(f"No samples found for FF={ff}%"), ""
        
        # Generate z-score candidates
        zcut_candidates = np.arange(-50, -2, 0.5)
        
        # Calculate ROC metrics
        roc_df = calculate_roc_metrics(
            df, sample_dir, method_type, output_type, ff, min_length, zcut_candidates
        )
        
        if roc_df is None or len(roc_df) == 0:
            return html.Div("No ROC data available"), ""
        
        # Create ROC plots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'ROC Curve',
                'Sensitivity vs Threshold',
                'PPV vs Threshold',
                'Sensitivity-PPV Trade-off'
            )
        )
        
        # ROC Curve
        fig.add_trace(
            go.Scatter(
                x=roc_df['fpr'],
                y=roc_df['tpr'],
                mode='lines+markers',
                name='ROC',
                hovertemplate=(
                    'FPR=%{x:.3f}<br>' +
                    'TPR=%{y:.3f}<br>' +
                    'Threshold=%{text}<extra></extra>'
                ),
                text=[f'{t:.1f}' for t in roc_df['threshold']]
            ),
            row=1, col=1
        )
        
        # Diagonal line
        fig.add_trace(
            go.Scatter(
                x=[0, 1], y=[0, 1],
                mode='lines',
                line=dict(dash='dash', color='gray'),
                showlegend=False
            ),
            row=1, col=1
        )
        
        # Sensitivity vs Threshold
        fig.add_trace(
            go.Scatter(
                x=roc_df['threshold'],
                y=roc_df['sensitivity'],
                mode='lines+markers',
                name='Sensitivity',
                line=dict(color='green')
            ),
            row=1, col=2
        )
        
        # PPV vs Threshold
        fig.add_trace(
            go.Scatter(
                x=roc_df['threshold'],
                y=roc_df['ppv'],
                mode='lines+markers',
                name='PPV',
                line=dict(color='blue')
            ),
            row=2, col=1
        )
        
        # Sensitivity vs PPV
        fig.add_trace(
            go.Scatter(
                x=roc_df['ppv'],
                y=roc_df['sensitivity'],
                mode='lines+markers',
                name='Trade-off',
                hovertemplate=(
                    'PPV=%{x:.3f}<br>' +
                    'Sensitivity=%{y:.3f}<br>' +
                    'Threshold=%{text}<extra></extra>'
                ),
                text=[f'{t:.1f}' for t in roc_df['threshold']],
                line=dict(color='purple')
            ),
            row=2, col=2
        )
        
        fig.update_xaxes(title_text="FPR", row=1, col=1)
        fig.update_yaxes(title_text="TPR", row=1, col=1)
        fig.update_xaxes(title_text="Z-score Threshold", row=1, col=2)
        fig.update_yaxes(title_text="Sensitivity", row=1, col=2)
        fig.update_xaxes(title_text="Z-score Threshold", row=2, col=1)
        fig.update_yaxes(title_text="PPV", row=2, col=1)
        fig.update_xaxes(title_text="PPV", row=2, col=2)
        fig.update_yaxes(title_text="Sensitivity", row=2, col=2)
        
        fig.update_layout(
            height=800,
            title_text=f'{method}, FF={ff}%, MinLen={min_length/1e6:.1f}Mb',
            showlegend=True,
            template='plotly_white'
        )
        
        # Find optimal points
        ppv_90_mask = roc_df['ppv'] >= 0.9
        ppv_80_mask = roc_df['ppv'] >= 0.8
        ppv_70_mask = roc_df['ppv'] >= 0.7
        
        optimal_info = []
        
        for ppv_target, mask, color in [(90, ppv_90_mask, '#e74c3c'), 
                                          (80, ppv_80_mask, '#f39c12'), 
                                          (70, ppv_70_mask, '#f1c40f')]:
            if mask.any():
                best_idx = roc_df[mask]['sensitivity'].idxmax()
                best_row = roc_df.loc[best_idx]
                optimal_info.append(html.Div([
                    html.Span(f"PPV ≥ {ppv_target}%: ", style={'fontWeight': 'bold', 'color': color}),
                    html.Span(f"z = {best_row['threshold']:.1f}, "),
                    html.Span(f"Sensitivity = {best_row['sensitivity']:.2%}, "),
                    html.Span(f"PPV = {best_row['ppv']:.2%}")
                ], style={'marginBottom': 5}))
        
        if not optimal_info:
            optimal_info = [html.Div("No thresholds meet PPV ≥ 70%", style={'color': '#e74c3c'})]
        
        return dcc.Graph(figure=fig), html.Div(optimal_info)
        
    except Exception as e:
        return html.Div(f"Error: {str(e)}", style={'color': '#e74c3c'}), ""


def calculate_sensitivity_for_group(
    df: pd.DataFrame,
    sample_dir: Path,
    ff_value: float,
    length_mb: float,
    mode: str,
    zcut: float,
    min_length_bp: int
) -> Dict:
    """Calculate sensitivity for a specific FF and length combination"""
    
    # Filter samples
    group_df = df[(df['ff'] == ff_value) & (df['deletion_length_mb'] == length_mb)].copy()
    
    if len(group_df) == 0:
        return {'n': 0, 'tp': 0, 'fn': 0, 'sensitivity': 0.0}
    
    tp = 0
    fn = 0
    
    for _, row in group_df.iterrows():
        sample_name = row['sample_name']
        sample_dir_path = sample_dir / sample_name
        
        expected_deletion = {
            'chr': str(row['expected_deletion_chr']).replace('chr', ''),
            'start': row['expected_deletion_start'],
            'end': row['expected_deletion_end']
        }
        
        detected = False
        
        # Check detection based on mode
        if mode == 'wc_orig':
            detected = check_detection(
                sample_dir_path, 'WC', 'orig', expected_deletion, zcut, min_length_bp
            )
        elif mode == 'wc_fetus':
            detected = check_detection(
                sample_dir_path, 'WC', 'fetus', expected_deletion, zcut, min_length_bp
            )
        elif mode == 'wcx_orig':
            detected = check_detection(
                sample_dir_path, 'WCX', 'orig', expected_deletion, zcut, min_length_bp
            )
        elif mode == 'wcx_fetus':
            detected = check_detection(
                sample_dir_path, 'WCX', 'fetus', expected_deletion, zcut, min_length_bp
            )
        elif mode == 'orig':
            # OR of wc_orig and wcx_orig
            detected = (
                check_detection(sample_dir_path, 'WC', 'orig', expected_deletion, zcut, min_length_bp) or
                check_detection(sample_dir_path, 'WCX', 'orig', expected_deletion, zcut, min_length_bp)
            )
        elif mode == 'fetus':
            # OR of wc_fetus and wcx_fetus
            detected = (
                check_detection(sample_dir_path, 'WC', 'fetus', expected_deletion, zcut, min_length_bp) or
                check_detection(sample_dir_path, 'WCX', 'fetus', expected_deletion, zcut, min_length_bp)
            )
        elif mode == 'any':
            # OR of all 4 methods
            detected = (
                check_detection(sample_dir_path, 'WC', 'orig', expected_deletion, zcut, min_length_bp) or
                check_detection(sample_dir_path, 'WC', 'fetus', expected_deletion, zcut, min_length_bp) or
                check_detection(sample_dir_path, 'WCX', 'orig', expected_deletion, zcut, min_length_bp) or
                check_detection(sample_dir_path, 'WCX', 'fetus', expected_deletion, zcut, min_length_bp)
            )
        
        if detected:
            tp += 1
        else:
            fn += 1
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    
    return {
        'n': len(group_df),
        'tp': tp,
        'fn': fn,
        'sensitivity': sensitivity
    }


def check_detection(
    sample_dir: Path,
    method: str,
    output_type: str,
    expected_deletion: Dict,
    zcut: float,
    min_length_bp: int
) -> bool:
    """Check if target deletion is detected"""
    
    # Get detected regions
    if method == 'WCX':
        bed_file = sample_dir / "results" / f"wcx_{output_type}_aberrations.bed"
        detected_regions = parse_wcx_bed_all_regions(bed_file)
    else:  # WC
        report_file = sample_dir / "results" / f"wc_{output_type}_report.txt"
        detected_regions = parse_wc_report_all_regions(report_file)
    
    # Filter by threshold
    filtered_regions = [r for r in detected_regions if r['zscore'] <= zcut]
    
    # Check if target is detected
    for region in filtered_regions:
        if region['chr'] == expected_deletion['chr']:
            if check_overlap(
                expected_deletion['start'], expected_deletion['end'],
                region['start'], region['end']
            ):
                # Detected! But check if it meets minimum length
                if region['length'] >= min_length_bp:
                    return True
    
    return False


@app.callback(
    [Output('heatmap-container', 'children'),
     Output('lineplot-container', 'children')],
    [Input('calc-perf-btn', 'n_clicks')],
    [State('loaded-data-store', 'data'),
     State('perf-method-selector', 'value'),
     State('perf-zcut-input', 'value'),
     State('perf-min-length-input', 'value')]
)
def update_performance_plots(n_clicks, data_store, mode, zcut, min_length_mb):
    """Update performance visualization"""
    if n_clicks == 0 or data_store is None:
        return html.Div("Load data and click 'Calculate Performance'"), ""
    
    try:
        # Reconstruct dataframe
        df = pd.DataFrame(data_store['df'])
        sample_dir = Path(data_store['sample_dir'])
        
        df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
        df['deletion_length_mb'] = pd.to_numeric(df['deletion_length_mb'], errors='coerce')
        
        min_length_bp = int(min_length_mb * 1_000_000)
        
        # Calculate sensitivity for each (FF, length) combination
        heatmap_data = []
        ff_values = sorted(df['ff'].dropna().unique())
        length_values = sorted(df['deletion_length_mb'].dropna().unique())
        
        logger.info(f"Calculating performance for mode={mode}, zcut={zcut}, min_length={min_length_mb}Mb")
        
        for ff in ff_values:
            for length in length_values:
                result = calculate_sensitivity_for_group(
                    df, sample_dir, ff, length, mode, zcut, min_length_bp
                )
                
                if result['n'] > 0:
                    heatmap_data.append({
                        'FF': ff,
                        'Length': length,
                        'Sensitivity': result['sensitivity'],
                        'N': result['n'],
                        'TP': result['tp'],
                        'FN': result['fn']
                    })
        
        if len(heatmap_data) == 0:
            return (
                html.Div("No data available for the selected parameters", style={'color': '#e74c3c'}),
                ""
            )
        
        heatmap_df = pd.DataFrame(heatmap_data)
        
        # Create pivot table for heatmap
        pivot = heatmap_df.pivot(index='Length', columns='FF', values='Sensitivity')
        
        # Create heatmap
        heatmap_fig = go.Figure(data=go.Heatmap(
            z=pivot.values,
            x=pivot.columns,
            y=pivot.index,
            colorscale='RdYlGn',
            text=pivot.values,
            texttemplate='%{text:.1%}',
            textfont={"size": 12},
            colorbar=dict(title='Sensitivity'),
            hovertemplate=(
                'FF=%{x}%<br>' +
                'Length=%{y}Mb<br>' +
                'Sensitivity=%{z:.1%}<extra></extra>'
            )
        ))
        
        heatmap_fig.update_layout(
            title=f'Sensitivity Heatmap: {mode} (z≤{zcut}, MinLen≥{min_length_mb}Mb)',
            xaxis_title='Fetal Fraction (%)',
            yaxis_title='Deletion Length (Mb)',
            height=500,
            template='plotly_white'
        )
        
        # Create line plot
        lineplot_fig = go.Figure()
        
        colors = {'5.0': '#e74c3c', '10.0': '#f39c12', '15.0': '#27ae60'}
        
        for ff in sorted(heatmap_df['FF'].unique()):
            ff_data = heatmap_df[heatmap_df['FF'] == ff].sort_values('Length')
            color = colors.get(str(ff), '#3498db')
            
            lineplot_fig.add_trace(go.Scatter(
                x=ff_data['Length'],
                y=ff_data['Sensitivity'],
                mode='lines+markers',
                name=f'FF={int(ff)}%',
                marker=dict(size=10),
                line=dict(color=color, width=2),
                hovertemplate=(
                    f'FF={int(ff)}%<br>' +
                    'Length=%{x}Mb<br>' +
                    'Sensitivity=%{y:.1%}<br>' +
                    '<extra></extra>'
                )
            ))
        
        lineplot_fig.update_layout(
            title=f'Sensitivity vs Length: {mode} (z≤{zcut}, MinLen≥{min_length_mb}Mb)',
            xaxis_title='Deletion Length (Mb)',
            yaxis_title='Sensitivity',
            height=500,
            yaxis_range=[-0.05, 1.05],
            template='plotly_white'
        )
        
        return (
            dcc.Graph(figure=heatmap_fig),
            dcc.Graph(figure=lineplot_fig)
        )
        
    except Exception as e:
        import traceback
        error_msg = f"Error: {str(e)}\n{traceback.format_exc()}"
        logger.error(error_msg)
        return html.Div(error_msg, style={'color': '#e74c3c', 'whiteSpace': 'pre-wrap'}), ""


def main():
    parser = argparse.ArgumentParser(
        description="Microdeletion Detection Performance Dashboard"
    )
    parser.add_argument(
        '--port',
        type=int,
        default=8050,
        help='Port to run dashboard (default: 8050)'
    )
    parser.add_argument(
        '--host',
        type=str,
        default='0.0.0.0',
        help='Host to run dashboard (default: 0.0.0.0)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Run in debug mode'
    )
    
    args = parser.parse_args()
    
    logger.info("="*80)
    logger.info("Starting Microdeletion Detection Performance Dashboard")
    logger.info("="*80)
    logger.info(f"URL: http://localhost:{args.port}")
    logger.info("Press Ctrl+C to stop")
    logger.info("="*80)
    
    app.run(
        host=args.host,
        port=args.port,
        debug=args.debug
    )


if __name__ == '__main__':
    main()

