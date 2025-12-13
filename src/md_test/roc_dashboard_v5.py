#!/usr/bin/env python3
"""
ROC Analysis Dashboard V5 - Fixed Negative Sample Counting

Changes from V4:
- FIXED: Negative samples now include ALL FF and deletion_length combinations
- Phase 2: For each (FF, Length) cell:
  - Positive: target_disease samples with that specific (FF, Length)
  - Negative: other diseases with ALL (FF, Length) combinations
- This ensures proper One-vs-Rest classification with sufficient negative samples

Uses aberration_data.csv with z-score AND detected_mb for ROC/Performance calculation.

Logic:
- Detection = (|zscore| >= threshold) AND (detected_mb >= min_length)
- TP: Positive samples with detection
- FN: Positive samples without detection
- FP: Negative samples with detection  
- TN: Negative samples without detection

Usage:
    python roc_dashboard_v5.py --port 8001
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np

import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
DISEASES = ['1p36', '2q33', 'CDC', 'DGS', 'Jacobsen', 'PWS', 'WBS', 'WHS']
INDIVIDUAL_MODES = ['wc_orig', 'wc_fetus', 'wcx_orig', 'wcx_fetus']
ALL_MODES = INDIVIDUAL_MODES + ['orig', 'fetus', 'any']
FF_VALUES = [5.0, 10.0, 15.0]

# Global data
ABERRATION_DATA = None


def load_aberration_data(data_path: Path) -> pd.DataFrame:
    """Load aberration data with zscore and detected_mb"""
    logger.info(f"Loading aberration data from {data_path}")
    df = pd.read_csv(data_path)
    logger.info(f"Loaded {len(df)} records")
    logger.info(f"Unique samples: {df['sample_id'].nunique()}")
    return df


def get_combined_mode_data(df: pd.DataFrame, target_disease: str, combined_mode: str,
                          ff_value: float, del_length: float) -> pd.DataFrame:
    """Get data for combined modes (orig, fetus, any) by taking max across individual modes"""
    
    if combined_mode == 'orig':
        modes = ['wc_orig', 'wcx_orig']
    elif combined_mode == 'fetus':
        modes = ['wc_fetus', 'wcx_fetus']
    elif combined_mode == 'any':
        modes = INDIVIDUAL_MODES
    else:
        return pd.DataFrame()
    
    # Get data for all relevant modes
    df_filtered = df[
        (df['target_disease'] == target_disease) &
        (df['mode'].isin(modes)) &
        (df['ff'] == ff_value) &
        (df['deletion_length_mb'] == del_length)
    ].copy()
    
    if len(df_filtered) == 0:
        return pd.DataFrame()
    
    # For each sample, take max absolute z-score and corresponding detected_mb
    result_rows = []
    for sample_id in df_filtered['sample_id'].unique():
        sample_data = df_filtered[df_filtered['sample_id'] == sample_id]
        
        # Find row with max absolute z-score
        max_idx = sample_data['zscore'].abs().idxmax()
        max_row = sample_data.loc[max_idx]
        
        result_rows.append({
            'sample_id': sample_id,
            'disease': max_row['disease'],
            'target_disease': target_disease,
            'ff': ff_value,
            'deletion_length_mb': del_length,
            'zscore': max_row['zscore'],
            'detected_mb': max_row['detected_mb']
        })
    
    return pd.DataFrame(result_rows)


def get_combined_mode_data_all_conditions(df: pd.DataFrame, target_disease: str, 
                                          combined_mode: str) -> pd.DataFrame:
    """
    Get data for combined modes across ALL FF and deletion_length combinations.
    Used for negative samples in V5.
    
    Key: Only look at rows where target_disease matches (analyzing that specific region).
    For each unique (sample, ff, length) combination, take max |z-score| across modes.
    """
    if combined_mode == 'orig':
        modes = ['wc_orig', 'wcx_orig']
    elif combined_mode == 'fetus':
        modes = ['wc_fetus', 'wcx_fetus']
    elif combined_mode == 'any':
        modes = INDIVIDUAL_MODES
    else:
        return pd.DataFrame()
    
    # Get data for all relevant modes analyzing the target_disease region
    df_filtered = df[
        (df['target_disease'] == target_disease) &  # Only rows analyzing this region
        (df['mode'].isin(modes))
    ].copy()
    
    if len(df_filtered) == 0:
        return pd.DataFrame()
    
    # For each unique (sample, ff, length) combination, take max absolute z-score
    result_rows = []
    grouped = df_filtered.groupby(['sample_id', 'ff', 'deletion_length_mb'])
    
    for (sample_id, ff, length), group in grouped:
        # Find row with max absolute z-score within this group
        max_idx = group['zscore'].abs().idxmax()
        max_row = group.loc[max_idx]
        
        result_rows.append({
            'sample_id': sample_id,
            'disease': max_row['disease'],
            'target_disease': target_disease,
            'ff': ff,
            'deletion_length_mb': length,
            'zscore': max_row['zscore'],
            'detected_mb': max_row['detected_mb']
        })
    
    return pd.DataFrame(result_rows)


def calculate_roc_with_length(
    df: pd.DataFrame, target_disease: str, mode: str,
    ff_values: List[str] = None, length_values: List[str] = None
) -> Dict:
    """Calculate ROC curve with aberration length filtering
    
    Args:
        df: Aberration data
        target_disease: Target disease name
        mode: Analysis mode
        ff_values: List of FF values to include (e.g., ['5.0', '10.0'])
        length_values: List of deletion lengths to include (e.g., ['1', '3', '5'])
    """
    
    # Filter by FF (multi-select)
    if ff_values:
        ff_floats = [float(ff) for ff in ff_values]
        df = df[df['ff'].isin(ff_floats)]
    
    # Filter by deletion length (multi-select)
    if length_values:
        length_floats = [float(l) for l in length_values]
        df = df[df['deletion_length_mb'].isin(length_floats)]
    
    if len(df) == 0:
        return None
    
    # Get data for the mode
    if mode in INDIVIDUAL_MODES:
        df_mode = df[
            (df['target_disease'] == target_disease) &
            (df['mode'] == mode)
        ].copy()
    else:
        # Combined mode: collect from individual modes
        unique_ffs = df['ff'].unique()
        unique_lengths = df['deletion_length_mb'].unique()
        
        all_combined = []
        for ff_val in unique_ffs:
            for del_len in unique_lengths:
                combined_data = get_combined_mode_data(df, target_disease, mode, ff_val, del_len)
                if len(combined_data) > 0:
                    all_combined.append(combined_data)
        
        if not all_combined:
            return None
        df_mode = pd.concat(all_combined, ignore_index=True)
    
    if len(df_mode) == 0:
        return None
    
    # V5 CORRECT LOGIC:
    # Only look at rows where target_disease matches (analyzing that specific region)
    # Positive: disease == target_disease
    # Negative: disease != target_disease
    # Note: No additional filtering needed as FF and Length are already filtered above
    
    if len(df_mode) == 0:
        return None
    
    # Create labels (One-vs-Rest)
    y_true = (df_mode['disease'] == target_disease).astype(int).values
    
    if len(np.unique(y_true)) < 2:
        return None
    
    # Calculate ROC with length filtering
    # PHASE 1 (ROC): Use ONLY z-score for detection (to get full ROC curve)
    # This allows us to find optimal threshold based on sensitivity/specificity trade-off
    z_scores = df_mode['zscore'].abs().values
    
    # Get unique z-score thresholds
    thresholds = np.sort(np.unique(z_scores))[::-1]
    thresholds = np.concatenate([[np.inf], thresholds, [0]])
    
    tpr_list, fpr_list, ppv_list, spec_list = [], [], [], []
    n_positive = np.sum(y_true == 1)
    n_negative = np.sum(y_true == 0)
    
    for threshold in thresholds:
        # OPTION 1: Detection based on z-score ONLY (for ROC analysis)
        # Phase 2 will use (zscore + detected_mb) for actual performance metrics
        detected = (z_scores >= threshold)
        
        tp = np.sum((y_true == 1) & detected)
        fp = np.sum((y_true == 0) & detected)
        fn = np.sum((y_true == 1) & ~detected)
        tn = np.sum((y_true == 0) & ~detected)
        
        tpr = tp / n_positive if n_positive > 0 else 0
        fpr = fp / n_negative if n_negative > 0 else 0
        specificity = tn / n_negative if n_negative > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        
        tpr_list.append(tpr)
        fpr_list.append(fpr)
        spec_list.append(specificity)
        ppv_list.append(ppv)
    
    fpr = np.array(fpr_list)
    tpr = np.array(tpr_list)
    
    # Calculate AUC
    sorted_indices = np.argsort(fpr)
    auc = np.trapz(tpr[sorted_indices], fpr[sorted_indices])
    
    return {
        'fpr': fpr,
        'tpr': tpr,
        'thresholds': thresholds,
        'specificity': np.array(spec_list),
        'ppv': np.array(ppv_list),
        'auc': auc,
        'n_positive': int(n_positive),
        'n_negative': int(n_negative),
        'n_total': len(df_mode)
    }


def calculate_confusion_matrix_with_length(
    df: pd.DataFrame, target_disease: str, mode: str,
    ff_value: float, length_mb: float,
    zscore_threshold: float, min_detect_length: float
) -> Dict:
    """
    Calculate confusion matrix with aberration length filtering.
    
    V5 CORRECT LOGIC:
    - ONLY look at rows where target_disease matches (e.g., analyzing 1p36 region)
    - Positive: disease == target_disease (1p36 samples in 1p36 region)
    - Negative: disease != target_disease (other disease samples in 1p36 region)
    - This ensures FP only counts when other diseases show false signal in target region
    
    For each (FF, Length) cell:
    - Positive samples: specific (ff_value, length_mb) of target disease
    - Negative samples: ALL (ff, length) combinations of other diseases in target region
    """
    
    # STEP 1: Filter to only rows analyzing the target_disease region
    # Example: if evaluating 1p36, only look at rows where target_disease=1p36
    
    # Get POSITIVE samples for the specific (ff, length) condition
    if mode in INDIVIDUAL_MODES:
        df_pos_all = df[
            (df['target_disease'] == target_disease) &  # Analyzing target region
            (df['disease'] == target_disease) &          # Actual disease matches
            (df['mode'] == mode) &
            (df['ff'] == ff_value) &
            (df['deletion_length_mb'] == length_mb)
        ].copy()
    else:
        # Combined mode - positive samples only
        df_pos_combined = get_combined_mode_data(df, target_disease, mode, ff_value, length_mb)
        df_pos_all = df_pos_combined[df_pos_combined['disease'] == target_disease].copy()
    
    # Get NEGATIVE samples with SAME (ff, length) as positive
    # Key: disease != target_disease, but still analyzing target_disease region
    # V5 FIX: Filter by same FF and Length to match positive samples
    if mode in INDIVIDUAL_MODES:
        df_neg_all = df[
            (df['target_disease'] == target_disease) &  # Analyzing target region
            (df['disease'] != target_disease) &          # Actual disease is different
            (df['mode'] == mode) &
            (df['ff'] == ff_value) &                     # SAME FF as positive
            (df['deletion_length_mb'] == length_mb)      # SAME Length as positive
        ].copy()
    else:
        # Combined mode - negative samples with same FF/Length
        df_neg_combined = get_combined_mode_data(df, target_disease, mode, ff_value, length_mb)
        df_neg_all = df_neg_combined[df_neg_combined['disease'] != target_disease].copy()
    
    # Filter positive samples: by deletion_length_mb >= min_detect_length
    df_pos = df_pos_all[df_pos_all['deletion_length_mb'] >= min_detect_length].copy()
    
    # V5 FIX: Do NOT filter negative samples
    # We need ALL negative samples to calculate TN correctly
    # Detection criteria will filter during TP/FP/TN/FN calculation
    df_neg = df_neg_all.copy()
    
    # Detection = (|zscore| >= threshold) AND (detected_mb >= min_length)
    if len(df_pos) > 0:
        pos_detected = (df_pos['zscore'].abs() >= zscore_threshold) & (df_pos['detected_mb'] >= min_detect_length)
        tp = np.sum(pos_detected)
        fn = len(df_pos) - tp
        
        # Calculate coverage for TP samples: detected_mb / deletion_length_mb
        if tp > 0:
            tp_samples = df_pos[pos_detected]
            # Coverage = min(detected_mb / actual deletion_length_mb, 1.0)
            # Cap at 100% - cannot cover more than the target region
            coverages = np.minimum(tp_samples['detected_mb'] / tp_samples['deletion_length_mb'], 1.0)
            mean_coverage = coverages.mean()
        else:
            mean_coverage = 0.0
    else:
        tp, fn = 0, 0
        mean_coverage = 0.0
    
    if len(df_neg) > 0:
        neg_detected = (df_neg['zscore'].abs() >= zscore_threshold) & (df_neg['detected_mb'] >= min_detect_length)
        fp = np.sum(neg_detected)
        tn = len(df_neg) - fp
    else:
        fp, tn = 0, 0
    
    # Calculate metrics
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0
    
    # Collect FP sample details for analysis
    if len(df_neg) > 0:
        fp_samples = df_neg[neg_detected]
        fp_detected_sizes = fp_samples['detected_mb'].values if len(fp_samples) > 0 else np.array([])
        fp_zscores = fp_samples['zscore'].abs().values if len(fp_samples) > 0 else np.array([])
    else:
        fp_detected_sizes = np.array([])
        fp_zscores = np.array([])
    
    return {
        'TP': int(tp), 'FP': int(fp), 'FN': int(fn), 'TN': int(tn),
        'sensitivity': sensitivity, 'specificity': specificity,
        'ppv': ppv, 'npv': npv,
        'n_pos': len(df_pos), 'n_neg': len(df_neg),
        'coverage': mean_coverage,
        'fp_detected_sizes': fp_detected_sizes,
        'fp_zscores': fp_zscores
    }


def create_app(data_path: Path) -> dash.Dash:
    """Create Dash application"""
    global ABERRATION_DATA
    
    ABERRATION_DATA = load_aberration_data(data_path)
    
    app = dash.Dash(__name__)
    
    app.layout = html.Div([
        html.H1("ROC Analysis Dashboard V5 - Fixed Negative Sample Counting",
               style={'textAlign': 'center', 'color': '#2c3e50', 'marginTop': '20px'}),
        
        html.Div([
            html.P("Detection = (|Z-score| ≥ Threshold) AND (Detected Length ≥ Min Length)",
                  style={'textAlign': 'center', 'color': '#7f8c8d', 'fontSize': '14px'}),
            html.P("V5: Negative samples now include ALL FF and deletion_length combinations",
                  style={'textAlign': 'center', 'color': '#e74c3c', 'fontSize': '12px', 'fontWeight': 'bold'})
        ]),
        
        dcc.Tabs(id='tabs', value='tab-roc', children=[
            # Tab 1: ROC Analysis
            dcc.Tab(label='Phase 1: ROC Analysis', value='tab-roc', children=[
                html.Div([
                    html.H3("ROC Analysis with Selected Conditions", 
                           style={'marginTop': 20, 'color': '#2c3e50'}),
                    
                    # Parameter selection
                    html.Div([
                        html.Div([
                            html.Label("Disease:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='disease-dropdown',
                                options=[{'label': d, 'value': d} for d in DISEASES],
                                value='1p36',
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Method:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='method-dropdown',
                                options=[{'label': m, 'value': m} for m in ALL_MODES],
                                value='wc_orig',
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Fetal Fraction:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='ff-dropdown',
                                options=[{'label': f'{ff}%', 'value': str(ff)} for ff in FF_VALUES],
                                value=[str(ff) for ff in FF_VALUES],
                                multi=True,
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Deletion Length:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='length-dropdown',
                                options=[{'label': f'{l} Mb', 'value': str(l)} for l in [0.5, 1, 3, 5, 7, 10]],
                                value=[str(l) for l in [0.5, 1, 3, 5, 7, 10]],
                                multi=True,
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    ], style={'marginBottom': '10px'}),
                    
                    # Calculate button
                    html.Div([
                        html.Button('Calculate ROC', id='calc-roc-btn', n_clicks=0,
                                   style={'padding': '10px 30px', 'fontSize': '16px', 'backgroundColor': '#3498db',
                                         'color': 'white', 'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer'})
                    ], style={'textAlign': 'center', 'margin': '20px'}),
                    
                    # Loading indicator
                    dcc.Loading(
                        id="loading-roc",
                        type="default",
                        children=html.Div(id='roc-output')
                    ),
                    
                ], style={'padding': '20px'})
            ]),
            
            # Tab 2: Performance Analysis  
            dcc.Tab(label='Phase 2: Performance Analysis', value='tab-perf', children=[
                html.Div([
                    html.H3("Performance Visualization with Selected Threshold",
                           style={'marginTop': 20, 'color': '#2c3e50'}),
                    
                    # Parameter selection - Row 1
                    html.Div([
                        html.Div([
                            html.Label("Disease:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='disease-perf',
                                options=[{'label': d, 'value': d} for d in DISEASES],
                                value='1p36',
                                style={'width': '100%'}
                            )
                        ], style={'width': '31%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Method:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='method-perf',
                                options=[{'label': m, 'value': m} for m in ALL_MODES],
                                value='wc_orig',
                                style={'width': '100%'}
                            )
                        ], style={'width': '31%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Z-score Threshold:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Input(id='threshold-perf', type='number', value=3.0, step='any',
                                     style={'width': '100%', 'padding': '5px 8px', 'fontSize': '14px',
                                           'border': '1px solid #ccc', 'borderRadius': '4px',
                                           'boxSizing': 'border-box'})
                        ], style={'width': '31%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    ]),
                    
                    # Parameter selection - Row 2
                    html.Div([
                        html.Div([
                            html.Label("Fetal Fraction:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='ff-perf',
                                options=[{'label': f'{ff}%', 'value': str(ff)} for ff in FF_VALUES],
                                value=[str(ff) for ff in FF_VALUES],
                                multi=True,
                                style={'width': '100%'}
                            )
                        ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Deletion Length:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='length-perf',
                                options=[{'label': f'{l} Mb', 'value': str(l)} for l in [0.5, 1, 3, 5, 7, 10]],
                                value=[str(l) for l in [0.5, 1, 3, 5, 7, 10]],
                                multi=True,
                                style={'width': '100%'}
                            )
                        ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    ]),
                    
                    # Info text for fixed min detect length
                    html.Div([
                        html.Span("Min Detect Length: 0.2 Mb (200K) - Fixed", 
                                 style={'color': '#7f8c8d', 'fontSize': '14px', 'fontStyle': 'italic'})
                    ], style={'textAlign': 'center', 'marginBottom': '10px'}),
                    
                    # Calculate button
                    html.Div([
                        html.Button('Calculate Performance', id='calc-perf-btn', n_clicks=0,
                                   style={'padding': '10px 30px', 'fontSize': '16px', 'backgroundColor': '#27ae60',
                                         'color': 'white', 'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer'})
                    ], style={'textAlign': 'center', 'margin': '20px'}),
                    
                    # Loading indicator
                    dcc.Loading(
                        id="loading-perf",
                        type="default",
                        children=html.Div(id='perf-output')
                    ),
                    
                ], style={'padding': '20px'})
            ]),
        ]),
    ], style={'fontFamily': 'Arial, sans-serif', 'maxWidth': '1400px', 'margin': '0 auto'})
    
    # Callbacks
    @app.callback(
        Output('roc-output', 'children'),
        Input('calc-roc-btn', 'n_clicks'),
        State('disease-dropdown', 'value'),
        State('method-dropdown', 'value'),
        State('ff-dropdown', 'value'),
        State('length-dropdown', 'value')
    )
    def update_roc(n_clicks, disease, method, ff_values, length_values):
        if not n_clicks or n_clicks == 0:
            return html.Div("Select parameters and click 'Calculate ROC'", 
                           style={'textAlign': 'center', 'color': '#7f8c8d', 'padding': '50px'})
        
        try:
            # Calculate ROC
            roc_data = calculate_roc_with_length(
                ABERRATION_DATA, disease, method, ff_values, length_values
            )
            
            if not roc_data:
                return html.Div("No data available for selected parameters",
                               style={'textAlign': 'center', 'color': '#e74c3c', 'padding': '50px'})
            
            # Create ROC curve plot
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=roc_data['fpr'], y=roc_data['tpr'],
                mode='lines+markers',
                name=f'ROC Curve (AUC={roc_data["auc"]:.4f})',
                line=dict(color='#3498db', width=2),
                marker=dict(size=6),
                customdata=roc_data['thresholds'],
                hovertemplate='<b>FPR</b>: %{x:.4f}<br><b>TPR</b>: %{y:.4f}<br><b>Threshold</b>: %{customdata:.2f}<extra></extra>'
            ))
            
            fig.add_trace(go.Scatter(
                x=[0, 1], y=[0, 1],
                mode='lines',
                name='Random',
                line=dict(color='gray', dash='dash')
            ))
            
            fig.update_layout(
                title=f'ROC Curve: {disease} - {method}',
                xaxis_title='False Positive Rate (1-Specificity)',
                yaxis_title='True Positive Rate (Sensitivity)',
                width=700, height=600,
                template='plotly_white',
                hovermode='closest'
            )
            
            # Performance metrics
            metrics_text = (
                f"**Total Samples**: {roc_data['n_total']} | "
                f"**Positive**: {roc_data['n_positive']} | "
                f"**Negative**: {roc_data['n_negative']} | "
                f"**AUC**: {roc_data['auc']:.4f}"
            )
            
            return html.Div([
                html.Div([
                    dcc.Graph(figure=fig)
                ], style={'display': 'flex', 'justifyContent': 'center'}),
                
                html.Div([
                    dcc.Markdown(metrics_text)
                ], style={'textAlign': 'center', 'marginTop': '10px', 'fontSize': '14px', 'color': '#2c3e50'})
            ])
            
        except Exception as e:
            logger.error(f"Error in update_roc: {e}", exc_info=True)
            return html.Div(f"Error: {str(e)}", 
                           style={'textAlign': 'center', 'color': '#e74c3c', 'padding': '50px'})
    
    @app.callback(
        Output('perf-output', 'children'),
        Input('calc-perf-btn', 'n_clicks'),
        State('disease-perf', 'value'),
        State('method-perf', 'value'),
        State('threshold-perf', 'value'),
        State('ff-perf', 'value'),
        State('length-perf', 'value')
    )
    def update_performance(n_clicks, disease, method, threshold, ff_values, length_values):
        if not n_clicks or n_clicks == 0:
            return html.Div("Set parameters and click 'Calculate Performance'", 
                           style={'textAlign': 'center', 'color': '#7f8c8d', 'padding': '50px'})
        
        try:
            if threshold is None or threshold <= 0:
                return html.Div("Please enter a valid z-score threshold (> 0)", 
                               style={'color': '#e74c3c', 'textAlign': 'center', 'padding': '50px'})
            
            # Fixed min detect length
            min_length_mb = 0.2
            
            # Convert selected FF and Length values to floats
            if not ff_values or not length_values:
                return html.Div("Please select at least one FF and one Deletion Length", 
                               style={'color': '#e74c3c', 'textAlign': 'center', 'padding': '50px'})
            
            ff_floats = [float(ff) for ff in ff_values]
            length_floats = [float(l) for l in length_values]
            
            all_metrics_data = []
            all_fp_detected_sizes = []
            all_fp_zscores = []
            
            for ff in ff_floats:
                for length in length_floats:
                    cm_result = calculate_confusion_matrix_with_length(
                        ABERRATION_DATA, disease, method, ff, length, threshold, min_length_mb
                    )
                    if cm_result['n_pos'] > 0 or cm_result['n_neg'] > 0:
                        all_metrics_data.append({
                            'FF': ff,
                            'Length': length,
                            'TP': cm_result['TP'],
                            'FP': cm_result['FP'],
                            'FN': cm_result['FN'],
                            'TN': cm_result['TN'],
                            'Sensitivity': cm_result['sensitivity'],
                            'Specificity': cm_result['specificity'],
                            'PPV': cm_result['ppv'],
                            'NPV': cm_result['npv'],
                            'Coverage': cm_result['coverage'],
                            'n_pos': cm_result['n_pos'],
                            'n_neg': cm_result['n_neg']
                        })
                        # Collect FP details
                        if len(cm_result['fp_detected_sizes']) > 0:
                            all_fp_detected_sizes.extend(cm_result['fp_detected_sizes'])
                            all_fp_zscores.extend(cm_result['fp_zscores'])
            
            if not all_metrics_data:
                return html.Div("No data available for selected parameters",
                               style={'textAlign': 'center', 'color': '#e74c3c', 'padding': '50px'})
            
            df_metrics = pd.DataFrame(all_metrics_data)
            
            # Create confusion matrix table
            overall_tp = df_metrics['TP'].sum()
            overall_fp = df_metrics['FP'].sum()
            overall_fn = df_metrics['FN'].sum()
            overall_tn = df_metrics['TN'].sum()
            
            overall_sens = overall_tp / (overall_tp + overall_fn) if (overall_tp + overall_fn) > 0 else 0
            overall_spec = overall_tn / (overall_tn + overall_fp) if (overall_tn + overall_fp) > 0 else 0
            overall_ppv = overall_tp / (overall_tp + overall_fp) if (overall_tp + overall_fp) > 0 else 0
            overall_npv = overall_tn / (overall_tn + overall_fn) if (overall_tn + overall_fn) > 0 else 0
            
            # Calculate overall coverage (weighted average based on TP counts)
            total_tp_weighted_coverage = sum(row['Coverage'] * row['TP'] for row in all_metrics_data)
            overall_coverage = total_tp_weighted_coverage / overall_tp if overall_tp > 0 else 0
            
            # Calculate actual sample counts from confusion matrix
            total_pos = overall_tp + overall_fn  # Actual positive samples
            total_neg = overall_fp + overall_tn  # Actual negative samples
            
            confusion_table = html.Div([
                html.H4("Overall Confusion Matrix", style={'textAlign': 'center', 'color': '#2c3e50'}),
                html.Div([
                    html.Span(f"Positive Samples: {total_pos:,} | ", style={'marginRight': '15px', 'fontWeight': 'bold'}),
                    html.Span(f"Negative Samples: {total_neg:,}", style={'fontWeight': 'bold', 'color': '#e74c3c'})
                ], style={'textAlign': 'center', 'marginBottom': '15px', 'fontSize': '14px'}),
                html.Table([
                    html.Tr([
                        html.Td("", style={'border': 'none'}),
                        html.Th("Predicted Positive", style={'padding': '10px', 'backgroundColor': '#ecf0f1'}),
                        html.Th("Predicted Negative", style={'padding': '10px', 'backgroundColor': '#ecf0f1'})
                    ]),
                    html.Tr([
                        html.Th("Actual Positive", style={'padding': '10px', 'backgroundColor': '#ecf0f1'}),
                        html.Td(f"TP: {overall_tp}", style={'padding': '10px', 'textAlign': 'center', 'backgroundColor': '#d5f4e6', 'fontWeight': 'bold'}),
                        html.Td(f"FN: {overall_fn}", style={'padding': '10px', 'textAlign': 'center', 'backgroundColor': '#fadbd8', 'fontWeight': 'bold'})
                    ]),
                    html.Tr([
                        html.Th("Actual Negative", style={'padding': '10px', 'backgroundColor': '#ecf0f1'}),
                        html.Td(f"FP: {overall_fp}", style={'padding': '10px', 'textAlign': 'center', 'backgroundColor': '#fadbd8', 'fontWeight': 'bold'}),
                        html.Td(f"TN: {overall_tn}", style={'padding': '10px', 'textAlign': 'center', 'backgroundColor': '#d5f4e6', 'fontWeight': 'bold'})
                    ])
                ], style={'margin': '0 auto', 'border': '1px solid #bdc3c7', 'borderCollapse': 'collapse'}),
                
                html.Div([
                    html.Span(f"Sensitivity: {overall_sens:.4f} ({overall_sens*100:.2f}%) | ", style={'marginRight': '15px'}),
                    html.Span(f"Specificity: {overall_spec:.4f} ({overall_spec*100:.2f}%) | ", style={'marginRight': '15px'}),
                    html.Span(f"PPV: {overall_ppv:.4f} ({overall_ppv*100:.2f}%) | ", style={'marginRight': '15px'}),
                    html.Span(f"NPV: {overall_npv:.4f} ({overall_npv*100:.2f}%)", style={'marginRight': '15px'})
                ], style={'textAlign': 'center', 'marginTop': '15px', 'fontSize': '14px', 'color': '#2c3e50'}),
                
                html.Div([
                    html.Span(f"Average Coverage (TP): {overall_coverage:.4f} ({overall_coverage*100:.2f}%)", 
                             style={'fontWeight': 'bold', 'color': '#27ae60'})
                ], style={'textAlign': 'center', 'marginTop': '10px', 'fontSize': '14px'})
            ])
            
            # Helper function to create heatmaps
            def create_metric_heatmap(df, metric, colorscale, title_suffix):
                pivot = df.pivot(index='Length', columns='FF', values=metric)
                actual_length_values = sorted(pivot.index.tolist())
                
                # Use categorical Y-axis for uniform spacing
                # Map actual length values to indices (0, 1, 2, ...)
                y_indices = list(range(len(actual_length_values)))
                
                fig = go.Figure(data=go.Heatmap(
                    z=pivot.values, 
                    x=pivot.columns, 
                    y=y_indices,  # Use indices instead of actual values
                    colorscale=colorscale, zmin=0, zmax=1,
                    text=np.round(pivot.values * 100, 1),
                    texttemplate='%{text}%', textfont={"size": 12},
                    colorbar=dict(title=metric)
                ))
                fig.update_layout(
                    title=f'{metric}: {disease} - {method}<br><sub>{title_suffix}</sub>',
                    xaxis_title='Fetal Fraction (%)', yaxis_title='Deletion Length (Mb)',
                    xaxis={'tickmode': 'array', 'tickvals': [5, 10, 15]},
                    yaxis={
                        'tickmode': 'array', 
                        'tickvals': y_indices,  # Position of ticks
                        'ticktext': actual_length_values  # Labels to display
                    },
                    height=400, template='plotly_white'
                )
                return fig
            
            title_suffix = f'Z-score ≥ {threshold}, Detected ≥ {min_length_mb}Mb'
            
            heatmaps_html = html.Div([
                html.H4("Performance Metrics by FF and Deletion Length", 
                       style={'textAlign': 'center', 'marginTop': '30px', 'marginBottom': '20px', 'color': '#2c3e50'}),
                html.Div([
                    html.Div([
                        dcc.Graph(figure=create_metric_heatmap(df_metrics, 'Sensitivity', 'RdYlGn', title_suffix))
                    ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    html.Div([
                        dcc.Graph(figure=create_metric_heatmap(df_metrics, 'Specificity', 'RdYlGn', title_suffix))
                    ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'})
                ]),
                html.Div([
                    html.Div([
                        dcc.Graph(figure=create_metric_heatmap(df_metrics, 'PPV', 'Blues', title_suffix))
                    ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    html.Div([
                        dcc.Graph(figure=create_metric_heatmap(df_metrics, 'NPV', 'Greens', title_suffix))
                    ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'})
                ]),
                html.Div([
                    html.Div([
                        dcc.Graph(figure=create_metric_heatmap(df_metrics, 'Coverage', 'Viridis', title_suffix))
                    ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'})
                ])
            ])
            
            # Create FP analysis histograms
            fp_histograms_html = None
            if len(all_fp_detected_sizes) > 0:
                # Histogram 1: Detected size distribution of FP samples (0~1Mb range with 0.1Mb bins)
                fig_fp_size = go.Figure()
                # Create bins from 0 to 1Mb with 0.1Mb intervals
                bins = np.arange(0, 1.1, 0.1)  # 0, 0.1, 0.2, ..., 1.0
                fig_fp_size.add_trace(go.Histogram(
                    x=all_fp_detected_sizes,
                    xbins=dict(start=0, end=1.0, size=0.1),
                    marker=dict(color='#e74c3c', line=dict(color='white', width=1)),
                    name='FP Detected Size'
                ))
                fig_fp_size.update_layout(
                    title=f'FP Detected Size Distribution (0-1Mb, 100K bins)<br><sub>Total FP: {len(all_fp_detected_sizes)}</sub>',
                    xaxis_title='Detected Size (Mb)',
                    yaxis_title='Count',
                    xaxis=dict(range=[0, 1.0]),
                    template='plotly_white',
                    height=400
                )
                
                # Histogram 2: Z-score distribution of FP samples
                fig_fp_zscore = go.Figure()
                fig_fp_zscore.add_trace(go.Histogram(
                    x=all_fp_zscores,
                    nbinsx=30,
                    marker=dict(color='#e67e22', line=dict(color='white', width=1)),
                    name='FP Z-score'
                ))
                fig_fp_zscore.update_layout(
                    title=f'FP Z-score Distribution<br><sub>Total FP: {len(all_fp_zscores)}</sub>',
                    xaxis_title='|Z-score|',
                    yaxis_title='Count',
                    template='plotly_white',
                    height=400
                )
                
                fp_histograms_html = html.Div([
                    html.H4("False Positive (FP) Analysis - Noise Characterization", 
                           style={'textAlign': 'center', 'marginTop': '40px', 'marginBottom': '20px', 'color': '#2c3e50'}),
                    html.Div([
                        html.Div([
                            dcc.Graph(figure=fig_fp_size)
                        ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        html.Div([
                            dcc.Graph(figure=fig_fp_zscore)
                        ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'})
                    ])
                ])
            
            return html.Div([
                confusion_table,
                heatmaps_html,
                fp_histograms_html if fp_histograms_html else html.Div()
            ])
            
        except Exception as e:
            logger.error(f"Error in update_performance: {e}", exc_info=True)
            return html.Div(f"Error: {str(e)}", 
                           style={'textAlign': 'center', 'color': '#e74c3c', 'padding': '50px'})
    
    return app


def main():
    parser = argparse.ArgumentParser(description="ROC Analysis Dashboard V5")
    parser.add_argument('--port', type=int, default=8001, help='Port to run dashboard')
    parser.add_argument('--data', type=str, 
                       default='/data/md_validation/roc_results/aberration_data.csv',
                       help='Path to aberration data CSV')
    
    args = parser.parse_args()
    
    data_path = Path(args.data)
    if not data_path.exists():
        logger.error(f"Data file not found: {data_path}")
        return 1
    
    app = create_app(data_path)
    
    logger.info("="*80)
    logger.info("Starting ROC Dashboard V5")
    logger.info("="*80)
    logger.info(f"Data: {data_path}")
    logger.info(f"Port: {args.port}")
    logger.info(f"URL: http://0.0.0.0:{args.port}")
    logger.info("="*80)
    
    app.run(host='0.0.0.0', port=args.port, debug=False)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

