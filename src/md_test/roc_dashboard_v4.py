#!/usr/bin/env python3
"""
ROC Analysis Dashboard V4 - With Aberration Length Filtering

Uses aberration_data.csv with z-score AND detected_mb for ROC/Performance calculation.

Logic:
- Detection = (|zscore| >= threshold) AND (detected_mb >= min_length)
- TP: Positive samples with detection
- FN: Positive samples without detection
- FP: Negative samples with detection  
- TN: Negative samples without detection

Usage:
    python roc_dashboard_v4.py --port 8001
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


def calculate_roc_with_length(
    df: pd.DataFrame, target_disease: str, mode: str,
    ff_filter: str = 'All', length_filter: str = 'All',
    custom_length: float = None, min_detect_length: float = 0.5
) -> Dict:
    """Calculate ROC curve with aberration length filtering"""
    
    # Filter by FF
    if ff_filter != 'All':
        df = df[df['ff'] == float(ff_filter)]
    
    # Filter by deletion length
    if length_filter == 'Custom' and custom_length is not None:
        df = df[df['deletion_length_mb'] >= custom_length]
    elif length_filter != 'All':
        df = df[df['deletion_length_mb'] >= int(length_filter)]
    
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
    
    # NEW LOGIC (Updated): Filter both positive and negative samples
    # Positive: by deletion_length_mb >= min_detect_length (designed deletion size)
    # Negative: by detected_mb >= min_detect_length (actual detected aberration size)
    df_pos = df_mode[df_mode['disease'] == target_disease].copy()
    df_neg = df_mode[df_mode['disease'] != target_disease].copy()
    
    # Filter positive samples: only include those with deletion_length_mb >= min_detect_length
    df_pos_filtered = df_pos[df_pos['deletion_length_mb'] >= min_detect_length].copy()
    
    # Filter negative samples: only include those with detected_mb >= min_detect_length
    # This ensures we only count "clinically significant" aberrations as FP
    df_neg_filtered = df_neg[df_neg['detected_mb'] >= min_detect_length].copy()
    
    # Combine filtered positive + filtered negative samples
    df_mode = pd.concat([df_pos_filtered, df_neg_filtered], ignore_index=True)
    
    if len(df_mode) == 0:
        return None
    
    # Create labels (One-vs-Rest)
    y_true = (df_mode['disease'] == target_disease).astype(int).values
    
    if len(np.unique(y_true)) < 2:
        return None
    
    # Calculate ROC with length filtering
    # For each threshold, count TP, FP, FN, TN considering BOTH zscore and detected_mb
    z_scores = df_mode['zscore'].abs().values
    detected_lengths = df_mode['detected_mb'].values
    
    # Get unique z-score thresholds
    thresholds = np.sort(np.unique(z_scores))[::-1]
    thresholds = np.concatenate([[np.inf], thresholds, [0]])
    
    tpr_list, fpr_list, ppv_list, spec_list = [], [], [], []
    n_positive = np.sum(y_true == 1)
    n_negative = np.sum(y_true == 0)
    
    for threshold in thresholds:
        # Detection = (|zscore| >= threshold) AND (detected_mb >= min_length)
        detected = (z_scores >= threshold) & (detected_lengths >= min_detect_length)
        
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
    """Calculate confusion matrix with aberration length filtering"""
    
    # Get data for the mode
    if mode in INDIVIDUAL_MODES:
        df_mode = df[
            (df['target_disease'] == target_disease) &
            (df['mode'] == mode) &
            (df['ff'] == ff_value) &
            (df['deletion_length_mb'] == length_mb)
        ].copy()
    else:
        # Combined mode
        df_mode = get_combined_mode_data(df, target_disease, mode, ff_value, length_mb)
    
    if len(df_mode) == 0:
        return {
            'TP': 0, 'FP': 0, 'FN': 0, 'TN': 0,
            'sensitivity': 0.0, 'specificity': 0.0,
            'ppv': 0.0, 'npv': 0.0,
            'n_pos': 0, 'n_neg': 0
        }
    
    # Separate positive and negative samples
    df_pos_all = df_mode[df_mode['disease'] == target_disease].copy()
    df_neg_all = df_mode[df_mode['disease'] != target_disease].copy()
    
    # NEW LOGIC (Updated): Filter both positive and negative samples
    # Positive: by deletion_length_mb >= min_detect_length
    df_pos = df_pos_all[df_pos_all['deletion_length_mb'] >= min_detect_length].copy()
    
    # Negative: by detected_mb >= min_detect_length
    # Only count "clinically significant" aberrations as FP
    df_neg = df_neg_all[df_neg_all['detected_mb'] >= min_detect_length].copy()
    
    # Detection = (|zscore| >= threshold) AND (detected_mb >= min_length)
    if len(df_pos) > 0:
        pos_detected = (df_pos['zscore'].abs() >= zscore_threshold) & (df_pos['detected_mb'] >= min_detect_length)
        tp = np.sum(pos_detected)
        fn = len(df_pos) - tp
    else:
        tp, fn = 0, 0
    
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
    
    return {
        'TP': int(tp), 'FP': int(fp), 'FN': int(fn), 'TN': int(tn),
        'sensitivity': sensitivity, 'specificity': specificity,
        'ppv': ppv, 'npv': npv,
        'n_pos': len(df_pos), 'n_neg': len(df_neg)
    }


def create_app(data_path: Path) -> dash.Dash:
    """Create Dash application"""
    global ABERRATION_DATA
    
    ABERRATION_DATA = load_aberration_data(data_path)
    
    app = dash.Dash(__name__)
    
    app.layout = html.Div([
        html.H1("ROC Analysis Dashboard V4 - With Aberration Length Filtering",
               style={'textAlign': 'center', 'color': '#2c3e50', 'marginTop': '20px'}),
        
        html.Div([
            html.P("Detection = (|Z-score| ≥ Threshold) AND (Detected Length ≥ Min Length)",
                  style={'textAlign': 'center', 'color': '#7f8c8d', 'fontSize': '14px'})
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
                                options=[{'label': 'All', 'value': 'All'}] + 
                                       [{'label': f'{ff}%', 'value': str(ff)} for ff in FF_VALUES],
                                value='All',
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Deletion Length:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='length-dropdown',
                                options=[{'label': 'All', 'value': 'All'}] + 
                                       [{'label': f'{l} Mb', 'value': str(l)} for l in [1, 3, 5, 7, 10]] +
                                       [{'label': 'Custom', 'value': 'Custom'}],
                                value='All',
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    ], style={'marginBottom': '10px'}),
                    
                    # Custom length and Min Detect Length inputs
                    html.Div([
                        html.Div([
                            html.Label("Custom Deletion Length (Mb):", style={'marginRight': '10px'}),
                            dcc.Input(id='custom-length-input', type='number', value=2, min=0.1, max=15, step=0.1, 
                                     style={'width': '100px'})
                        ], id='custom-length-div', style={'display': 'none', 'padding': '10px', 'textAlign': 'center'}),
                        
                        html.Div([
                            html.Label("Min Detect Length (Mb):", style={'marginRight': '10px', 'fontWeight': 'bold'}),
                            dcc.Input(id='min-detect-length-input', type='number', value=0.5, min=0.1, max=10, step=0.1, 
                                     style={'width': '100px'}),
                            html.Span(" (Aberration must be ≥ this length)", 
                                     style={'marginLeft': '10px', 'color': '#7f8c8d', 'fontSize': '12px'})
                        ], style={'padding': '10px', 'textAlign': 'center'}),
                    ]),
                    
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
                    
                    # Parameter selection
                    html.Div([
                        html.Div([
                            html.Label("Disease:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='disease-perf',
                                options=[{'label': d, 'value': d} for d in DISEASES],
                                value='1p36',
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Method:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Dropdown(
                                id='method-perf',
                                options=[{'label': m, 'value': m} for m in ALL_MODES],
                                value='wc_orig',
                                style={'width': '100%'}
                            )
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Z-score Threshold:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Input(id='threshold-perf', type='number', value=3.0, step='any',
                                     style={'width': '100%', 'padding': '5px 8px', 'fontSize': '14px',
                                           'border': '1px solid #ccc', 'borderRadius': '4px',
                                           'boxSizing': 'border-box'})
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                        
                        html.Div([
                            html.Label("Min Detect Length (Mb):", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                            dcc.Input(id='minlen-perf', type='number', value=1.0, step='any',
                                     style={'width': '100%', 'padding': '5px 8px', 'fontSize': '14px',
                                           'border': '1px solid #ccc', 'borderRadius': '4px',
                                           'boxSizing': 'border-box'})
                        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    ]),
                    
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
        Output('custom-length-div', 'style'),
        Input('length-dropdown', 'value')
    )
    def toggle_custom_length(length_value):
        if length_value == 'Custom':
            return {'display': 'block', 'padding': '10px', 'textAlign': 'center'}
        return {'display': 'none'}
    
    @app.callback(
        Output('roc-output', 'children'),
        Input('calc-roc-btn', 'n_clicks'),
        State('disease-dropdown', 'value'),
        State('method-dropdown', 'value'),
        State('ff-dropdown', 'value'),
        State('length-dropdown', 'value'),
        State('custom-length-input', 'value'),
        State('min-detect-length-input', 'value')
    )
    def update_roc(n_clicks, disease, method, ff, length, custom_length, min_detect_length):
        if not n_clicks or n_clicks == 0:
            return html.Div("Select parameters and click 'Calculate ROC'", 
                           style={'textAlign': 'center', 'color': '#7f8c8d', 'padding': '50px'})
        
        try:
            # Calculate ROC
            roc_data = calculate_roc_with_length(
                ABERRATION_DATA, disease, method, ff, length, custom_length, min_detect_length
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
                title=f'ROC Curve: {disease} - {method}<br><sub>Min Detect Length: {min_detect_length}Mb</sub>',
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
        State('minlen-perf', 'value')
    )
    def update_performance(n_clicks, disease, method, threshold, min_length_mb):
        if not n_clicks or n_clicks == 0:
            return html.Div("Set parameters and click 'Calculate Performance'", 
                           style={'textAlign': 'center', 'color': '#7f8c8d', 'padding': '50px'})
        
        try:
            if threshold is None or threshold <= 0:
                return html.Div("Please enter a valid z-score threshold (> 0)", 
                               style={'color': '#e74c3c', 'textAlign': 'center', 'padding': '50px'})
            
            if min_length_mb is None or min_length_mb <= 0:
                return html.Div("Please enter a valid min detect length (> 0)", 
                               style={'color': '#e74c3c', 'textAlign': 'center', 'padding': '50px'})
            
            # Fixed FF and Length values for heatmap axes
            ff_values = [5.0, 10.0, 15.0]
            all_length_values = [0.5, 1, 3, 5, 7, 10]
            
            # Filter length values based on min_length_mb
            if min_length_mb is not None and min_length_mb > 0:
                length_values = [l for l in all_length_values if l >= min_length_mb]
            else:
                length_values = all_length_values
            
            if not length_values:
                return html.Div("No length values satisfy the min detect length criterion", 
                               style={'color': '#e74c3c', 'textAlign': 'center', 'padding': '50px'})
            
            all_metrics_data = []
            
            for ff in ff_values:
                for length in length_values:
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
                            'NPV': cm_result['npv']
                        })
            
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
            
            confusion_table = html.Div([
                html.H4("Overall Confusion Matrix", style={'textAlign': 'center', 'color': '#2c3e50'}),
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
                ], style={'textAlign': 'center', 'marginTop': '15px', 'fontSize': '14px', 'color': '#2c3e50'})
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
                ])
            ])
            
            return html.Div([
                confusion_table,
                heatmaps_html
            ])
            
        except Exception as e:
            logger.error(f"Error in update_performance: {e}", exc_info=True)
            return html.Div(f"Error: {str(e)}", 
                           style={'textAlign': 'center', 'color': '#e74c3c', 'padding': '50px'})
    
    return app


def main():
    parser = argparse.ArgumentParser(description="ROC Analysis Dashboard V4")
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
    logger.info("Starting ROC Dashboard V4")
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


