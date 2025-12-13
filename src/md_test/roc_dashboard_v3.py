#!/usr/bin/env python3
"""
ROC Analysis Dashboard V3 - With Interactive Features and Performance Analysis

Interactive dashboard with:
- Phase 1: ROC Analysis with clickable curve (shows threshold and metrics)
- Phase 2: Performance Visualization (Heatmap & Line plots)

Usage:
    python roc_dashboard_v3.py --port 8001
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
from plotly.subplots import make_subplots

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
DISEASES = ['1p36', '2q33', 'CDC', 'DGS', 'Jacobsen', 'PWS', 'WBS', 'WHS']
INDIVIDUAL_MODES = ['wc_orig', 'wc_fetus', 'wcx_orig', 'wcx_fetus']
GROUP_MODES = ['orig', 'fetus', 'any']
ALL_MODES = INDIVIDUAL_MODES + GROUP_MODES
FF_VALUES = [5.0, 10.0, 15.0]
LENGTH_VALUES = [1, 3, 5, 7, 10]

DISEASE_CHROMOSOMES = {
    '1p36': '1', '2q33': '2', 'CDC': '5', 'DGS': '22',
    'Jacobsen': '11', 'PWS': '15', 'WBS': '7', 'WHS': '4'
}

# Global data
COLLECTED_DATA = None


def load_collected_data(data_path: Path) -> pd.DataFrame:
    """Load pre-calculated collected data"""
    logger.info(f"Loading data from {data_path}")
    df = pd.read_csv(data_path)
    logger.info(f"Loaded {len(df)} records")
    return df


def calculate_roc_curve_manual(y_true: np.ndarray, y_score: np.ndarray) -> Tuple:
    """Calculate ROC curve manually"""
    thresholds = np.sort(np.unique(y_score))[::-1]
    thresholds = np.concatenate([[np.inf], thresholds, [-np.inf]])
    
    tpr_list, fpr_list = [], []
    n_positive = np.sum(y_true == 1)
    n_negative = np.sum(y_true == 0)
    
    for threshold in thresholds:
        y_pred = (y_score >= threshold).astype(int)
        tp = np.sum((y_true == 1) & (y_pred == 1))
        fp = np.sum((y_true == 0) & (y_pred == 1))
        
        tpr = tp / n_positive if n_positive > 0 else 0
        fpr = fp / n_negative if n_negative > 0 else 0
        
        tpr_list.append(tpr)
        fpr_list.append(fpr)
    
    return np.array(fpr_list), np.array(tpr_list), thresholds


def calculate_auc(fpr: np.ndarray, tpr: np.ndarray) -> float:
    """Calculate AUC"""
    sorted_indices = np.argsort(fpr)
    return np.trapz(tpr[sorted_indices], fpr[sorted_indices])


def filter_and_calculate_roc(
    df: pd.DataFrame, target_disease: str, mode: str,
    ff_filter: str = 'All', length_filter: str = 'All', custom_length: float = None
) -> Dict:
    """Filter data and calculate ROC curve"""
    df_target = df[df['target_disease'] == target_disease].copy()
    
    if len(df_target) == 0:
        return None
    
    # Apply filters
    if ff_filter != 'All':
        df_target = df_target[df_target['ff'] == float(ff_filter)]
    
    if length_filter == 'Custom' and custom_length is not None:
        df_target = df_target[df_target['deletion_length_mb'] >= custom_length]
    elif length_filter != 'All':
        df_target = df_target[df_target['deletion_length_mb'] >= int(length_filter)]
    
    if len(df_target) == 0:
        return None
    
    # Labels
    y_true = (df_target['disease'] == target_disease).astype(int).values
    
    # Scores
    if mode in INDIVIDUAL_MODES:
        if mode not in df_target.columns:
            return None
        y_score = df_target[mode].values
    elif mode == 'orig':
        cols = [c for c in ['wc_orig', 'wcx_orig'] if c in df_target.columns]
        y_score = df_target[cols].max(axis=1).values if cols else None
    elif mode == 'fetus':
        cols = [c for c in ['wc_fetus', 'wcx_fetus'] if c in df_target.columns]
        y_score = df_target[cols].max(axis=1).values if cols else None
    elif mode == 'any':
        cols = [c for c in INDIVIDUAL_MODES if c in df_target.columns]
        y_score = df_target[cols].max(axis=1).values if cols else None
    else:
        return None
    
    if y_score is None or len(np.unique(y_true)) < 2:
        return None
    
    fpr, tpr, thresholds = calculate_roc_curve_manual(y_true, y_score)
    roc_auc = calculate_auc(fpr, tpr)
    
    # Calculate Specificity, PPV for each threshold
    specificities = 1 - fpr
    ppv_list = []
    
    for threshold in thresholds:
        y_pred = (y_score >= threshold).astype(int)
        tp = np.sum((y_true == 1) & (y_pred == 1))
        fp = np.sum((y_true == 0) & (y_pred == 1))
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        ppv_list.append(ppv)
    
    return {
        'fpr': fpr, 'tpr': tpr, 'thresholds': thresholds,
        'specificity': specificities, 'ppv': np.array(ppv_list),
        'auc': roc_auc,
        'n_positive': int(y_true.sum()),
        'n_negative': int((1 - y_true).sum()),
        'n_total': len(df_target)
    }


def calculate_sensitivity_for_group(
    df: pd.DataFrame, target_disease: str, mode: str,
    ff_value: float, length_mb: float, threshold: float
) -> Dict:
    """Calculate sensitivity for specific FF and length"""
    # Filter samples
    group_df = df[
        (df['target_disease'] == target_disease) &
        (df['disease'] == target_disease) &
        (df['ff'] == ff_value) &
        (df['deletion_length_mb'] == length_mb)
    ].copy()
    
    if len(group_df) == 0:
        return {'n': 0, 'detected': 0, 'sensitivity': 0.0}
    
    # Get scores
    if mode in INDIVIDUAL_MODES:
        if mode not in group_df.columns:
            return {'n': 0, 'detected': 0, 'sensitivity': 0.0}
        scores = group_df[mode].fillna(0).values
    elif mode == 'orig':
        cols = [c for c in ['wc_orig', 'wcx_orig'] if c in group_df.columns]
        if not cols:
            return {'n': 0, 'detected': 0, 'sensitivity': 0.0}
        scores = group_df[cols].fillna(0).max(axis=1).values
    elif mode == 'fetus':
        cols = [c for c in ['wc_fetus', 'wcx_fetus'] if c in group_df.columns]
        if not cols:
            return {'n': 0, 'detected': 0, 'sensitivity': 0.0}
        scores = group_df[cols].fillna(0).max(axis=1).values
    elif mode == 'any':
        cols = [c for c in INDIVIDUAL_MODES if c in group_df.columns]
        if not cols:
            return {'n': 0, 'detected': 0, 'sensitivity': 0.0}
        scores = group_df[cols].fillna(0).max(axis=1).values
    else:
        return {'n': 0, 'detected': 0, 'sensitivity': 0.0}
    
    # Ensure no None/NaN values
    scores = np.nan_to_num(scores, nan=0.0)
    
    detected = np.sum(scores >= threshold)
    total = len(scores)
    
    return {
        'n': total,
        'detected': detected,
        'sensitivity': detected / total if total > 0 else 0.0
    }


def get_scores_for_samples(df, target_disease, mode, ff_value, length_mb, is_positive=True):
    """Get z-scores for positive or negative samples"""
    if is_positive:
        # Positive samples: the target disease
        group_df = df[
            (df['target_disease'] == target_disease) &
            (df['disease'] == target_disease) &
            (df['ff'] == ff_value) &
            (df['deletion_length_mb'] == length_mb)
        ].copy()
    else:
        # Negative samples: other diseases (One-vs-Rest)
        group_df = df[
            (df['target_disease'] == target_disease) &
            (df['disease'] != target_disease) &
            (df['ff'] == ff_value) &
            (df['deletion_length_mb'] == length_mb)
        ].copy()
    
    if len(group_df) == 0:
        return np.array([])
    
    # Get scores
    if mode in INDIVIDUAL_MODES:
        if mode not in group_df.columns:
            return np.array([])
        scores = group_df[mode].fillna(0).values
    elif mode == 'orig':
        cols = [c for c in ['wc_orig', 'wcx_orig'] if c in group_df.columns]
        if not cols:
            return np.array([])
        scores = group_df[cols].fillna(0).abs().max(axis=1).values
    elif mode == 'fetus':
        cols = [c for c in ['wc_fetus', 'wcx_fetus'] if c in group_df.columns]
        if not cols:
            return np.array([])
        scores = group_df[cols].fillna(0).abs().max(axis=1).values
    elif mode == 'any':
        cols = [c for c in INDIVIDUAL_MODES if c in group_df.columns]
        if not cols:
            return np.array([])
        scores = group_df[cols].fillna(0).abs().max(axis=1).values
    else:
        return np.array([])
    
    # Ensure no None/NaN values
    scores = np.nan_to_num(scores, nan=0.0)
    return np.abs(scores)


def calculate_confusion_matrix(df, target_disease, mode, ff_value, length_mb, threshold):
    """Calculate confusion matrix and all metrics for a specific group"""
    # Get positive samples (the target disease)
    pos_scores = get_scores_for_samples(df, target_disease, mode, ff_value, length_mb, is_positive=True)
    # Get negative samples (other diseases)
    neg_scores = get_scores_for_samples(df, target_disease, mode, ff_value, length_mb, is_positive=False)
    
    if len(pos_scores) == 0 and len(neg_scores) == 0:
        return {
            'TP': 0, 'FP': 0, 'FN': 0, 'TN': 0,
            'sensitivity': 0.0, 'specificity': 0.0,
            'ppv': 0.0, 'npv': 0.0,
            'n_pos': 0, 'n_neg': 0
        }
    
    # Calculate confusion matrix
    TP = np.sum(pos_scores >= threshold) if len(pos_scores) > 0 else 0
    FN = np.sum(pos_scores < threshold) if len(pos_scores) > 0 else 0
    FP = np.sum(neg_scores >= threshold) if len(neg_scores) > 0 else 0
    TN = np.sum(neg_scores < threshold) if len(neg_scores) > 0 else 0
    
    # Calculate metrics
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0.0
    ppv = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    npv = TN / (TN + FN) if (TN + FN) > 0 else 0.0
    
    return {
        'TP': int(TP), 'FP': int(FP), 'FN': int(FN), 'TN': int(TN),
        'sensitivity': sensitivity, 'specificity': specificity,
        'ppv': ppv, 'npv': npv,
        'n_pos': len(pos_scores), 'n_neg': len(neg_scores)
    }


# Initialize Dash app
app = dash.Dash(__name__, suppress_callback_exceptions=True)

app.layout = html.Div([
    html.Div([
        html.H1("ROC Analysis Dashboard V3", 
                style={'textAlign': 'center', 'color': '#2c3e50', 'marginBottom': 20}),
    
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
                            options=[{'label': f'{d} (chr{DISEASE_CHROMOSOMES[d]})', 'value': d} for d in DISEASES],
                            value='1p36', clearable=False
                        )
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    
                    html.Div([
                        html.Label("Method:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                        dcc.Dropdown(
                            id='method-dropdown',
                            options=[
                                {'label': 'WC Original', 'value': 'wc_orig'},
                                {'label': 'WC Fetus', 'value': 'wc_fetus'},
                                {'label': 'WCX Original', 'value': 'wcx_orig'},
                                {'label': 'WCX Fetus', 'value': 'wcx_fetus'},
                                {'label': 'Original (WC+WCX)', 'value': 'orig'},
                                {'label': 'Fetus (WC+WCX)', 'value': 'fetus'},
                                {'label': 'Any (All 4)', 'value': 'any'},
                            ],
                            value='wc_orig', clearable=False
                        )
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    
                    html.Div([
                        html.Label("Fetal Fraction:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                        dcc.Dropdown(
                            id='ff-dropdown',
                            options=[
                                {'label': 'All FF', 'value': 'All'},
                                {'label': 'FF = 5%', 'value': '5'},
                                {'label': 'FF = 10%', 'value': '10'},
                                {'label': 'FF = 15%', 'value': '15'}
                            ],
                            value='All', clearable=False
                        )
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    
                    html.Div([
                        html.Label("Min Deletion Length:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                        dcc.Dropdown(
                            id='length-dropdown',
                            options=[
                                {'label': 'All Lengths', 'value': 'All'},
                                {'label': '≥ 1 Mb', 'value': '1'},
                                {'label': '≥ 3 Mb', 'value': '3'},
                                {'label': '≥ 5 Mb', 'value': '5'},
                                {'label': '≥ 7 Mb', 'value': '7'},
                                {'label': '≥ 10 Mb', 'value': '10'},
                                {'label': 'Custom', 'value': 'Custom'}
                            ],
                            value='All', clearable=False
                        )
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                ], style={'backgroundColor': '#ecf0f1', 'padding': '20px', 'borderRadius': '10px'}),
                
                html.Div([
                    html.Label("Custom Min Length (Mb):", style={'fontWeight': 'bold', 'marginRight': '10px'}),
                    dcc.Input(id='custom-length-input', type='number', value=2, min=0.1, max=15, step=0.1, 
                             style={'width': '100px'})
                ], id='custom-length-div', style={'display': 'none', 'padding': '10px', 'textAlign': 'center'}),
                
                html.Div([
                    html.Button('Calculate ROC', id='calculate-button', n_clicks=0,
                               style={'fontSize': '16px', 'padding': '10px 30px', 'backgroundColor': '#3498db',
                                     'color': 'white', 'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer'})
                ], style={'textAlign': 'center', 'padding': '20px'}),
                
                # ROC plot and stats with single loading indicator
                dcc.Loading(id="loading-roc", type="default", children=[
                    html.Div([
                        html.Div(id='roc-output', style={'padding': '20px', 'display': 'flex', 'justifyContent': 'center'}),
                        html.Div(id='stats-output', style={'padding': '20px'})
                    ])
                ])
            ])
        ]),
        
        # Tab 2: Performance Analysis
        dcc.Tab(label='Phase 2: Performance Analysis', value='tab-perf', children=[
            html.Div([
                html.H3("Performance Visualization with Selected Threshold", 
                       style={'marginTop': 20, 'color': '#2c3e50'}),
                
                html.Div([
                    html.Div([
                        html.Label("Disease:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                        dcc.Dropdown(id='disease-perf', options=[{'label': d, 'value': d} for d in DISEASES],
                                   value='1p36', clearable=False)
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    
                    html.Div([
                        html.Label("Method:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                        dcc.Dropdown(
                            id='method-perf',
                            options=[
                                {'label': 'WC Original', 'value': 'wc_orig'},
                                {'label': 'WC Fetus', 'value': 'wc_fetus'},
                                {'label': 'WCX Original', 'value': 'wcx_orig'},
                                {'label': 'WCX Fetus', 'value': 'wcx_fetus'},
                                {'label': 'Original (WC+WCX)', 'value': 'orig'},
                                {'label': 'Fetus (WC+WCX)', 'value': 'fetus'},
                                {'label': 'Any (All 4)', 'value': 'any'},
                            ],
                            value='wc_orig', clearable=False
                        )
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    
                    html.Div([
                        html.Label("Z-score Threshold:", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                        dcc.Input(id='threshold-perf', type='number', value=5.0, 
                                 step='any',
                                 style={'width': '100%', 'padding': '5px 8px', 'fontSize': '14px',
                                       'border': '1px solid #ccc', 'borderRadius': '4px',
                                       'boxSizing': 'border-box', 'height': '38px'})
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                    
                    html.Div([
                        html.Label("Min Detect Length (Mb):", style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
                        dcc.Input(id='minlen-perf', type='number', value=1.0, 
                                 step='any',
                                 style={'width': '100%', 'padding': '5px 8px', 'fontSize': '14px',
                                       'border': '1px solid #ccc', 'borderRadius': '4px',
                                       'boxSizing': 'border-box', 'height': '38px'})
                    ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                ], style={'backgroundColor': '#ecf0f1', 'padding': '20px', 'borderRadius': '10px'}),
                
                html.Div([
                    html.Button('Calculate Performance', id='calc-perf-btn', n_clicks=0,
                               style={'fontSize': '16px', 'padding': '10px 30px', 'backgroundColor': '#27ae60',
                                     'color': 'white', 'border': 'none', 'borderRadius': '5px', 'cursor': 'pointer'})
                ], style={'textAlign': 'center', 'padding': '20px'}),
                
                html.Div([
                    html.P("ℹ️ Displays Confusion Matrix and performance metrics at fixed FF (5%, 10%, 15%) and Length (1, 3, 5, 7, 10 Mb) grid", 
                          style={'textAlign': 'center', 'color': '#7f8c8d', 'fontSize': '13px', 'fontStyle': 'italic'})
                ], style={'padding': '0 20px'}),
                
                dcc.Loading(id="loading-perf", type="default", children=[
                    html.Div(id='confusion-matrix-container', style={'marginTop': 20}),
                    html.Div(id='metrics-heatmaps-container', style={'marginTop': 20})
                ])
            ], style={'padding': '20px'})
        ])
    ])
    ], style={
        'maxWidth': '1400px',
        'margin': '0 auto',
        'backgroundColor': 'white',
        'padding': '20px',
        'boxShadow': '0 0 10px rgba(0,0,0,0.1)'
    })
], style={'backgroundColor': '#ecf0f1', 'minHeight': '100vh', 'padding': '20px'})


@app.callback(
    Output('custom-length-div', 'style'),
    Input('length-dropdown', 'value')
)
def toggle_custom_length(length_value):
    if length_value == 'Custom':
        return {'display': 'block', 'padding': '10px', 'textAlign': 'center', 'backgroundColor': '#f8f9fa'}
    return {'display': 'none'}


@app.callback(
    [Output('roc-output', 'children'),
     Output('stats-output', 'children')],
    Input('calculate-button', 'n_clicks'),
    [State('disease-dropdown', 'value'),
     State('method-dropdown', 'value'),
     State('ff-dropdown', 'value'),
     State('length-dropdown', 'value'),
     State('custom-length-input', 'value')]
)
def update_roc(n_clicks, disease, method, ff_filter, length_filter, custom_length):
    if not n_clicks or n_clicks == 0:
        return (html.Div("Select parameters and click 'Calculate ROC'", 
                        style={'textAlign': 'center', 'color': '#7f8c8d'}), "")
    
    try:
        roc_data = filter_and_calculate_roc(
            COLLECTED_DATA, disease, method, ff_filter, length_filter,
            custom_length if length_filter == 'Custom' else None
        )
        
        if roc_data is None:
            return (html.Div("No data available", style={'color': '#e74c3c', 'textAlign': 'center'}), "")
        
        # Create ROC plot
        fig = go.Figure()
        
        # Prepare customdata (threshold values for each point)
        hover_text = [f"Threshold: {t:.2f}" for t in roc_data['thresholds']]
        
        fig.add_trace(go.Scatter(
            x=roc_data['fpr'], 
            y=roc_data['tpr'],
            mode='lines+markers',
            name=f'ROC (AUC={roc_data["auc"]:.4f})',
            line=dict(color='#3498db', width=3),
            marker=dict(size=4),
            text=hover_text,
            customdata=roc_data['thresholds'],
            hovertemplate='<b>FPR</b>: %{x:.4f}<br>' +
                         '<b>TPR</b>: %{y:.4f}<br>' +
                         '<b>Z-score Threshold</b>: %{customdata:.2f}<br>' +
                         '<extra></extra>'
        ))
        
        fig.add_trace(go.Scatter(
            x=[0, 1], y=[0, 1], mode='lines',
            name='Random', line=dict(color='gray', width=2, dash='dash')
        ))
        
        filter_desc = []
        if ff_filter != 'All':
            filter_desc.append(f"FF={ff_filter}%")
        else:
            filter_desc.append("All FF")
        if length_filter == 'Custom':
            filter_desc.append(f"MinLen≥{custom_length}Mb")
        elif length_filter != 'All':
            filter_desc.append(f"MinLen≥{length_filter}Mb")
        else:
            filter_desc.append("All Lengths")
        
        fig.update_layout(
            title=f'ROC: {disease} (chr{DISEASE_CHROMOSOMES[disease]}) - {method}<br>{", ".join(filter_desc)}',
            xaxis_title='False Positive Rate (1-Specificity)',
            yaxis_title='True Positive Rate (Sensitivity)',
            height=550, 
            width=700,
            template='plotly_white',
            xaxis=dict(range=[0, 1]), 
            yaxis=dict(range=[0, 1.05]),
            showlegend=True, 
            legend=dict(x=0.6, y=0.1), 
            hovermode='closest',
            margin=dict(l=80, r=80, t=100, b=80)
        )
        
        # Stats (compact single line)
        stats_html = html.Div([
            html.Div([
                html.Span("AUC: ", style={'fontWeight': 'bold', 'color': '#7f8c8d'}),
                html.Span(f"{roc_data['auc']:.4f}", style={'color': '#27ae60', 'fontWeight': 'bold', 'fontSize': '16px'}),
                html.Span(" | ", style={'margin': '0 10px', 'color': '#bdc3c7'}),
                
                html.Span(f"Positive ({disease}): ", style={'fontWeight': 'bold', 'color': '#7f8c8d'}),
                html.Span(f"{roc_data['n_positive']}", style={'color': '#3498db', 'fontWeight': 'bold'}),
                html.Span(" | ", style={'margin': '0 10px', 'color': '#bdc3c7'}),
                
                html.Span("Negative (Others): ", style={'fontWeight': 'bold', 'color': '#7f8c8d'}),
                html.Span(f"{roc_data['n_negative']}", style={'color': '#e74c3c', 'fontWeight': 'bold'}),
                html.Span(" | ", style={'margin': '0 10px', 'color': '#bdc3c7'}),
                
                html.Span("Total: ", style={'fontWeight': 'bold', 'color': '#7f8c8d'}),
                html.Span(f"{roc_data['n_total']}", style={'color': '#9b59b6', 'fontWeight': 'bold'}),
            ], style={'backgroundColor': '#ecf0f1', 'borderRadius': '10px', 'padding': '15px', 'textAlign': 'center'})
        ])
        
        return (dcc.Graph(figure=fig, id='roc-graph', 
                         config={'displayModeBar': True, 'displaylogo': False}), 
                stats_html)
        
    except Exception as e:
        import traceback
        return (html.Div(f"Error: {str(e)}\n{traceback.format_exc()}", 
                        style={'color': '#e74c3c', 'whiteSpace': 'pre-wrap'}), "")


# Callback removed - click info functionality removed per user request


@app.callback(
    [Output('confusion-matrix-container', 'children'),
     Output('metrics-heatmaps-container', 'children')],
    Input('calc-perf-btn', 'n_clicks'),
    [State('disease-perf', 'value'),
     State('method-perf', 'value'),
     State('threshold-perf', 'value'),
     State('minlen-perf', 'value')]
)
def update_performance(n_clicks, disease, method, threshold, min_length_mb):
    if not n_clicks or n_clicks == 0:
        return (html.Div("Set parameters and click 'Calculate Performance'", 
                       style={'textAlign': 'center', 'color': '#7f8c8d'}), "")
    
    try:
        # Validate threshold input
        if threshold is None or threshold <= 0:
            return (html.Div("Please enter a valid z-score threshold (> 0)", 
                           style={'color': '#e74c3c', 'textAlign': 'center'}), "")
        
        # Fixed FF and Length values
        ff_values = [5.0, 10.0, 15.0]
        all_length_values = [0.5, 1, 3, 5, 7, 10]
        
        # Filter length values based on min_length_mb
        if min_length_mb is not None and min_length_mb > 0:
            length_values = [l for l in all_length_values if l >= min_length_mb]
        else:
            length_values = all_length_values
        
        all_metrics_data = []
        
        for ff in ff_values:
            for length in length_values:
                    
                cm_result = calculate_confusion_matrix(
                    COLLECTED_DATA, disease, method, ff, length, threshold
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
                        'n_pos': cm_result['n_pos'],
                        'n_neg': cm_result['n_neg']
                    })
        
        if len(all_metrics_data) == 0:
            return (html.Div("No data available", style={'color': '#e74c3c', 'textAlign': 'center'}), "")
        
        df_metrics = pd.DataFrame(all_metrics_data)
        
        # Calculate overall confusion matrix
        total_TP = df_metrics['TP'].sum()
        total_FP = df_metrics['FP'].sum()
        total_FN = df_metrics['FN'].sum()
        total_TN = df_metrics['TN'].sum()
        total_sensitivity = total_TP / (total_TP + total_FN) if (total_TP + total_FN) > 0 else 0
        total_specificity = total_TN / (total_TN + total_FP) if (total_TN + total_FP) > 0 else 0
        total_ppv = total_TP / (total_TP + total_FP) if (total_TP + total_FP) > 0 else 0
        total_npv = total_TN / (total_TN + total_FN) if (total_TN + total_FN) > 0 else 0
        
        # Create confusion matrix table
        cm_table = html.Div([
            html.H4(f"Overall Confusion Matrix: {disease} - {method} (Threshold={threshold})", 
                   style={'textAlign': 'center', 'marginBottom': '20px', 'color': '#2c3e50'}),
            html.Div([
                html.Table([
                    html.Thead(html.Tr([
                        html.Th("", style={'border': '2px solid #34495e', 'padding': '12px', 'backgroundColor': '#ecf0f1'}),
                        html.Th("Has the disease", style={'border': '2px solid #34495e', 'padding': '12px', 'backgroundColor': '#d5f4e6', 'fontWeight': 'bold'}),
                        html.Th("Does not have the disease", style={'border': '2px solid #34495e', 'padding': '12px', 'backgroundColor': '#fadbd8', 'fontWeight': 'bold'})
                    ])),
                    html.Tbody([
                        html.Tr([
                            html.Td("Test: Positive", style={'border': '2px solid #34495e', 'padding': '12px', 'backgroundColor': '#d5f4e6', 'fontWeight': 'bold'}),
                            html.Td(f"TP = {total_TP}", style={'border': '2px solid #34495e', 'padding': '12px', 'textAlign': 'center', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#27ae60'}),
                            html.Td(f"FP = {total_FP}", style={'border': '2px solid #34495e', 'padding': '12px', 'textAlign': 'center', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#e74c3c'})
                        ]),
                        html.Tr([
                            html.Td("Test: Negative", style={'border': '2px solid #34495e', 'padding': '12px', 'backgroundColor': '#fadbd8', 'fontWeight': 'bold'}),
                            html.Td(f"FN = {total_FN}", style={'border': '2px solid #34495e', 'padding': '12px', 'textAlign': 'center', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#e67e22'}),
                            html.Td(f"TN = {total_TN}", style={'border': '2px solid #34495e', 'padding': '12px', 'textAlign': 'center', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#3498db'})
                        ])
                    ])
                ], style={'margin': '0 auto', 'borderCollapse': 'collapse'}),
            ], style={'textAlign': 'center', 'marginBottom': '20px'}),
            
            html.Div([
                html.Span(f"Sensitivity = {total_sensitivity:.2%}", style={'margin': '0 15px', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#27ae60'}),
                html.Span(" | ", style={'margin': '0 10px', 'color': '#bdc3c7'}),
                html.Span(f"Specificity = {total_specificity:.2%}", style={'margin': '0 15px', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#3498db'}),
                html.Span(" | ", style={'margin': '0 10px', 'color': '#bdc3c7'}),
                html.Span(f"PPV = {total_ppv:.2%}", style={'margin': '0 15px', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#9b59b6'}),
                html.Span(" | ", style={'margin': '0 10px', 'color': '#bdc3c7'}),
                html.Span(f"NPV = {total_npv:.2%}", style={'margin': '0 15px', 'fontSize': '18px', 'fontWeight': 'bold', 'color': '#e67e22'})
            ], style={'textAlign': 'center', 'backgroundColor': '#ecf0f1', 'padding': '15px', 'borderRadius': '10px'})
        ], style={'marginBottom': '30px', 'padding': '20px', 'backgroundColor': 'white', 'borderRadius': '10px', 'boxShadow': '0 0 10px rgba(0,0,0,0.1)'})
        
        # Helper function to create metric heatmap
        def create_metric_heatmap(df, metric, colorscale, disease, method):
            pivot = df.pivot(index='Length', columns='FF', values=metric)
            # Get actual length values from the data
            actual_length_values = sorted(pivot.index.tolist())
            
            fig = go.Figure(data=go.Heatmap(
                z=pivot.values, x=pivot.columns, y=pivot.index,
                colorscale=colorscale, zmin=0, zmax=1,
                text=np.round(pivot.values * 100, 1),
                texttemplate='%{text}%', textfont={"size": 12},
                colorbar=dict(title=metric)
            ))
            fig.update_layout(
                title=f'{metric}: {disease} - {method}',
                xaxis_title='Fetal Fraction (%)', yaxis_title='Deletion Length (Mb)',
                xaxis={'tickmode': 'array', 'tickvals': [5, 10, 15]},
                yaxis={'tickmode': 'array', 'tickvals': actual_length_values},
                height=400, template='plotly_white'
            )
            return fig
        
        # Create 4 heatmaps
        heatmaps_html = html.Div([
            html.H4("Performance Metrics by FF and Deletion Length", 
                   style={'textAlign': 'center', 'marginTop': '20px', 'marginBottom': '20px', 'color': '#2c3e50'}),
            html.Div([
                html.Div([
                    dcc.Graph(figure=create_metric_heatmap(df_metrics, 'Sensitivity', 'RdYlGn', disease, method))
                ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                html.Div([
                    dcc.Graph(figure=create_metric_heatmap(df_metrics, 'Specificity', 'RdYlGn', disease, method))
                ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'})
            ]),
            html.Div([
                html.Div([
                    dcc.Graph(figure=create_metric_heatmap(df_metrics, 'PPV', 'Blues', disease, method))
                ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'}),
                html.Div([
                    dcc.Graph(figure=create_metric_heatmap(df_metrics, 'NPV', 'Greens', disease, method))
                ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px', 'verticalAlign': 'top'})
            ])
        ])
        
        return cm_table, heatmaps_html
        
    except Exception as e:
        import traceback
        return (html.Div(f"Error: {str(e)}\n{traceback.format_exc()}", 
                       style={'color': '#e74c3c', 'whiteSpace': 'pre-wrap'}), "")


def main():
    parser = argparse.ArgumentParser(description="ROC Dashboard V3")
    parser.add_argument('--data_file', type=str, 
                       default='/data/md_validation/roc_results/collected_data.csv')
    parser.add_argument('--port', type=int, default=8001)
    parser.add_argument('--host', type=str, default='0.0.0.0')
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    global COLLECTED_DATA
    data_path = Path(args.data_file)
    
    if not data_path.exists():
        logger.error(f"Data file not found: {data_path}")
        return 1
    
    COLLECTED_DATA = load_collected_data(data_path)
    
    logger.info("="*80)
    logger.info("Starting ROC Dashboard V3")
    logger.info("="*80)
    logger.info(f"URL: http://localhost:{args.port}")
    logger.info("="*80)
    
    app.run(host=args.host, port=args.port, debug=args.debug)
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

