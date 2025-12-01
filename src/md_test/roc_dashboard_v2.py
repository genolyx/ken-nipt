#!/usr/bin/env python3
"""
ROC Analysis Dashboard V2 - Using Pre-calculated Results

Interactive dashboard for exploring ROC analysis results with flexible filtering.

Features:
- 7 detection modes (wc_orig, wc_fetus, wcx_orig, wcx_fetus, orig, fetus, any)
- 8 diseases (One-vs-Rest)
- FF filtering (5%, 10%, 15%, All)
- Min Length filtering (1, 3, 5, 7, 10 Mb, All, Custom)
- Real-time ROC curve recalculation

Usage:
    python roc_dashboard_v2.py --port 8051
    
Then open browser: http://localhost:8051
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
import plotly.express as px
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
    '1p36': '1',
    '2q33': '2',
    'CDC': '5',
    'DGS': '22',
    'Jacobsen': '11',
    'PWS': '15',
    'WBS': '7',
    'WHS': '4'
}

# Global data storage
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
    
    tpr_list = []
    fpr_list = []
    
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
    """Calculate AUC using trapezoidal rule"""
    sorted_indices = np.argsort(fpr)
    fpr_sorted = fpr[sorted_indices]
    tpr_sorted = tpr[sorted_indices]
    return np.trapz(tpr_sorted, fpr_sorted)


def filter_and_calculate_roc(
    df: pd.DataFrame,
    target_disease: str,
    mode: str,
    ff_filter: str = 'All',
    length_filter: str = 'All',
    custom_length: float = None
) -> Dict:
    """Filter data and calculate ROC curve
    
    Args:
        df: Wide format dataframe with all samples
        target_disease: Target disease
        mode: Detection mode (individual or group)
        ff_filter: FF filter ('5', '10', '15', 'All')
        length_filter: Length filter ('1', '3', '5', '7', '10', 'All', 'Custom')
        custom_length: Custom minimum length in Mb
    
    Returns:
        Dictionary with fpr, tpr, thresholds, auc, n_samples
    """
    # Filter by target disease
    df_target = df[df['target_disease'] == target_disease].copy()
    
    if len(df_target) == 0:
        return None
    
    # Filter by FF
    if ff_filter != 'All':
        ff_value = float(ff_filter)
        df_target = df_target[df_target['ff'] == ff_value]
    
    # Filter by length
    if length_filter == 'Custom' and custom_length is not None:
        df_target = df_target[df_target['deletion_length_mb'] >= custom_length]
    elif length_filter != 'All':
        length_value = int(length_filter)
        df_target = df_target[df_target['deletion_length_mb'] >= length_value]
    
    if len(df_target) == 0:
        logger.warning(f"No samples after filtering")
        return None
    
    # Create labels
    y_true = (df_target['disease'] == target_disease).astype(int).values
    
    # Get scores based on mode
    if mode in INDIVIDUAL_MODES:
        if mode not in df_target.columns:
            return None
        y_score = df_target[mode].values
        
    elif mode == 'orig':
        cols = [c for c in ['wc_orig', 'wcx_orig'] if c in df_target.columns]
        if not cols:
            return None
        y_score = df_target[cols].max(axis=1).values
        
    elif mode == 'fetus':
        cols = [c for c in ['wc_fetus', 'wcx_fetus'] if c in df_target.columns]
        if not cols:
            return None
        y_score = df_target[cols].max(axis=1).values
        
    elif mode == 'any':
        cols = [c for c in INDIVIDUAL_MODES if c in df_target.columns]
        if not cols:
            return None
        y_score = df_target[cols].max(axis=1).values
    else:
        return None
    
    # Check if we have both classes
    if len(np.unique(y_true)) < 2:
        return None
    
    # Calculate ROC
    fpr, tpr, thresholds = calculate_roc_curve_manual(y_true, y_score)
    roc_auc = calculate_auc(fpr, tpr)
    
    return {
        'fpr': fpr,
        'tpr': tpr,
        'thresholds': thresholds,
        'auc': roc_auc,
        'n_positive': int(y_true.sum()),
        'n_negative': int((1 - y_true).sum()),
        'n_total': len(df_target)
    }


# Initialize Dash app
app = dash.Dash(__name__, suppress_callback_exceptions=True)

app.layout = html.Div([
    html.H1("ROC Analysis Dashboard V2", style={'textAlign': 'center', 'color': '#2c3e50'}),
    
    html.Div([
        html.Div([
            html.Label("Disease:", style={'fontWeight': 'bold'}),
            dcc.Dropdown(
                id='disease-dropdown',
                options=[{'label': f'{d} (chr{DISEASE_CHROMOSOMES[d]})', 'value': d} for d in DISEASES],
                value='1p36',
                clearable=False
            )
        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px'}),
        
        html.Div([
            html.Label("Method:", style={'fontWeight': 'bold'}),
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
                value='wc_orig',
                clearable=False
            )
        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px'}),
        
        html.Div([
            html.Label("Fetal Fraction:", style={'fontWeight': 'bold'}),
            dcc.Dropdown(
                id='ff-dropdown',
                options=[
                    {'label': 'All FF', 'value': 'All'},
                    {'label': 'FF = 5%', 'value': '5'},
                    {'label': 'FF = 10%', 'value': '10'},
                    {'label': 'FF = 15%', 'value': '15'}
                ],
                value='All',
                clearable=False
            )
        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px'}),
        
        html.Div([
            html.Label("Min Deletion Length:", style={'fontWeight': 'bold'}),
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
                value='All',
                clearable=False
            )
        ], style={'width': '23%', 'display': 'inline-block', 'padding': '10px'}),
    ], style={'backgroundColor': '#ecf0f1', 'padding': '20px', 'borderRadius': '10px'}),
    
    # Custom length input (hidden by default)
    html.Div([
        html.Label("Custom Min Length (Mb):", style={'fontWeight': 'bold', 'marginRight': '10px'}),
        dcc.Input(
            id='custom-length-input',
            type='number',
            value=2,
            min=0.1,
            max=15,
            step=0.1,
            style={'width': '100px'}
        )
    ], id='custom-length-div', style={'display': 'none', 'padding': '10px', 'textAlign': 'center'}),
    
    html.Br(),
    
    html.Div([
        html.Button(
            'Calculate ROC',
            id='calculate-button',
            n_clicks=0,
            style={
                'fontSize': '16px',
                'padding': '10px 30px',
                'backgroundColor': '#3498db',
                'color': 'white',
                'border': 'none',
                'borderRadius': '5px',
                'cursor': 'pointer'
            }
        )
    ], style={'textAlign': 'center', 'padding': '20px'}),
    
    dcc.Loading(
        id="loading-roc",
        type="default",
        children=[
            html.Div(id='roc-output', style={'padding': '20px'})
        ]
    ),
    
    dcc.Loading(
        id="loading-stats",
        type="circle",
        children=[
            html.Div(id='stats-output', style={'padding': '20px'})
        ]
    )
])


@app.callback(
    Output('custom-length-div', 'style'),
    Input('length-dropdown', 'value')
)
def toggle_custom_length(length_value):
    """Show/hide custom length input"""
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
    """Update ROC curve based on filters"""
    if n_clicks == 0:
        return html.Div("Select parameters and click 'Calculate ROC'", 
                       style={'textAlign': 'center', 'color': '#7f8c8d'}), ""
    
    try:
        # Calculate ROC
        roc_data = filter_and_calculate_roc(
            COLLECTED_DATA,
            disease,
            method,
            ff_filter,
            length_filter,
            custom_length if length_filter == 'Custom' else None
        )
        
        if roc_data is None:
            return html.Div("No data available for selected parameters", 
                           style={'color': '#e74c3c', 'textAlign': 'center'}), ""
        
        # Create ROC plot
        fig = go.Figure()
        
        # ROC curve
        fig.add_trace(go.Scatter(
            x=roc_data['fpr'],
            y=roc_data['tpr'],
            mode='lines',
            name=f'ROC (AUC={roc_data["auc"]:.4f})',
            line=dict(color='#3498db', width=3)
        ))
        
        # Diagonal line
        fig.add_trace(go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode='lines',
            name='Random',
            line=dict(color='gray', width=2, dash='dash')
        ))
        
        # Filter description
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
        
        filter_text = ", ".join(filter_desc)
        
        fig.update_layout(
            title=f'ROC Curve: {disease} (chr{DISEASE_CHROMOSOMES[disease]}) - {method}<br>{filter_text}',
            xaxis_title='False Positive Rate (1-Specificity)',
            yaxis_title='True Positive Rate (Sensitivity)',
            height=600,
            template='plotly_white',
            xaxis=dict(range=[0, 1]),
            yaxis=dict(range=[0, 1.05]),
            showlegend=True,
            legend=dict(x=0.6, y=0.1)
        )
        
        # Statistics
        stats_html = html.Div([
            html.H3("Performance Metrics", style={'color': '#2c3e50'}),
            html.Div([
                html.Div([
                    html.H4(f"{roc_data['auc']:.4f}", style={'color': '#27ae60', 'margin': '0'}),
                    html.P("AUC", style={'margin': '0', 'color': '#7f8c8d'})
                ], style={'display': 'inline-block', 'padding': '20px', 'textAlign': 'center', 'width': '20%'}),
                
                html.Div([
                    html.H4(f"{roc_data['n_positive']}", style={'color': '#3498db', 'margin': '0'}),
                    html.P(f"Positive ({disease})", style={'margin': '0', 'color': '#7f8c8d'})
                ], style={'display': 'inline-block', 'padding': '20px', 'textAlign': 'center', 'width': '25%'}),
                
                html.Div([
                    html.H4(f"{roc_data['n_negative']}", style={'color': '#e74c3c', 'margin': '0'}),
                    html.P("Negative (Others)", style={'margin': '0', 'color': '#7f8c8d'})
                ], style={'display': 'inline-block', 'padding': '20px', 'textAlign': 'center', 'width': '25%'}),
                
                html.Div([
                    html.H4(f"{roc_data['n_total']}", style={'color': '#9b59b6', 'margin': '0'}),
                    html.P("Total Samples", style={'margin': '0', 'color': '#7f8c8d'})
                ], style={'display': 'inline-block', 'padding': '20px', 'textAlign': 'center', 'width': '20%'}),
            ], style={'backgroundColor': '#ecf0f1', 'borderRadius': '10px', 'padding': '20px'})
        ])
        
        return dcc.Graph(figure=fig), stats_html
        
    except Exception as e:
        import traceback
        error_msg = f"Error: {str(e)}\n{traceback.format_exc()}"
        logger.error(error_msg)
        return html.Div(error_msg, style={'color': '#e74c3c', 'whiteSpace': 'pre-wrap'}), ""


def main():
    parser = argparse.ArgumentParser(
        description="ROC Analysis Dashboard V2"
    )
    parser.add_argument(
        '--data_file',
        type=str,
        default='/data/md_validation/roc_results/collected_data.csv',
        help='Path to collected_data.csv'
    )
    parser.add_argument(
        '--port',
        type=int,
        default=8051,
        help='Port to run dashboard (default: 8051)'
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
    
    # Load data
    global COLLECTED_DATA
    data_path = Path(args.data_file)
    
    if not data_path.exists():
        logger.error(f"Data file not found: {data_path}")
        logger.error("Please run roc_analysis.py first to generate the data")
        return 1
    
    COLLECTED_DATA = load_collected_data(data_path)
    
    logger.info("="*80)
    logger.info("Starting ROC Analysis Dashboard V2")
    logger.info("="*80)
    logger.info(f"URL: http://localhost:{args.port}")
    logger.info("Press Ctrl+C to stop")
    logger.info("="*80)
    
    app.run(
        host=args.host,
        port=args.port,
        debug=args.debug
    )
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

