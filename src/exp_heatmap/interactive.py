"""
Interactive visualization module for ExP Heatmap.

This module provides Plotly-based interactive heatmap visualizations
that enable zooming, panning, and hover tooltips for detailed exploration
of cross-population genomic data.
"""

import os
import pandas as pd
import numpy as np
from typing import Optional, List, Tuple, Union

# Check for plotly availability
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

from exp_heatmap.plot import (
    populations_1000genomes, 
    superpopulations,
    superpopulation_colors,
    pop_to_superpop,
    create_plot_input
)


def check_plotly():
    """Check if plotly is available and raise informative error if not."""
    if not PLOTLY_AVAILABLE:
        raise ImportError(
            "Plotly is required for interactive visualizations. "
            "Install it with: pip install plotly"
        )


def plot_interactive_heatmap(
    input_df: pd.DataFrame,
    start: int,
    end: int,
    title: Optional[str] = None,
    output: str = "ExP_heatmap_interactive",
    populations: Union[str, tuple] = "1000Genomes",
    colorscale: str = "Blues",
    zmin: Optional[float] = None,
    zmax: Optional[float] = None,
    display_limit: Optional[float] = None,
    display_values: str = "higher",
    height: int = 800,
    width: int = 1200,
    show_superpop_annotations: bool = True
) -> 'go.Figure':
    """
    Generate an interactive ExP heatmap using Plotly.
    
    This function creates an HTML-based interactive visualization that allows
    users to zoom, pan, and hover over data points for detailed information.
    
    Parameters
    ----------
    input_df : pandas.DataFrame
        Input data containing pairwise statistics or rank scores to visualize.
        Rows should be population pairs (e.g., "ACB_ASW"), columns should be
        genomic positions.
    start : int
        Start genomic position for the x-axis (region to display).
    end : int
        End genomic position for the x-axis (region to display).
    title : str, optional
        Title of the heatmap.
    output : str, optional
        Output filename without extension (default: "ExP_heatmap_interactive").
        The file will be saved as '{output}.html'.
    populations : str or tuple, optional
        Population set to use. Default "1000Genomes" uses standard 26 populations.
    colorscale : str, optional
        Plotly colorscale name (default: "Blues"). Options include:
        'Viridis', 'Plasma', 'Inferno', 'Magma', 'Blues', 'Reds', 'YlOrRd', etc.
    zmin, zmax : float, optional
        Minimum and maximum values for the color scale.
    display_limit : float, optional
        If set, values above or below this threshold are masked.
    display_values : str, optional
        'higher' or 'lower' - which values to retain when using display_limit.
    height : int, optional
        Height of the figure in pixels (default: 800).
    width : int, optional
        Width of the figure in pixels (default: 1200).
    show_superpop_annotations : bool, optional
        If True, add superpopulation group annotations on y-axis.
        
    Returns
    -------
    plotly.graph_objects.Figure
        Interactive Plotly figure object. Also saves to HTML file.
    """
    check_plotly()
    
    input_df = input_df.copy()
    
    # Filter to specified region
    columns = [c for c in input_df.columns if start <= c <= end]
    if not columns:
        raise ValueError(f"No data found in range {start}-{end}")
    input_df = input_df[columns]
    
    # Determine populations
    is_1000genomes = populations == "1000Genomes"
    if is_1000genomes:
        populations = populations_1000genomes
    
    # Apply display limit
    if display_limit is not None:
        if display_values == "higher":
            input_df = input_df.where(input_df >= display_limit, 0)
        elif display_values == "lower":
            input_df = input_df.where(input_df <= display_limit, 0)
    
    # Set color scale limits
    if zmin is None:
        zmin = 1.301 if is_1000genomes else input_df.min().min()
    if zmax is None:
        zmax = 4.833 if is_1000genomes else input_df.max().max()
    
    # Create custom hover text
    hover_text = []
    for row_idx, row_name in enumerate(input_df.index):
        row_hover = []
        pop1, pop2 = row_name.split("_")
        superpop1 = pop_to_superpop.get(pop1, "Unknown")
        superpop2 = pop_to_superpop.get(pop2, "Unknown")
        for col_idx, pos in enumerate(input_df.columns):
            value = input_df.iloc[row_idx, col_idx]
            row_hover.append(
                f"Position: {pos:,}<br>"
                f"Population pair: {row_name}<br>"
                f"Superpopulations: {superpop1} vs {superpop2}<br>"
                f"Score: {value:.3f}"
            )
        hover_text.append(row_hover)
    
    # Create the heatmap
    fig = go.Figure()
    
    fig.add_trace(go.Heatmap(
        z=input_df.values,
        x=input_df.columns,
        y=input_df.index,
        colorscale=colorscale,
        zmin=zmin,
        zmax=zmax,
        hovertext=hover_text,
        hoverinfo='text',
        colorbar=dict(
            title=dict(text="-log10(rank score)", side="right"),
            thickness=15,
            len=0.9
        )
    ))
    
    # Add superpopulation annotations if enabled and using 1000 Genomes
    if show_superpop_annotations and is_1000genomes:
        n_pops = len(populations)
        n_pairs_per_pop = n_pops - 1
        
        # Add horizontal lines between superpopulation groups
        superpop_boundaries = []
        current_pos = 0
        for superpop in ["AFR", "SAS", "EAS", "EUR", "AMR"]:
            pops_in_superpop = len(superpopulations[superpop])
            superpop_size = pops_in_superpop * n_pairs_per_pop
            superpop_boundaries.append({
                'name': superpop,
                'start': current_pos,
                'end': current_pos + superpop_size,
                'mid': current_pos + superpop_size / 2
            })
            current_pos += superpop_size
        
        # Add annotations for superpopulation groups
        annotations = []
        for boundary in superpop_boundaries:
            annotations.append(dict(
                x=-0.02,
                y=boundary['mid'] / len(input_df),
                xref='paper',
                yref='paper',
                text=boundary['name'],
                showarrow=False,
                font=dict(size=10, color=superpopulation_colors[boundary['name']]),
                textangle=-90
            ))
        
        fig.update_layout(annotations=annotations)
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=title or "ExP Heatmap (Interactive)",
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title=f"Genomic Position ({start:,} - {end:,})",
            tickformat=",d",
            showgrid=False
        ),
        yaxis=dict(
            title="Population Pairs",
            showgrid=False,
            tickfont=dict(size=6) if len(input_df) > 100 else dict(size=8)
        ),
        height=height,
        width=width,
        margin=dict(l=100, r=50, t=80, b=80)
    )
    
    # Add range slider for x-axis
    fig.update_xaxes(
        rangeslider=dict(visible=True, thickness=0.05),
        type="linear"
    )
    
    # Save to HTML
    output_file = f"{output}.html"
    fig.write_html(output_file, include_plotlyjs=True, full_html=True)
    print(f"Interactive heatmap saved to: {output_file}")
    
    return fig


def plot_interactive(
    input_dir: str,
    start: int,
    end: int,
    title: Optional[str] = None,
    output: str = "ExP_heatmap_interactive",
    **kwargs
) -> 'go.Figure':
    """
    High-level function to create interactive heatmap from computed results.
    
    This is the main entry point for creating interactive visualizations
    from a directory of TSV files produced by the compute step.
    
    Parameters
    ----------
    input_dir : str
        Path to directory containing TSV files from 'exp_heatmap compute'.
    start : int
        Start genomic position.
    end : int
        End genomic position.
    title : str, optional
        Title for the heatmap.
    output : str, optional
        Output filename without extension.
    **kwargs : dict
        Additional arguments passed to plot_interactive_heatmap().
        
    Returns
    -------
    plotly.graph_objects.Figure
        Interactive Plotly figure object.
    """
    check_plotly()
    
    # Load and prepare data
    plot_input = create_plot_input(input_dir, start=start, end=end)
    
    return plot_interactive_heatmap(
        plot_input,
        start=plot_input.columns[0],
        end=plot_input.columns[-1],
        title=title,
        output=output,
        **kwargs
    )


def create_comparison_view(
    input_df: pd.DataFrame,
    region1: Tuple[int, int],
    region2: Tuple[int, int],
    title: Optional[str] = None,
    output: str = "ExP_comparison",
    **kwargs
) -> 'go.Figure':
    """
    Create a side-by-side comparison of two genomic regions.
    
    Useful for comparing selection signals between different genomic regions
    or genes.
    
    Parameters
    ----------
    input_df : pandas.DataFrame
        Input data with population pairs as rows and positions as columns.
    region1 : tuple
        (start, end) positions for first region.
    region2 : tuple
        (start, end) positions for second region.
    title : str, optional
        Title for the comparison figure.
    output : str, optional
        Output filename without extension.
    **kwargs : dict
        Additional arguments for heatmap styling.
        
    Returns
    -------
    plotly.graph_objects.Figure
        Interactive comparison figure.
    """
    check_plotly()
    
    # Extract regions
    cols1 = [c for c in input_df.columns if region1[0] <= c <= region1[1]]
    cols2 = [c for c in input_df.columns if region2[0] <= c <= region2[1]]
    
    df1 = input_df[cols1]
    df2 = input_df[cols2]
    
    # Create subplots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[
            f"Region 1: {region1[0]:,}-{region1[1]:,}",
            f"Region 2: {region2[0]:,}-{region2[1]:,}"
        ],
        horizontal_spacing=0.1
    )
    
    colorscale = kwargs.get('colorscale', 'Blues')
    zmin = kwargs.get('zmin', min(df1.min().min(), df2.min().min()))
    zmax = kwargs.get('zmax', max(df1.max().max(), df2.max().max()))
    
    # Add first heatmap
    fig.add_trace(
        go.Heatmap(
            z=df1.values,
            x=df1.columns,
            y=df1.index,
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            showscale=False
        ),
        row=1, col=1
    )
    
    # Add second heatmap with colorbar
    fig.add_trace(
        go.Heatmap(
            z=df2.values,
            x=df2.columns,
            y=df2.index,
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            colorbar=dict(title="-log10(rank score)")
        ),
        row=1, col=2
    )
    
    fig.update_layout(
        title=dict(text=title or "Region Comparison", x=0.5),
        height=kwargs.get('height', 800),
        width=kwargs.get('width', 1600)
    )
    
    # Save
    output_file = f"{output}.html"
    fig.write_html(output_file)
    print(f"Comparison view saved to: {output_file}")
    
    return fig


def create_population_focus_view(
    input_df: pd.DataFrame,
    focus_population: str,
    start: int,
    end: int,
    title: Optional[str] = None,
    output: str = "ExP_population_focus",
    **kwargs
) -> 'go.Figure':
    """
    Create a focused view showing all comparisons involving a specific population.
    
    Parameters
    ----------
    input_df : pandas.DataFrame
        Input data with population pairs as rows.
    focus_population : str
        Population code to focus on (e.g., "CEU", "YRI").
    start : int
        Start genomic position.
    end : int
        End genomic position.
    title : str, optional
        Title for the figure.
    output : str, optional
        Output filename without extension.
    **kwargs : dict
        Additional styling arguments.
        
    Returns
    -------
    plotly.graph_objects.Figure
        Interactive figure focused on one population.
    """
    check_plotly()
    
    # Filter to rows containing the focus population
    mask = input_df.index.str.contains(focus_population)
    filtered_df = input_df[mask]
    
    if filtered_df.empty:
        raise ValueError(f"No data found for population: {focus_population}")
    
    # Filter columns to region
    cols = [c for c in filtered_df.columns if start <= c <= end]
    filtered_df = filtered_df[cols]
    
    # Create heatmap
    fig = go.Figure()
    
    fig.add_trace(go.Heatmap(
        z=filtered_df.values,
        x=filtered_df.columns,
        y=filtered_df.index,
        colorscale=kwargs.get('colorscale', 'Blues'),
        colorbar=dict(title="-log10(rank score)")
    ))
    
    superpop = pop_to_superpop.get(focus_population, "Unknown")
    
    fig.update_layout(
        title=dict(
            text=title or f"Selection signals involving {focus_population} ({superpop})",
            x=0.5
        ),
        xaxis=dict(title=f"Genomic Position ({start:,} - {end:,})", tickformat=",d"),
        yaxis=dict(title="Population Pairs"),
        height=kwargs.get('height', 600),
        width=kwargs.get('width', 1200)
    )
    
    # Add range slider
    fig.update_xaxes(rangeslider=dict(visible=True, thickness=0.05))
    
    output_file = f"{output}.html"
    fig.write_html(output_file)
    print(f"Population focus view saved to: {output_file}")
    
    return fig


