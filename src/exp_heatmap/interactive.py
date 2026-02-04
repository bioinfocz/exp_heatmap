"""
Interactive visualization module for ExP Heatmap.

This module provides Plotly-based interactive heatmap visualizations
that enable zooming, panning, and hover tooltips for detailed exploration
of cross-population genomic data.
"""

import pandas as pd
from typing import Optional, Tuple, Union

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from exp_heatmap.plot import (
    populations_1000genomes, 
    superpopulations,
    pop_to_superpop,
    create_plot_input
)
from exp_heatmap.logging import get_logger

logger = get_logger(__name__)


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
    height: int = 600,
    width: int = 1500,
    show_superpop_annotations: bool = True,
    max_columns: int = 30000
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
        Height of the figure in pixels (default: 600).
    width : int, optional
        Width of the figure in pixels (default: 1500).
    show_superpop_annotations : bool, optional
        If True, add superpopulation group annotations on y-axis.
    max_columns : int, optional
        Maximum number of genomic positions (columns) to display (default: 30000).
        If the data has more columns, it will be downsampled by taking the maximum
        value within each bin. This keeps the file size manageable for browsers.
        Set to None to disable downsampling (not recommended for large regions).
        
    Returns
    -------
    plotly.graph_objects.Figure
        Interactive Plotly figure object. Also saves to HTML file.
    """
    input_df = input_df.copy()
    
    # Filter to specified region
    columns = [c for c in input_df.columns if start <= c <= end]
    if not columns:
        raise ValueError(f"No data found in range {start}-{end}")
    input_df = input_df[columns]
    
    # Downsample if too many variants
    n_cols = len(input_df.columns)
    if max_columns is not None and n_cols > max_columns:
        bin_size = n_cols // max_columns
        downsampled_data = []
        downsampled_positions = []
        
        for i in range(max_columns):
            bin_start = i * bin_size
            bin_end = min((i + 1) * bin_size, n_cols)
            bin_data = input_df.iloc[:, bin_start:bin_end]
            # Take max across the bin for each row (preserves selection signals)
            downsampled_data.append(bin_data.max(axis=1))
            # Use the middle position of the bin as the representative position
            mid_idx = (bin_start + bin_end) // 2
            downsampled_positions.append(input_df.columns[mid_idx])
        
        input_df = pd.DataFrame(downsampled_data, index=downsampled_positions).T
        input_df.columns = downsampled_positions
        
        logger.debug(f"Downsampled from {n_cols:,} to {len(input_df.columns):,} positions for interactive view")
    
    is_1000genomes = populations == "1000Genomes"
    if is_1000genomes:
        populations = populations_1000genomes
    
    if display_limit is not None:
        if display_values == "higher":
            input_df = input_df.where(input_df >= display_limit, 0)
        elif display_values == "lower":
            input_df = input_df.where(input_df <= display_limit, 0)
    
    if zmin is None:
        zmin = 1.301 if is_1000genomes else input_df.min().min()
    if zmax is None:
        zmax = 4.833 if is_1000genomes else input_df.max().max()
    
    # Create the heatmap
    fig = go.Figure()
    
    fig.add_trace(go.Heatmap(
        z=input_df.values,
        x=input_df.columns,
        y=input_df.index,
        colorscale=colorscale,
        zmin=zmin,
        zmax=zmax,
        hovertemplate='Position: %{x:,}<br>Population pair: %{y}<br>Score: %{z:.3f}<extra></extra>',
        colorbar=dict(
            title=dict(text="-log10(rank score)", side="right"),
            thickness=15,
            len=0.9
        )
    ))
    
    # Add superpopulation annotations
    if show_superpop_annotations and is_1000genomes:
        n_pops = len(populations)
        n_pairs_per_pop = n_pops - 1
        total_rows = len(input_df)
        
        superpop_order = ["AFR", "SAS", "EAS", "EUR", "AMR"]
        superpop_boundaries = []
        current_row = 0
        for superpop in superpop_order:
            pops_in_superpop = len(superpopulations[superpop])
            superpop_rows = pops_in_superpop * n_pairs_per_pop
            superpop_boundaries.append({
                'name': superpop,
                'start_row': current_row,
                'end_row': current_row + superpop_rows,
                'mid_row': current_row + superpop_rows / 2
            })
            current_row += superpop_rows
        
        annotations = []
        for boundary in superpop_boundaries:
            y_paper = 1.0 - (boundary['mid_row'] / total_rows)
            annotations.append(dict(
                x=1.02,
                y=y_paper,
                xref='paper',
                yref='paper',
                text=boundary['name'],
                showarrow=False,
                font=dict(size=10, color='black'),
                textangle=90
            ))
        
        shapes = []
        for i, boundary in enumerate(superpop_boundaries[:-1]):
            boundary_row_idx = boundary['end_row']
            if boundary_row_idx < total_rows:
                boundary_row_name = input_df.index[boundary_row_idx]
                shapes.append(dict(
                    type='line',
                    x0=0,
                    x1=1,
                    y0=boundary_row_name,
                    y1=boundary_row_name,
                    xref='paper',
                    yref='y',
                    line=dict(color='#CCCCCC', width=1)
                ))
        
        fig.update_layout(annotations=annotations, shapes=shapes)
    
    # Update layout for 1k genomes
    if is_1000genomes:
        n_pops = len(populations)
        n_pairs_per_pop = n_pops - 1
        
        tickvals = []
        ticktext = []
        for i, pop in enumerate(populations):
            center_row_idx = int(i * n_pairs_per_pop + (n_pairs_per_pop) / 2)
            if center_row_idx < len(input_df):
                tickvals.append(input_df.index[center_row_idx])
                ticktext.append(pop)
        
        yaxis_config = dict(
            title="",
            showgrid=False,
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext,
            tickfont=dict(size=6),
            side='left',
            autorange='reversed'
        )
    else:
        yaxis_config = dict(
            title="",
            showgrid=False,
            tickfont=dict(size=6) if len(input_df) > 100 else dict(size=8),
            autorange='reversed'
        )
    
    fig.update_layout(
        title=dict(
            text=title or "ExP Heatmap (Interactive)",
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title=f"{start:,} - {end:,}",
            tickformat=",d",
            showgrid=False
        ),
        yaxis=yaxis_config,
        height=height,
        width=width,
        margin=dict(l=100, r=50, t=80, b=80),
        dragmode='zoom'
    )
    
    # Add range slider for x-axis
    fig.update_xaxes(
        rangeslider=dict(visible=True, thickness=0.05),
        type="linear"
    )
    
    # Save to HTML
    output_file = f"{output}.html"
    centering_css = """
    <style>
        body {
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
            margin: 0;
            padding: 20px;
            box-sizing: border-box;
            background-color: #f5f5f5;
        }
        .plotly-graph-div {
            margin: auto;
        }
    </style>
    """
    
    html_content = fig.to_html(include_plotlyjs=True, full_html=True)
    html_content = html_content.replace('</head>', f'{centering_css}</head>')
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    logger.info(f"Interactive heatmap saved to: {output_file}")
    
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
    logger.info(f"Comparison view saved to: {output_file}")
    
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
    logger.info(f"Population focus view saved to: {output_file}")
    
    return fig

