import os
import sys
import matplotlib
import seaborn as sns
import pandas as pd
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
from bisect import bisect_left
from tqdm import tqdm

from exp_heatmap.logging import get_logger

logger = get_logger(__name__)

# 1000 Genomes Project data
populations_1000genomes = ("ACB","ASW","ESN","GWD","LWK","MSL","YRI","BEB","GIH","ITU","PJL","STU","CDX","CHB","CHS","JPT","KHV","CEU","FIN","GBR","IBS","TSI","CLM","MXL","PEL","PUR")
superpopulations = {"AFR": ("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"),
                    "SAS": ("BEB", "GIH", "ITU", "PJL", "STU"),
                    "EAS": ("CDX", "CHB", "CHS", "JPT", "KHV"),
                    "EUR": ("CEU", "FIN", "GBR", "IBS", "TSI"),
                    "AMR": ("CLM", "MXL", "PEL", "PUR",)}


# Mapping from population to superpopulation
pop_to_superpop = {}
for superpop, pops in superpopulations.items():
    for pop in pops:
        pop_to_superpop[pop] = superpop

RANK_SCORE_MODE_ALIASES = {
    "ascending": "ascending",
    "descending": "descending",
    "2-tailed": "directional",
    "directional": "directional",
}


def _extract_pop_pairs_from_files(segment_files):
    pop_pairs = []
    for segment_file in sorted(segment_files):
        pop_pair = os.path.splitext(os.path.basename(segment_file))[0].split(".")[0]
        parts = pop_pair.split("_")
        if len(parts) != 2:
            logger.debug(f"Skipping unrecognized TSV filename format: {segment_file}")
            continue
        pop_pairs.append(tuple(parts))
    return pop_pairs


def infer_populations_from_pairs(pop_pairs):
    seen = []
    for p1, p2 in pop_pairs:
        for pop in (p1, p2):
            if pop not in seen:
                seen.append(pop)
    return tuple(seen)


def resolve_population_configuration(populations, segment_files):
    if populations != "1000Genomes":
        return tuple(populations)

    pop_pairs = _extract_pop_pairs_from_files(segment_files)
    inferred = infer_populations_from_pairs(pop_pairs)
    if (
        len(inferred) == len(populations_1000genomes)
        and set(inferred) == set(populations_1000genomes)
    ):
        return "1000Genomes"

    return inferred


def normalize_rank_score_mode(rank_scores):
    try:
        return RANK_SCORE_MODE_ALIASES[rank_scores]
    except KeyError as exc:
        valid_modes = "', '".join(RANK_SCORE_MODE_ALIASES)
        raise ValueError(
            f"Unknown value for 'rank_scores'. Expected one of '{valid_modes}', got '{rank_scores}'"
        ) from exc


def auto_max_columns(figsize, dpi):
    """
    Estimate a sensible static heatmap column budget from the export width.
    """
    return max(1000, int(figsize[0] * dpi))


def downsample_heatmap_columns(input_df, max_columns, aggregation="max"):
    """
    Downsample wide genomic windows into a fixed number of representative columns.

    The representative x-position for each bin is the midpoint genomic coordinate of
    the source columns. Aggregation defaults to 'max' to preserve sharp local signals.
    """
    if max_columns is None or max_columns <= 0 or input_df.shape[1] <= max_columns:
        return input_df, False

    aggregation = aggregation.lower()
    reducers = {
        "max": lambda frame: frame.max(axis=1, skipna=True),
        "mean": lambda frame: frame.mean(axis=1, skipna=True),
        "median": lambda frame: frame.median(axis=1, skipna=True),
    }
    if aggregation not in reducers:
        raise ValueError(
            f"Unknown aggregation '{aggregation}'. Expected one of: {', '.join(reducers)}"
        )

    reducer = reducers[aggregation]
    bins = np.array_split(np.arange(input_df.shape[1]), max_columns)
    downsampled_rows = []
    representative_positions = []

    for column_indices in bins:
        if len(column_indices) == 0:
            continue
        bin_frame = input_df.iloc[:, column_indices]
        bin_positions = input_df.columns[column_indices]
        representative_positions.append(int(bin_positions[len(bin_positions) // 2]))
        downsampled_rows.append(reducer(bin_frame))

    downsampled_df = pd.DataFrame(downsampled_rows, index=representative_positions).T
    downsampled_df.columns = representative_positions
    return downsampled_df, True

def create_plot_input(input_dir, start, end, populations="1000Genomes", rank_scores="directional"):
    """
    Generate a pandas DataFrame for plotting from a directory of pairwise population .tsv files.

    Parameters:
        input_dir (str): Directory containing .tsv files named as 'POP1_POP2.*.tsv' for each population pair.
        start (int): Start genomic position (X-axis lower bound).
        end (int): End genomic position (X-axis upper bound).
        populations (iterable or "1000Genomes"): List of population names to include and order in the heatmap. If not using 1000 Genomes, provide a custom iterable.
        rank_scores (str): Determines which rank-score mode to use:
            - "ascending": Use ascending rank scores (lowest test values ranked first).
            - "descending": Use descending rank scores (highest test values ranked first).
            - "directional" / legacy "2-tailed": For pop1_pop2, use descending; for pop2_pop1, use ascending.
            
            Note: TSV column names contain "p_value" for backward compatibility, but these are
            empirical rank scores (genome-wide percentile ranks transformed to -log10 scale),
            not classical p-values.

    Returns:
        pd.DataFrame: Formatted data for heatmap plotting.
    """
    df_list = []
    pop_id_list = []
    segment_files = sorted(glob.glob(os.path.join(input_dir, "*.tsv")))
    populations = resolve_population_configuration(populations, segment_files)
    population_sorter = populations_1000genomes if populations == "1000Genomes" else populations
    index = 1
    
    # Load population pair files between start and end and check if the populations are in the population_sorter
    _first_file_range = None
    for segment_file in tqdm(segment_files, desc="Loading population pair files"):
        pop_pair = os.path.splitext(os.path.basename(segment_file))[0].split(".")[0]
        p1, p2 = pop_pair.split("_")
        
        if all(x in population_sorter for x in (p1, p2)):
            pop_id_list.append(pop_pair)

            segments = pd.read_csv(segment_file, sep="\t")
            if _first_file_range is None:
                _first_file_range = (int(segments.variant_pos.min()), int(segments.variant_pos.max()))
            segments = segments[(segments.variant_pos >= start) & (segments.variant_pos <= end)]

            df_list.append(segments)

        else:
            logger.debug(f"[{index}/{len(segment_files)}] Skipping {pop_pair} from {segment_file}. {p1} or {p2} not in provided population list.")
            
        index += 1

    if not df_list:
        raise ValueError(f"No compatible TSV files were found in '{input_dir}'.")
            
    # Validate the dimensions and variant_pos in each dataframe
    df_shape = df_list[0].shape
    variant_positions = df_list[0].variant_pos.values

    for i in range(len(df_list)):
        if df_list[i].shape != df_shape:
            raise ValueError(f"the shapes dont match in df {i}")
        if not np.array_equal(df_list[i].variant_pos.values, variant_positions):
            raise ValueError(f"the variant_positions dont match in df {i}")

    # Select variant_pos and rank score column, then transpose each df
    transp_list = []
    rank_mode = normalize_rank_score_mode(rank_scores)

    for df, pop_pair in zip(df_list, pop_id_list):
        # Determine which columns to use based on rank_scores parameter
        # Note: TSV column names contain "p_value" for backward compatibility
        if rank_mode == "ascending":
            score_column1 = "-log10_p_value_ascending"
            score_column2 = "-log10_p_value_ascending"
        
        elif rank_mode == "descending":
            score_column1 = "-log10_p_value_descending"
            score_column2 = "-log10_p_value_descending"

        elif rank_mode == "directional":
            score_column1 = "-log10_p_value_descending"
            score_column2 = "-log10_p_value_ascending"
            
        # Extract and transpose rank scores for pop1_pop2 (pop1 under selection)
        left_df = df[["variant_pos", score_column1]].copy()
        left_df.rename(columns={score_column1: pop_pair}, inplace=True)
        left_df = left_df.set_index("variant_pos").T
        transp_list.append(left_df)

        # Extract and transpose rank scores for pop2_pop1 (pop2 under selection)
        reverse_pop_pair = "_".join(pop_pair.split("_")[::-1])
        right_df = df[["variant_pos", score_column2]].copy()
        right_df.rename(columns={score_column2: reverse_pop_pair}, inplace=True)
        right_df = right_df.set_index("variant_pos").T
        transp_list.append(right_df)

    # Concatenate all transposed DataFrames into a single DataFrame
    concat_df = pd.concat(transp_list, ignore_index=False)
    
    # Validate that data was found in the requested region
    if concat_df.empty or len(concat_df.columns) == 0:
        raise ValueError(
            f"No data found in the requested genomic region ({start:,} - {end:,}). "
            f"The data in '{input_dir}' contains positions from {_first_file_range[0]:,} to {_first_file_range[1]:,}. "
            f"Please specify a region that overlaps with this range using --start/--end."
        )
    
    # Get actual data range after filtering
    actual_start = int(concat_df.columns[0])
    actual_end = int(concat_df.columns[-1])
    num_variants = len(concat_df.columns)
    
    logger.debug(f"{'─' * 60}")
    logger.debug("Region Selection Summary")
    logger.debug(f"{'─' * 60}")
    logger.debug(f"  Requested region:    {start:,} - {end:,}")
    logger.debug(f"  Available data:      {_first_file_range[0]:,} - {_first_file_range[1]:,}")
    logger.debug(f"  Selected region:     {actual_start:,} - {actual_end:,}")
    logger.debug(f"  Variants in region:  {num_variants:,}")
    
    adjustments = []
    if start < _first_file_range[0]:
        adjustments.append(f"start adjusted from {start:,} to {actual_start:,} (requested start precedes available data)")
    elif actual_start != start:
        adjustments.append(f"start snapped to nearest variant at {actual_start:,}")
    
    if end > _first_file_range[1]:
        adjustments.append(f"end adjusted from {end:,} to {actual_end:,} (requested end exceeds available data)")
    elif actual_end != end:
        adjustments.append(f"end snapped to nearest variant at {actual_end:,}")
    
    if adjustments:
        logger.debug(f"  Note: Region boundaries adjusted to match available variant positions:")
        for adj in adjustments:
            logger.debug(f"    • {adj}")
    logger.debug(f"{'─' * 60}")
    
    pop_labels = concat_df.index.values
    first_pop = [pop.split("_")[0] for pop in pop_labels]
    second_pop = [pop.split("_")[1] for pop in pop_labels]
    concat_df["first_pop"] = first_pop
    concat_df["second_pop"] = second_pop

    concat_df.first_pop = concat_df.first_pop.astype("category")
    concat_df.first_pop = concat_df.first_pop.cat.set_categories(population_sorter)
    concat_df.second_pop = concat_df.second_pop.astype("category")
    concat_df.second_pop = concat_df.second_pop.cat.set_categories(population_sorter)

    # Clean up the DataFrame
    concat_df.sort_values(["first_pop", "second_pop"], inplace=True)
    concat_df.drop(["first_pop", "second_pop"], axis=1, inplace=True)
    concat_df.index.name = "pop_pairs"
    concat_df.attrs["populations"] = tuple(population_sorter)
    concat_df.attrs["population_mode"] = populations

    return concat_df

def summarize_by_superpopulation(input_df, populations="1000Genomes", agg_func='mean'):
    """
    Collapse population pair rows to superpopulation-level summaries.
    
    This function aggregates the data by superpopulation pairs, reducing the
    number of rows from n*(n-1) population pairs to m*(m-1) superpopulation pairs,
    where m is the number of superpopulations (5 for 1000 Genomes).
    
    Parameters
    ----------
    input_df : pandas.DataFrame
        Input DataFrame with population pairs as rows (index format: "POP1_POP2").
    populations : str or iterable, optional
        Population set. Default "1000Genomes" uses the standard 26 populations.
    agg_func : str or callable, optional
        Aggregation function: 'mean' (default), 'median', 'max', 'min', or callable.
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with superpopulation pairs as rows, aggregated values as data.
    """
    if populations != "1000Genomes":
        raise ValueError("summarize_by_superpopulation currently only supports 1000Genomes populations")
    
    # Add superpopulation columns based on row index
    df = input_df.copy()
    pop_pairs = df.index.tolist()
    
    superpop1 = []
    superpop2 = []
    for pair in pop_pairs:
        p1, p2 = pair.split("_")
        superpop1.append(pop_to_superpop.get(p1, "Unknown"))
        superpop2.append(pop_to_superpop.get(p2, "Unknown"))
    
    df['superpop1'] = superpop1
    df['superpop2'] = superpop2
    
    # Group by superpopulation pairs and aggregate
    grouped = df.groupby(['superpop1', 'superpop2']).agg(agg_func)
    
    # Create new index in format "SUPERPOP1_SUPERPOP2"
    grouped.index = [f"{sp1}_{sp2}" for sp1, sp2 in grouped.index]
    grouped.index.name = "superpop_pairs"
    
    return grouped


def extract_top_regions(input_df, n_top=100, window_size=10000, min_gap=5000):
    """
    Extract the top-scoring genomic regions from the heatmap data.
    
    This function identifies genomic windows with the highest aggregate scores
    across all population pairs, useful for automated discovery of selection signals.
    
    Parameters
    ----------
    input_df : pandas.DataFrame
        Input DataFrame with population pairs as rows and genomic positions as columns.
    n_top : int, optional
        Number of top regions to return (default: 100).
    window_size : int, optional
        Size of genomic windows in base pairs (default: 10000).
    min_gap : int, optional
        Minimum gap between reported regions to avoid overlapping hits (default: 5000).
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'start', 'end', 'center', 'mean_score', 'max_score',
        'top_population_pair', sorted by mean_score descending.
    """
    positions = input_df.columns.tolist()
    position_scores = input_df.mean(axis=0)
    sorted_positions = position_scores.sort_values(ascending=False)
    
    results = []
    used_positions = set()
    
    for pos in sorted_positions.index:
        # Skip if too close to already-selected region
        if any(abs(pos - used) < min_gap for used in used_positions):
            continue
            
        # Define window around this position
        window_start = pos - window_size // 2
        window_end = pos + window_size // 2
        
        # Get positions within window
        window_positions = [p for p in positions if window_start <= p <= window_end]
        if not window_positions:
            continue
            
        # Calculate statistics for this window
        window_data = input_df[window_positions]
        mean_score = window_data.values.mean()
        max_score = window_data.values.max()
        
        # Find top population pair in this window
        pair_means = window_data.mean(axis=1)
        top_pair = pair_means.idxmax()
        
        results.append({
            'start': min(window_positions),
            'end': max(window_positions),
            'center': pos,
            'mean_score': mean_score,
            'max_score': max_score,
            'top_population_pair': top_pair,
            'n_variants': len(window_positions)
        })
        
        used_positions.add(pos)
        
        if len(results) >= n_top:
            break
    
    return pd.DataFrame(results)


def plot_exp_heatmap(
    input_df,
    start,
    end,
    title=None,
    output="ExP_heatmap",
    output_suffix="png",
    cmap="Blues",
    populations="1000Genomes",
    vertical_line=True,
    cbar_vmin=None,
    cbar_vmax=None,
    xlabel=None,
    cbar_ticks=None,
    display_limit=None,
    display_values="higher",
    row_order=None,
    max_columns=None,
    column_aggregation="max",
    dpi=400,
):
    """
    Generate an ExP heatmap from a pandas DataFrame of pairwise statistics.

    Parameters
    ----------
    input_df : pandas.DataFrame
        Input data containing pairwise statistics or rank scores to visualize.
    start : int
        Start genomic position for the x-axis (region to display).
    end : int
        End genomic position for the x-axis (region to display).
    title : str, optional
        Title of the heatmap.
    output : str, optional
        Output filename (default: ExP_heatmap).
    output_suffix : str, optional
        File extension for the output (default: "png").
    cmap : str, optional
        Colormap for the heatmap (default Blues).
    populations : str or iterable, optional
        Population set to use. By default, expects "1000Genomes" (26 populations).
        For custom populations, provide an iterable of names. The input_df
        should contain all pairwise combinations in both directions.
    vertical_line : bool or list, optional
        If True, draws a vertical line at the center of the x-axis.
        If a list of (position, label) tuples is provided, draws multiple vertical lines with labels.
    cbar_vmin, cbar_vmax : float, optional
        Minimum and maximum values for the colorbar.
    xlabel : str, optional
        Custom x-axis label.
    cbar_ticks : list, optional
        Custom ticks for the colorbar.
    display_limit : float, optional
        If set, values above or below this threshold are set to zero, depending on `display_values`.
        Useful for highlighting only the most extreme results (e.g., top 1000 SNPs).
    display_values : {'higher', 'lower'}, optional
        Determines which values to retain when applying `display_limit`.
        Use "higher" for rank scores or distances (default), "lower" for metrics where smaller values are more extreme.
    row_order : list, optional
        Custom row order as a list of population pair names (e.g., ['ACB_ASW', 'ACB_BEB', ...]).
        If provided, rows will be reordered accordingly.
    max_columns : int or None, optional
        Maximum number of x-axis columns to render in the static figure. If None,
        an automatic budget is derived from the export width and DPI.
    column_aggregation : {'max', 'mean', 'median'}, optional
        Aggregation used when downsampling wide regions for static rendering.
    dpi : int, optional
        Output DPI for the saved figure.

    Returns
    -------
    matplotlib.axes.Axes
        The Matplotlib axis object. Also saves the heatmap to file.
    """

    # Helper functions to handle cases where start/end positions are not in the data
    def take_closest_start(sorted_list, target):
        """
        Given a sorted list, return the value equal to target if present,
        otherwise return the smallest value in the list that is >= target.
        If target is greater than all elements, raise ValueError.
        """
        if not sorted_list:
            raise ValueError("Input list is empty")
        pos = bisect_left(sorted_list, target)
        if pos < len(sorted_list) and sorted_list[pos] == target:
            return target
        if pos == len(sorted_list):
            raise ValueError("Your 'start' position index was higher than the range of the input data")
        return sorted_list[pos]

    def take_closest_end(sorted_list, target):
        """
        Given a sorted list, return the value equal to target if present,
        otherwise return the largest value in the list that is <= target.
        If target is less than all elements, raise ValueError.
        """
        if not sorted_list:
            raise ValueError("Input list is empty")
        pos = bisect_left(sorted_list, target)
        if pos < len(sorted_list) and sorted_list[pos] == target:
            return target
        if pos == 0:
            raise ValueError("Your 'end' position index was lower than the range of the input data")
        return sorted_list[pos - 1]
    
    input_df = input_df.copy()
    
    # Crop the input_df according to user defined range
    try: # Given values are in the data (column index)
        input_df = input_df.loc[:, start:end]
    except: # Given values are not in the column index, choose the new closest ones
        sorted_columns = sorted(list(input_df.columns))
        new_start = take_closest_start(sorted_columns, start)
        new_end = take_closest_end(sorted_columns, end)
        input_df = input_df.loc[:, new_start:new_end]
        
        logger.warning(f"Given 'start' and 'end' datapoints are not found in the input data. New closest datapoints were selected: start={new_start}, end={new_end}")
        
    # Check the input data for number of populations and input_df shape
    if populations == "1000Genomes":
        logger.debug("Input data: 1000 Genomes Project, phase 3 (26 populations expected)")
        logger.debug(f"Expecting 650 population pairwise comparisons: {input_df.shape[0]} found")
        if input_df.shape[0] != 650:
            raise ValueError(f"With selected populations='1000Genomes' option, the input_df was expected to have 650 rows, actual shape was: {input_df.shape[0]} rows, {input_df.shape[1]} columns")

    else:
        expected_pop_pairs = len(populations) * (len(populations) - 1)
        logger.debug(f"Input data: {len(populations)} custom populations")
        logger.debug(f"Expecting {expected_pop_pairs} population pairs: {input_df.shape[0]} found")
        if input_df.shape[0] != expected_pop_pairs:
            raise ValueError(f"With selected populations={populations} option, the input_df was expected to have {expected_pop_pairs} rows, actual shape was: {input_df.shape[0]} rows, {input_df.shape[1]} columns")

    # Apply the display_limit if selected
    if display_limit:
        if display_values == "higher":
            logger.debug(f"Displaying only values above the given display_limit: {display_limit}")
            input_df[input_df < display_limit] = 0
            
        elif display_values == "lower":
            logger.debug(f"Displaying only values below the given display_limit: {display_limit}")
            input_df[input_df > display_limit] = 0
            
        else:
            raise ValueError(f"plot_exp_heatmap() parameter 'display_values' has unknown value '{display_values}'. The only expected options are 'higher' or 'lower'.")
    
    # Apply custom row ordering
    if row_order is not None:
        try:
            input_df = input_df.loc[row_order]
            logger.debug(f"Applied custom row ordering with {len(row_order)} rows")
        except KeyError as e:
            raise ValueError(f"Some row names in row_order not found in data: {e}")
    
    # Custom expheatmap colormap
    if cmap == "expheatmap":
        from matplotlib import cm
        import matplotlib as mpl

        cmap = cm.gist_ncar_r(np.arange(256))
        cmap[0] = [1.0, 1.0, 1.0, 1.0]
        cmap[-1] = [0, 0, 0.302, 1.0]
        cmap = mpl.colors.ListedColormap(cmap, name='expheatmap_cmap', N=cmap.shape[0])

    # Determine if we're using 1000 Genomes populations
    is_1000genomes = populations == "1000Genomes"
    if is_1000genomes:
        populations = populations_1000genomes

    original_start = int(input_df.columns[0])
    original_end = int(input_df.columns[-1])

    # Calculate figure size
    if is_1000genomes:
        figsize = (15, 5)
    else:
        base_width = 15 + (input_df.shape[1] // 1000)
        base_height = 5 + (input_df.shape[0] // 900)
        figsize = (base_width, base_height)

    original_column_count = input_df.shape[1]
    if max_columns is None:
        max_columns = auto_max_columns(figsize, dpi)
    input_df, was_downsampled = downsample_heatmap_columns(
        input_df,
        max_columns=max_columns,
        aggregation=column_aggregation,
    )
    if was_downsampled:
        logger.info(
            "Static heatmap downsampled from "
            f"{original_column_count} to {input_df.shape[1]} columns "
            f"using {column_aggregation} aggregation"
        )

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot the main heatmap
    cbar_kws_base = {"shrink": 0.9, "aspect": 30}
    if is_1000genomes:
        cbar_kws = {**cbar_kws_base, "ticks": cbar_ticks} if cbar_ticks else cbar_kws_base
        hm = sns.heatmap(
            input_df,
            yticklabels=False,
            xticklabels=False,
            vmin=1.301 if cbar_vmin is None else cbar_vmin,
            vmax=4.833 if cbar_vmax is None else cbar_vmax,
            ax=ax,
            cmap=cmap,
            cbar_kws=cbar_kws,
        )
    else:
        cbar_kws = {**cbar_kws_base, "ticks": cbar_ticks}
        hm = sns.heatmap(
            input_df,
            yticklabels=False,
            xticklabels=False,
            vmin=cbar_vmin,
            vmax=cbar_vmax,
            cbar_kws=cbar_kws,
            ax=ax,
            cmap=cmap,
        )
    
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(length=0, labelsize=7)

    ax.set_title(title)
    
    # Set x-axis ticks and labels
    x_extent = input_df.shape[1]
    actual_middle = (original_start + original_end) // 2
    ax.set_xticks([0, x_extent // 2, x_extent - 1])
    ax.set_xticklabels([f"{original_start:,}", f"{actual_middle:,}", f"{original_end:,}"], fontsize=7)
    ax.tick_params(axis="x", length=0, pad=8)
    ax.set_xlabel(xlabel if xlabel else "")
        
    # Set the y-axis ticks and labels
    y_axis_len = len(populations) * (len(populations) - 1)
    n_pops = len(populations)
    ax.set_yticks([])
    label_positions = np.arange((n_pops - 1) / 2, y_axis_len, step=(n_pops - 1))
    ax.set_yticks(label_positions, minor=True)
    ax.tick_params(axis="y", which="minor", length=0)
    fontsize = max(4, min(8, 200 // n_pops))
    ax.set_yticklabels(populations, minor=True, fontsize=fontsize)
    ax.set_ylabel('')
    
    # Add superpopulation labels for 1k Genomes data
    if is_1000genomes:
        superpop_order = ["AFR", "SAS", "EAS", "EUR", "AMR"]
        superpop_sizes = [len(superpopulations[sp]) for sp in superpop_order]

        x_extent = input_df.shape[1]
        
        cumulative_pos = 0
        for i, (sp, sp_size) in enumerate(zip(superpop_order, superpop_sizes)):
            sp_rows = sp_size * (n_pops - 1)
            center_pos = cumulative_pos + sp_rows / 2
            
            ax.text(x_extent * 1.005, center_pos, sp, 
                    ha='left', va='center', fontsize=7,
                    rotation=270, clip_on=False)
            
            if i < len(superpop_order) - 1:
                boundary_pos = cumulative_pos + sp_rows
                ax.axhline(y=boundary_pos, color='#CCCCCC', linewidth=0.3, zorder=1)
            
            cumulative_pos += sp_rows

    # Add vertical line
    if vertical_line is True:
        middle = int(input_df.shape[1] / 2)
        ax.axvline(x=middle, linewidth=0.3, color='#CCCCCC', zorder=1)
        
    elif vertical_line and hasattr(vertical_line, '__iter__') and not isinstance(vertical_line, str):
        list_of_columns = input_df.columns.to_list()
        positions_indices = []
        labels = []
        
        for item in vertical_line:
            if not isinstance(item, (tuple, list)) or len(item) != 2:
                raise ValueError(f"Each vertical line item must be a (position, label) tuple, got: {item}")
            
            pos, label = item
            try:
                col_index = list_of_columns.index(pos)
            except ValueError:
                if not list_of_columns:
                    logger.warning(f"Position {pos} not found in data columns, skipping this vertical line")
                    continue
                col_index = min(range(len(list_of_columns)), key=lambda idx: abs(list_of_columns[idx] - pos))
                logger.warning(
                    f"Position {pos} not found exactly after binning; using nearest rendered position "
                    f"{list_of_columns[col_index]}"
                )

            positions_indices.append(col_index)
            labels.append(label)
            ax.axvline(x=col_index, linewidth=0.3, color='#CCCCCC', zorder=1)
        
        if positions_indices:
            ax.set_xticks(positions_indices)
            ax.set_xticklabels(labels)
    
    ax.figure.savefig(f"{output}.{output_suffix}", dpi=dpi, bbox_inches="tight")
    logger.info(f"ExP heatmap saved to {output}.{output_suffix}")
    return ax


def prepare_cbar_params(data_df, n_cbar_ticks=4):
    """
    Calculate optimal colorbar parameters for heatmap visualization.
    
    This function analyzes the data range and automatically determines appropriate
    minimum and maximum values for the colorbar, along with evenly spaced tick marks.
    The colorbar bounds are intelligently chosen to provide good visual contrast
    while encompassing the full data range.
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        Input DataFrame containing the data values to be visualized in the heatmap.
        All numeric values in the DataFrame are considered for range calculation.
    n_cbar_ticks : int, optional
        Number of tick marks to display on the colorbar (default: 4).
        Must be at least 2 to show minimum and maximum values.
    
    Returns
    -------
    tuple of (float, float, list)
        A tuple containing:
        - cbar_vmin (float): Minimum value for the colorbar scale
        - cbar_vmax (float): Maximum value for the colorbar scale  
        - cbar_ticks (list): List of evenly spaced tick values for the colorbar
   
    """
    from math import floor
    
    cbar_min = 0
    cbar_max = 0
    data_min = data_df.min().min()
    data_max = data_df.max().max()
    
    # Determine the colorbar minimum and maximum bounds
    if data_min < 0.5:
        cbar_min = 0    
    elif data_min < 1:
        cbar_min = 0.5
    else:
        cbar_min = floor(data_min)
        
    if data_max < 1:
        cbar_max = 1
    else:
        cbar_max = floor(data_max + 1)
        
    # Generate evenly spaced tick marks for the colorbar
    cbar_ticks = np.arange(cbar_min, cbar_max + 0.001, step=(cbar_max - cbar_min)/(n_cbar_ticks-1))
    
    return cbar_min, cbar_max, list(cbar_ticks)



def plot(
    input_dir,
    start,
    end,
    title,
    output="ExP_heatmap",
    cmap="Blues",
    max_columns=None,
    column_aggregation="max",
    dpi=400,
    rank_scores="directional",
):
    """
    Generate and save an ExP heatmap from XP-EHH analysis results.

    This function serves as a high-level wrapper that processes a directory of 
    XP-EHH (Cross-Population Extended Haplotype Homozygosity) results and creates 
    a publication-ready heatmap visualization. It combines data loading, processing, 
    and plotting into a single convenient function call.

    Parameters
    ----------
    input_dir : str
        Path to directory containing XP-EHH results as .tsv files.
        Files should be named in the format 'POP1_POP2.*.tsv' where POP1 and POP2 
        are population identifiers from the 1000 Genomes Project.
    start : int
        Start genomic position (inclusive) for the region to visualize.
        Positions are typically in base pairs along a chromosome.
    end : int
        End genomic position (inclusive) for the region to visualize.
        Must be greater than start position.
    title : str
        Main title for the heatmap plot.
    output : str
        Output filename (without extension) for saving the heatmap.
        The plot will be saved as '{output}.png' by default.
    cmap : str, optional
        Matplotlib colormap name for the heatmap visualization (default: "Blues").
        Can also use "expheatmap" for a custom colormap optimized for this analysis.
    max_columns : int or None, optional
        Maximum number of rendered x-axis columns. If None, choose automatically
        from figure width and DPI.
    column_aggregation : {'max', 'mean', 'median'}, optional
        Aggregation used when downsampling static output.
    dpi : int, optional
        DPI used when saving the static figure.
    rank_scores : {'directional', '2-tailed', 'ascending', 'descending'}, optional
        Which empirical rank-score mode to visualize.

    Returns
    -------
    matplotlib.axes.Axes
        The Matplotlib axis object.

    Notes
    -----
    This function uses 1000 Genomes-specific layout only when the full canonical
    26-population set is detected. Other datasets are treated as custom population
    panels inferred from the compute output filenames.
    """
    try:
        plot_input = create_plot_input(input_dir, start=start, end=end, rank_scores=rank_scores)
        populations = plot_input.attrs.get("population_mode", "1000Genomes")
        
        # Use actual data range for xlabel (not requested range)
        actual_start = int(plot_input.columns[0])
        actual_end = int(plot_input.columns[-1])
        
        return plot_exp_heatmap(
            plot_input, start=actual_start,
            end=actual_end,
            title=title,
            cmap=cmap,
            output=output,
            populations=populations,
            max_columns=max_columns,
            column_aggregation=column_aggregation,
            dpi=dpi,
        )
    except KeyboardInterrupt:
        logger.info("")
        sys.exit(1)
