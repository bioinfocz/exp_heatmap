import os
import sys

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from bisect import bisect_left


# 1000 Genomes populations
populations_1000genomes = (
    "ACB",  # AFR (superpopulations)
    "ASW",
    "ESN",
    "GWD",
    "LWK",
    "MSL",
    "YRI",
    "BEB",  # SAS
    "GIH",
    "ITU",
    "PJL",
    "STU",
    "CDX",  # EAS
    "CHB",
    "CHS",
    "JPT",
    "KHV",
    "CEU",  # EUR
    "FIN",
    "GBR",
    "IBS",
    "TSI",
    "CLM",  # AMR
    "MXL",
    "PEL",
    "PUR",
)

# 1000 Genomes super-populations
superpopulations = {"AFR": ("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"),
                    "SAS": ("BEB", "GIH", "ITU", "PJL", "STU"),
                    "EAS": ("CDX", "CHB", "CHS", "JPT", "KHV"),
                    "EUR": ("CEU", "FIN", "GBR", "IBS", "TSI"),
                    "AMR": ("CLM", "MXL", "PEL", "PUR",)}


def create_plot_input(input_dir, start, end, populations="1000Genomes", rank_pvalues="2-tailed"):
    """
    Generate a pandas DataFrame for plotting from a directory of pairwise population .tsv files.

    Parameters:
        input_dir (str): Directory containing .tsv files named as 'POP1_POP2.*.tsv' for each population pair.
        start (int): Start genomic position (X-axis lower bound).
        end (int): End genomic position (X-axis upper bound).
        populations (iterable or "1000Genomes"): List of population names to include and order in the heatmap. If not using 1000 Genomes, provide a custom iterable.
        rank_pvalues (str): Determines which p-value column to use:
            - "ascending": Use "-log10_p_value_ascending".
            - "descending": Use "-log10_p_value_descending" (default).
            - "2-tailed": For pop1_pop2, use ascending; for pop2_pop1, use descending.

    Returns:
        pd.DataFrame: Formatted data for heatmap plotting.
    """
    
    
    df_list = []
    pop_id_list = []
    different_dfs = False


    # if not using 1000Genomes, use the custom populations' list to sort the data
    if populations=="1000Genomes":
        population_sorter = populations_1000genomes
        
    else:
        population_sorter = populations
    
    
    
    segment_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    

    ##############################################

    # reading the input files, saving only the regions between START and END to process further
    index = 1
    for segment_file in segment_files:
        # segment_files is something like ACB_KHV.tsv or ACB_KHV.some_more_info.tsv
        pop_pair = os.path.splitext(os.path.basename(segment_file))[0].split(".")[0]
        
        # test, if file names (p1, p2) are in the 'populations' a.k.a 'population_sorter'
        p1, p2 = pop_pair.split("_")
        
        if all(x in population_sorter for x in (p1,p2)):
            pop_id_list.append(pop_pair)

            print(
                "[{}/{}] Loading {} from {}".format(
                    index, len(segment_files), pop_pair, segment_file
                )
            )

            segments = pd.read_csv(segment_file, sep="\t")
            segments = segments[
                (segments.variant_pos >= start) & (segments.variant_pos <= end)
            ]

            df_list.append(segments)

            index += 1
            
        else:
            print(
                "[{}/{}] ERROR Loading {} from {}. {} or {} not in provided 'populations' list.".format(
                    index, len(segment_files), pop_pair, segment_file, p1, p2
                )
            )
                
            return p1, p2, population_sorter
            

    # check that they all have the same dimensions AND variant_pos
    df_shape = df_list[0].shape
    variant_positions = df_list[0].variant_pos.values

    print("Transforming data matrix in preparation to plot heatmap")

    for i in range(len(df_list)):
        if df_list[i].shape != df_shape:
            print("the shapes dont match in df " + str(i))
            different_dfs = True
            break
        if not np.array_equal(df_list[i].variant_pos.values, variant_positions):
            print("the variant_positions dont match in df " + str(i))
            different_dfs = True
            break

    if different_dfs:
        sys.exit(1)

    # select only variant_pos and -log10_p_value and transpose each df
    transp_list = []

    for df, pop_pair in zip(df_list, pop_id_list):
        if rank_pvalues == "ascending":
            pvalues_column1 = "-log10_p_value_ascending"
            pvalues_column2 = "-log10_p_value_ascending"
        
        elif rank_pvalues == "descending":
            pvalues_column1 = "-log10_p_value_descending"
            pvalues_column2 = "-log10_p_value_descending"

        elif rank_pvalues == "2-tailed":
            pvalues_column1 = "-log10_p_value_descending"
            pvalues_column2 = "-log10_p_value_ascending"
        
        else:
            raise ValueError(f"Unknown value for 'rank_pvalues' parameter in create_plot_input(). Expected values are: 'ascending', 'descending' or '2-tailed', got '{rank_pvalues}'")
            
        
        
        # select the appropriate ranks that are significant for pop1_pop2 (pop1 is under selection)
        left_df = df[["variant_pos", pvalues_column1]].copy()
        left_df.rename(columns={pvalues_column1: pop_pair}, inplace=True)
        left_df = left_df.set_index("variant_pos").T
        transp_list.append(left_df)

        # select the apropriate ranks that are significant for pop2_pop1 (pop2 is under selection)
        reverse_pop_pair = "_".join(
            pop_pair.split("_")[::-1]
        )  # change name pop1_pop2 to pop2_pop1

        right_df = df[["variant_pos", pvalues_column2]].copy()
        right_df.rename(
            columns={pvalues_column2: reverse_pop_pair}, inplace=True
        )
        right_df = right_df.set_index("variant_pos").T
        transp_list.append(right_df)

    # concatenate all the dfs together
    big_df = pd.concat(transp_list, ignore_index=False)

    print("Sorting data by super populations")

    # add temporary columns with pop1 and pop2, I am gonna sort the df according to those
    pop_labels = big_df.index.values  # select the pop1_pop2 names

    first_pop = [pop.split("_")[0] for pop in pop_labels]  # pop1
    second_pop = [pop.split("_")[1] for pop in pop_labels]  # pop2

    big_df["first_pop"] = first_pop
    big_df["second_pop"] = second_pop

    # set pop1 to be a categorical column with value order defined by sorter
    big_df.first_pop = big_df.first_pop.astype("category")
    big_df.first_pop = big_df.first_pop.cat.set_categories(population_sorter)

    # set pop2 to be a categorical column with value order defined by sorter
    big_df.second_pop = big_df.second_pop.astype("category")
    big_df.second_pop = big_df.second_pop.cat.set_categories(population_sorter)

    # sort df by pop1 and withing pop1 by pop2
    big_df.sort_values(["first_pop", "second_pop"], inplace=True)

    # drop the temporary columns
    big_df.drop(["first_pop", "second_pop"], axis=1, inplace=True)
    
    # set index name
    big_df.index.name = "pop_pairs"
    
    # label it just by the pop1 (which is gonna be printed with plot ticks)
    #pop_labels = big_df.index.values
    #pop_labels = [pop.split("_")[0] for pop in pop_labels]
    #big_df.index = pop_labels

    return big_df



def plot_exp_heatmap(
    input_df,
    start,
    end,
    title,
    output=None,
    output_suffix="png",
    cmap=None,
    populations="1000Genomes",
    vertical_line=True,
    cbar_vmin=None,
    cbar_vmax=None,
    ylabel=False,
    xlabel=False,
    cbar_ticks=None,
    display_limit=None,
    display_values="higher"
):
    """
    Generate an ExP heatmap from a pandas DataFrame of pairwise statistics.

    Parameters
    ----------
    input_df : pandas.DataFrame
        Input data containing pairwise statistics or p-values to visualize.
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
    ylabel, xlabel : str, optional
        Custom y-axis/x-axis labels.
    cbar_ticks : list, optional
        Custom ticks for the colorbar.
    display_limit : float, optional
        If set, values above or below this threshold are set to zero, depending on `display_values`.
        Useful for highlighting only the most significant results (e.g., top 1000 SNPs).
    display_values : {'higher', 'lower'}, optional
        Determines which values to retain when applying `display_limit`.
        Use "higher" for rank p-values or distances (default), "lower" for classical p-values.

    Returns
    -------
    The Matplotlib axis object and saves the heatmap to file
    """

    
    # functions to solve the situation where given start and end indexes are not in the data
    def take_closest_start(myList, myNumber):
        """
        Given a sorted list, return the value equal to target if present,
        otherwise return the smallest value in the list that is >= target.
        If target is greater than all elements, raise ValueError.
        """

        if myNumber in myList:
            return myNumber

        pos = bisect_left(myList, myNumber)
        if pos == 0:
            return myList[0]
        if pos == len(myList):
            raise ValueError("Your 'start' position index was higher than the range of the input data")

        after = myList[pos]

        return after


    def take_closest_end(myList, myNumber):
        """
        Given a sorted list, return the value equal to target if present,
        otherwise return the largest value in the list that is <= target.
        If target is less than all elements, raise ValueError.
        """

        if myNumber in myList:
            return myNumber

        pos = bisect_left(myList, myNumber)
        if pos == 0:
            raise ValueError("Your 'end' position index was lower than the range of the input data")
        if pos == len(myList):
            return myList[-1]

        before = myList[pos - 1]

        return before
    
    
    
    input_df = input_df.copy()
    
    print("Checking input")
    
    # cropping the input_df according to user defined range
    try: # given values are in the data (column index)
        input_df = input_df.loc[:, start:end]
    
    except: # given values are not in the column index, choose the new closest ones
        sorted_columns = sorted(list(input_df.columns)) # sort columns (just to be sure)

        new_start = take_closest_start(sorted_columns, start) # take the closest value to the left from given start point
        new_end = take_closest_end(sorted_columns, end) # take the closest values the the right from given end point
        
        
        input_df = input_df.loc[:, new_start:new_end]
        
        print(f"WARNING: Given 'start' and 'end' datapoints are not found in the input data. New closest datapoints were selected.\nstart={new_start}\nend={new_end}")
        
        
    # check the input data for number of populations and input_df shape
    if populations == "1000Genomes":

        print("- expecting 1000 Genomes Project, phase 3 input data, 26 populations")
        print("- expecting 650 population pairs...", end="")

        if input_df.shape[0] == 650:
            print("CHECK\n")

        else:
            print("ERROR")
            raise ValueError(
                "With selected populations='1000Genomes' option, the input_df was expected to have 650 rows, actual shape was: {} rows, {} columns".format(
                    input_df.shape[0], input_df.shape[1]
                )
            )

    else:

        n_populations = len(populations)
        print("- custom {} populations entered:".format(str(n_populations)))
        for i in populations:
            print(i, end="  ")

        print()
        print(
            "- expecting {} population pairs...".format(
                str(n_populations * (n_populations - 1))
            ),
            end="",
        )

        if input_df.shape[0] == (n_populations * (n_populations - 1)):
            print("CHECK\n")

        else:
            print("ERROR")
            raise ValueError(
                "With selected populations={} option, the input_df was expected to have {} rows, actual shape was: {} rows, {} columns".format(
                    populations,
                    n_populations * (n_populations - 1),
                    input_df.shape[0],
                    input_df.shape[1],
                )
            )

    
    # Apply the display_limit
    if display_limit:
        if display_values == "higher":
            print()
            print(f"Displaying only values above the given display_limit: {display_limit}")
            print()
            
            input_df[input_df < display_limit] = 0

        elif display_values == "lower":
            print()
            print(f"Displaying only values below the given display_limit: {display_limit}")
            print()
            
            input_df[input_df > display_limit] = 0

        else:
            raise ValueError(f"plot_exp_heatmap() parameter 'display_values' has unknown value '{display_values}'. The only expected options are 'higher' or 'lower'.")
    
    

    ########################
    # Color map definition #
    ########################

    # custom colormap assembly
    if cmap == "expheatmap":
        
        from matplotlib import cm
        import matplotlib as mpl

        cmap = cm.gist_ncar_r(np.arange(256))  # just a np array from cmap
        cmap[0] = [1.0, 1.0, 1.0, 1.0]  # change the lowest values in colormap to white
        cmap[-1] = [0, 0, 0.302, 1.0]
        # create cmap object from a list of colors (RGB)
        cmap = mpl.colors.ListedColormap(cmap, name='expheatmap_cmap',
                                             N=cmap.shape[0])

    # default colormap
    if not cmap:
        cmap = "Blues"

    

    #########################
    # create the ExP figure #
    #########################
    print("Creating heatmap")

    # draw default exp heatmap with 26 populations from 1000 Genomes Project
    if populations == "1000Genomes":
        populations = populations_1000genomes
        
        fig, ax = plt.subplots(figsize=(15, 5))
        
        sns.heatmap(
            input_df,
            yticklabels=populations,
            xticklabels=False,
            vmin=1 if cbar_vmin==None else cbar_vmin,
            vmax=4.853 if cbar_vmax==None else cbar_vmax,
            #cbar_kws={"ticks": [1.3, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]},
            ax=ax,
            cmap=cmap,
        )

        if not title:
            title = "{} - {}".format(start, end)

        if not output:
            output = title

        ax.set_title(title)
        ax.set_ylabel(
            "population pairings\n\n    AMR   |     EUR    |     EAS    |    SAS     |       AFR          "
        )

        # custom or default x-axis label
        if xlabel:
            ax.set_xlabel(xlabel)
        else:
            ax.set_xlabel("{:,} - {:,}".format(start, end))

    # draw custom exp heatmap with user-defined populations (number of pops, labels)
    else:
        
        # need to set up figure size for large population sample-sets
        # here the default size (15,5) is increased by 1 inch for every 1000 SNPs (x axis) and by 1 inch for every 900 population pairs (y axis) in input_df
        fig, ax = plt.subplots(figsize=(15 + (input_df.shape[1] // 1000), 5 + (input_df.shape[0] // 900)))
        
        sns.heatmap(
            input_df,
            yticklabels=populations,
            xticklabels=False,
            vmin=cbar_vmin,
            vmax=cbar_vmax,
            cbar_kws={"ticks": cbar_ticks},
            ax=ax,
            cmap=cmap,
        )

        if not title:
            title = "{} - {}".format(start, end)

        if not output:
            output = title

        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        
    # set the y-axis tics and labels
    y_axis_len = len(populations) * (len(populations) - 1)  # length of the input data, number of df row
    y_labels_pos = list(np.arange(0, y_axis_len, step=(len(populations)-1)))  # arange positions in steps
    y_labels_pos.append(y_axis_len) # add also the last one
    ax.set_yticks(y_labels_pos)
    ax.set_yticks(np.arange(y_labels_pos[0] + ((len(populations)-1) / 2), y_labels_pos[-1], step=(len(populations)-1)),
                  minor=True) # position of minor yticks - just between the major one
    ax.tick_params(axis="y", which="minor", length=0) # set minor yaxis ticks to zero length

    ax.set_yticklabels(populations, minor=True)

    

    # optionally add vertical line in the middle of the figure, or else as defined
    try:  # test if vertical line is iterable --> draw more vertical lines
        iterator = iter(vertical_line)

        # vertical_line=([pos1, label1], [pos2, label2], [pos3, label3]...)

        try:
            list_of_columns = input_df.columns.to_list()
            
            for (pos, label) in vertical_line:
                ax.axvline(x=list_of_columns.index(pos), label=label, linewidth=1, color="grey")
                
            # I will get the list of positions, but the heatmap x-axis is indexed from 0
            # so I need to turn the positions (non-consecutive, as they are SNPs!!) into indices of column labels
            ax.set_xticks([list_of_columns.index(i[0]) for i in vertical_line])  # what column index is the user-defined x position of the vline?
            ax.set_xticklabels([i[1] for i in vertical_line])  # labels

        except:
            print("Could not read 'vertical_line', was expecting this 'vertical_line=([x1, label1], [x2, label2], [x3, label3]...)'.")
            print("Vertical lines might be out of range of displayed graph are ('start', 'end'), please double-check")
            print(f"Got this input for 'vertical_line': {vertical_line}")
            print("---")
            print("No vertical line will be displayed")
            print()



    
    except TypeError:
        # not iterable
        if vertical_line: # just one verticle line in the middle

            middle = int(input_df.shape[1] / 2)
            ax.axvline(x=middle, linewidth=1, color="grey")

        else:
            print("Could not read 'vertical_line', was expecting this 'vertical_line=([x1, label1], [x2, label2], [x3, label3]...)'.")
            print("No vertical line will be displayed")
            print()


    

    print("Savig heatmap")


    if output:
        print()
        print(f"ExP heatmap saved into {output}.{output_suffix}")
        
        ax.figure.savefig(f"{output}.{output_suffix}", dpi=400, bbox_inches="tight")
        
    else:
        plt.show()
        
    
    #plt.close(fig)

    
    
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
    
    import numpy as np
    from math import floor
    
    
    # target min max cbar values
    cbar_min = 0
    cbar_max = 0
    
    # min max values in data
    data_min = data_df.min().min()
    data_max = data_df.max().max()
    
    
    # deciding min cbar values
    if data_min < 0.5:
        cbar_min = 0
        
    elif data_min < 1:
        cbar_min = 0.5
        
    else:
        cbar_min = floor(data_min)
        
        
    # deciding max cbar values
    if data_max < 1:
        cbar_max = 1
        
    else:
        cbar_max = floor(data_max + 1)
        
    # cbar ticks; adjust the cmax value by minimal amount to get the clear cmax value from np.arange function
    cbar_ticks = np.arange(cbar_min, cbar_max + 0.001, step=(cbar_max - cbar_min)/(n_cbar_ticks-1))
    
    return cbar_min, cbar_max, list(cbar_ticks)



def plot(input_dir, start, end, title, output, cmap="Blues"):
    """
    Generate and save an ExP heatmap.

    This function serves as a high-level wrapper that processes a directory of 
    XP-EHH (Cross-Population Extended Haplotype Homozygosity) results and creates 
    a publication-ready heatmap visualization. It combines data loading, processing, 
    and plotting into a single convenient function call.

    Parameters
    ----------
    input_dir : str
        Path to directory containing XP-EHH results as .tsv files.
        Files should be named in the format 'POP1_POP2.*.tsv'
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

    Returns
    -------
    None
        The function saves the heatmap.

    Notes
    -----
    This function assumes 1000 Genomes Project population structure and expects
    650 pairwise population comparisons. The function will automatically handle
    cases where the exact start/end positions are not present in the data by
    selecting the closest available positions.
    """
    try:
        plot_input = create_plot_input(input_dir, start=start, end=end)
        
        plot_exp_heatmap(
            plot_input, start=plot_input.columns[0],
            end=plot_input.columns[-1],
            title=title,
            cmap=cmap,
            output=output,
            xlabel="{:,} - {:,}".format(start, end)
        )
    except KeyboardInterrupt:
        print("")
        sys.exit(1)