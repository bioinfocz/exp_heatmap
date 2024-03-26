import os
import sys

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob


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


def create_plot_input(input_dir, begin, end, populations="1000Genomes", rank_pvalues="descending"):
    """
    This function creates an input pandas DataFrame that will subsequently be used as input for ExP heatmap plotting function

    input_dir -- directory, where the function expects the data to be found in a serie of *.tsv files for each population pair,
                 the names of these files should start with 'POP1_POP2' string

    begin, end -- limit the X-axis (position), displayed area

    populations -- by default, the '1000Genomes' option, 26 populations from 1000 Genomes Project are expected. If your population/classes set is different,
                   provide an iterable with population/classes names. For example, with populations=['pop1', 'pop2', 'pop3'] this function will search the input directory (input_dir)
                   for *.tsv files with all combinations of pairwise values, named pop1_pop2.*.tsv, pop1_pop3.*.tsv, pop2_pop1.*.tsv etc.
                   
                   Populations in the final image output will be sorted according the their order in this iterable ('populations')

                   The output 'big_df' DataFrame with following 6 rows of pairwise values will be created:
                   pop1 vs pop2
                   pop1 vs pop3
                   pop2 vs pop1
                   pop2 vs pop3
                   pop3 vs pop1
                   pop3 vs pop2

    rank_pvalues -- which test results/pvalues are the most significant. In selection tests, usually the highest test values are the most interesting ones.
                    In this situation the results' ranks for rank pvalues are computed in descending order (highest first, lowest ranks).
                    Possible values: "ascending" -> "-log10_p_value_ascending" column will be used to create the ExP heatmap
                                     "descending" -> "-log10_p_value_descending" column will be used to create the Exp heatmap (default)
    
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

    # reading the input files, saving only the regions between BEGIN and END to process further
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
                (segments.variant_pos >= begin) & (segments.variant_pos <= end)
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
            pvalues_column = "-log10_p_value_ascending"
        
        elif rank_pvalues == "descending":
            pvalues_column = "-log10_p_value_descending"
        
        else:
            raise ValueError(f"Unknown value for 'rank_pvalues' parameter in create_plot_input(). Expected values are: 'ascending' or 'descending', got '{rank_pvalues}'")
            
        
        
        # select the appropriate ranks that are significant for pop1_pop2 (pop1 is under selection)
        left_df = df[["variant_pos", pvalues_column]].copy()
        left_df.rename(columns={pvalues_column: pop_pair}, inplace=True)
        left_df = left_df.set_index("variant_pos").T
        transp_list.append(left_df)

        # select the apropriate ranks that are significant for pop2_pop1 (pop2 is under selection)
        reverse_pop_pair = "_".join(
            pop_pair.split("_")[::-1]
        )  # change name pop1_pop2 to pop2_pop1

        right_df = df[["variant_pos", pvalues_column]].copy()
        right_df.rename(
            columns={pvalues_column: reverse_pop_pair}, inplace=True
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
    begin,
    end,
    title,
    output,
    cmap=None,
    populations="1000Genomes",
    vertical_line=True,
    cbar_vmin=None,
    cbar_vmax=None,
    ylabel=False,
    xlabel=False,
    cbar_ticks=None,
):
    """
    Read input DataFrame and create the ExP heatmap accordingly.
    input_df -- input pandas DataFrame, data to display
    begin, end -- limit the X-axis (position), displayed area
    title -- title of the graph
    cmap -- seaborn colormap, 'Blues' or 'crest' work well
    vertical_line -- if True, displays one line in the middle of the x-axis
                     if vertical_line=([x1, label1], [x2, label2], [x3, label3]...), display multiple vlines with different labels
    output -- ouput file, will be saved with *.png suffix
    populations -- by default, the '1000Genomes' option, 26 populations from 1000 Genomes Project are expected. If your population/classes set is different,
                   provide an iterable with population/classes names. For example, with populations=['pop1', 'pop2', 'pop3'] this function expects input_df
                   with 6 rows of pairwise values:
                   pop1 vs pop2
                   pop1 vs pop3
                   pop2 vs pop1
                   pop2 vs pop3
                   pop3 vs pop1
                   pop3 vs pop2
    """

    print("Checking input")
    
    # cropping the input_df according to user defined range
    input_df = input_df.loc[:, begin:end]
    
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
            vmin=1,
            vmax=4.853,
            cbar_kws={"ticks": [1.3, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]},
            ax=ax,
            cmap=cmap,
        )

        if not title:
            title = "{} - {}".format(begin, end)

        if not output:
            output = title

        ax.set_title(title)
        ax.set_ylabel(
            "population pairings\n\n    AMR   |     EUR    |     EAS    |    SAS     |       AFR          "
        )
        ax.set_xlabel("{:,} - {:,}".format(begin, end))

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
            vmax=cbar_vmin,
            cbar_kws={"ticks": cbar_ticks},
            ax=ax,
            cmap=cmap,
        )

        if not title:
            title = "{} - {}".format(begin, end)

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
            print("Vertical lines might be out of range of displayed graph are ('begin', 'end'), please double-check")
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
        print(f"ExP heatmap saved into {output}.png")
        
        ax.figure.savefig(output, dpi=400, bbox_inches="tight")
        
    else:
        plt.show()
        
    
    #plt.close(fig)

    
    
    return ax
    

    
def prepare_cbar_params(data_df, n_cbar_ticks=4):
    """
    Gets the pandas.DataFrame (i.e. data_to_plot) with the data you want to display with your custom plot_exp_heatmap function,
    returns cbar_vmin, cbar_vmax and cbar_ticks.
    """
    
    import numpy as np
    from math import floor
    
    
    # target min max cbar values
    cmin = 0
    cmax = 0
    
    # min max values in data
    dmin = data_df.min().min()
    dmax = data_df.max().max()
    
    
    # deciding min cbar values
    if dmin < 0.5:
        cmin = 0
        
    elif dmin < 1:
        cmin = 0.5
        
    else:
        cmin = floor(dmin)
        
        
    # deciding max cbar values
    if dmax < 1:
        cmax = 1
        
    else:
        cmax = floor(dmax + 1)
        
    # cbar ticks; adjust the cmax value by minimal amount to get the clear cmax value from np.arange function
    cbar_ticks = np.arange(cmin, cmax + 0.001, step=(cmax - cmin)/(n_cbar_ticks-1))
    
    return cmin, cmax, list(cbar_ticks)



def plot(xpehh_dir, begin, end, title, output, cmap="Blues"):
    """
    create the plot function input data and print/save them
    """

    data_to_plot = create_plot_input(xpehh_dir, begin=begin, end=end)
    
    plot_exp_heatmap(
        data_to_plot, begin=data_to_plot_example.columns[0], end=data_to_plot_example.columns[-1], title=title, cmap=cmap, output=output
    )
