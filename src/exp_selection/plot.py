import os
import sys

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

population_sorter = (
    "ACB",
    "ASW",
    "ESN",
    "GWD",
    "LWK",
    "MSL",
    "YRI",  # AFR
    "BEB",
    "GIH",
    "ITU",
    "PJL",
    "STU",  # SAS
    "CDX",
    "CHB",
    "CHS",
    "JPT",
    "KHV",  # EAS
    "CEU",
    "FIN",
    "GBR",
    "IBS",
    "TSI",  # EUR
    "CLM",
    "MXL",
    "PEL",
    "PUR",
)  # AMR


def create_plot_input(input_dir, begin, end, populations="1000Genomes"):
    """
    This function creates an input pandas DataFrame that will subsequently be used as input for ExP heatmap plotting function

    input_dir -- directory, where the function expects the data to be found in a serie of *.tsv files for each population pair,
                 the names of these files should start with 'POP1_POP2' string

    begin, end -- limit the X-axis (position), displayed area

    populations -- by default, the '1000Genomes' option, 26 populations from 1000 Genomes Project are expected. If your population/classes set is different,
                   provide an iterable with population/classes names. For example, with populations=['pop1', 'pop2', 'pop3'] this function will search the input directory (input_dir)
                   for *.tsv files with all combinations of pairwise values, named pop1_pop2.*.tsv, pop1_pop3.*.tsv, pop2_pop1.*.tsv etc.

                   The output 'big_df' DataFrame with following 6 rows of pairwise values will be created:
                   pop1 vs pop2
                   pop1 vs pop3
                   pop2 vs pop1
                   pop2 vs pop3
                   pop3 vs pop1
                   pop3 vs pop2

    """
    df_list = []
    pop_id_list = []
    different_dfs = False

    segment_files = glob.glob(os.path.join(input_dir, "*.tsv"))

    # reading the input files, saving only the regions between BEGIN and END to process further
    index = 1
    for segment_file in segment_files:
        # segment_files is something like ACB_KHV.tsv or ACB_KHV.some_more_info.tsv
        pop_pair = os.path.splitext(os.path.basename(segment_file))[0].split(".")[0]
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
        # select the descending ranks that are significant for pop1_pop2 (pop1 is under selection)
        left_df = df[["variant_pos", "-log10_p_value_descending"]].copy()
        left_df.rename(columns={"-log10_p_value_descending": pop_pair}, inplace=True)
        left_df = left_df.set_index("variant_pos").T
        transp_list.append(left_df)

        # select the ascending ranks that are significant for pop2_pop1 (pop2 is under selection)
        reverse_pop_pair = "_".join(
            pop_pair.split("_")[::-1]
        )  # change name pop1_pop2 to pop2_pop1

        right_df = df[["variant_pos", "-log10_p_value_ascending"]].copy()
        right_df.rename(
            columns={"-log10_p_value_ascending": reverse_pop_pair}, inplace=True
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
    big_df.first_pop.cat.set_categories(population_sorter, inplace=True)

    # set pop2 to be a categorical column with value order defined by sorter
    big_df.second_pop = big_df.second_pop.astype("category")
    big_df.second_pop.cat.set_categories(population_sorter, inplace=True)

    # sort df by pop1 and withing pop1 by pop2
    big_df.sort_values(["first_pop", "second_pop"], inplace=True)

    # drop the temporary columns
    big_df.drop(["first_pop", "second_pop"], axis=1, inplace=True)

    # label it just by the pop1 (which is gonna be printed with plot ticks)
    pop_labels = big_df.index.values
    pop_labels = [pop.split("_")[0] for pop in pop_labels]
    big_df.index = pop_labels

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
    cbar_vmin=1,
    cbar_vmax=4.853,
    ylabel=False,
    xlabel=False,
    cbar_ticks=False,
):
    """
    Read input DataFrame and create the ExP heatmap accordingly.

    input_df -- input pandas DataFrame, data to display

    begin, end -- limit the X-axis (position), displayed area

    title -- title of the graph

    cmap -- seaborn colormap, 'Blues' or 'crest' work well

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

    #########################
    # create the ExP figure #
    #########################
    print("Creating heatmap")

    fig, ax = plt.subplots(figsize=(15, 5))

    if not cmap:
        cmap = "Blues"

    # draw default exp heatmap with 26 populations from 1000 Genomes Project
    if populations == "1000Genomes":
        sns.heatmap(
            input_df,
            yticklabels="auto",
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
            "population pairings\n\nAMR  |    EUR     |     EAS    |    SAS     |       AFR  "
        )
        ax.set_xlabel("{:,} - {:,}".format(begin, end))

    # draw custom exp heatmap with user-defined populations (number of pops, labels)
    #
    else:
        sns.heatmap(
            input_df,
            yticklabels="auto",
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

    # optionally add vertical line in the middle of the figure
    if vertical_line:

        middle = int(input_df.shape[1] / 2)
        ax.axvline(x=middle, linewidth=1, color="grey")

    print("Savig heatmap")

    ax.figure.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"Saved into {output}.png")


def plot(xpehh_dir, begin, end, title, cmap, output):
    data_to_plot = create_plot_input(xpehh_dir, begin=begin, end=end)
    plot_exp_heatmap(
        data_to_plot, begin=begin, end=end, title=title, cmap=cmap, output=output
    )
