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


def plot(xpehh_dir, begin, end, title, cmap):
    df_list = []
    pop_id_list = []
    different_dfs = False

    segment_files = glob.glob(os.path.join(xpehh_dir, "*.tsv"))

    for segment_file in segment_files:
        # segment_files is something like ACB_KHV.tsv
        pop_pair = os.path.splitext(os.path.basename(segment_file))[0]
        pop_id_list.append(pop_pair)

        print("Loading {} from {}".format(pop_pair, segment_file))

        segments = pd.read_csv(segment_file, sep="\t")
        segments = segments[
            (segments.variant_pos >= begin) & (segments.variant_pos <= end)
        ]

        df_list.append(segments)

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

    print("Creating heatmap")

    fig, ax = plt.subplots(figsize=(15, 5))

    if not cmap:
        cmap = "Blues"

    sns.heatmap(
        big_df,
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

    ax.set_title(title)
    ax.set_ylabel(
        "population pairings\n\nAMR  |    EUR     |     EAS    |    SAS     |       AFR  "
    )
    ax.set_xlabel("-500 kbp" + (" " * 176) + "+500kbp \nsnp position relative to ")
    middle = int(big_df.shape[1] / 2)
    ax.axvline(x=middle, linewidth=2, color="grey")

    print("Savig heatmap")

    ax.figure.savefig(title, dpi=400, bbox_inches="tight")
    plt.close(fig)
