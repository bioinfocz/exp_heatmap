import numpy as np


def create_pop_pairs(panel):
    populations = np.unique(panel["pop"].values)

    pop_pairs = []
    for i in range(len(populations)):
        pop1 = populations[i]
        for j in range(i + 1, len(populations)):
            pop2 = populations[j]
            pop_pairs.append((pop1, pop2))

    return pop_pairs



def get_haplotypes(gt_array, panel, pop):
    # get the indices of samples which belong to given population
    indices_pop = panel.index[panel["pop"] == pop]

    # get genotype data belonging only to given population
    gt_pop = gt_array.take(indices_pop, axis=1)

    return gt_pop.to_haplotypes()



def get_pop_allele_counts(gt, panel, pop):
    """
    Returns allele counts for given population
    """
    
    # get the indices of samples (individuals) which belong to pop
    indices_pop = panel.index[panel["pop"] == pop]

    # get genotype data belonging only to pop
    gt_pop = gt.take(indices_pop, axis=1)

    # get the allel counts for population (input for pbs)
    ac = gt_pop.count_alleles()
    return ac
