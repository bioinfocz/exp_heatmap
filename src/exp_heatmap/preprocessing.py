import allel


def filter_by_AF(callset, af_t):
    """
    returns genotype array with variants
    with alternate alllel frequency > af_t

    and array of positions without filtered variants
    """

    # acess alternate allele frequencies
    af = callset["variants/AF"][:]

    loc_variant_selection = af[:, 0] > af_t

    # acces the genotype data from zarr
    gt_zarr = callset["calldata/GT"]

    # if big, load the genotype as chunked array
    # gt = allel.GenotypeChunkedArray(gt_zarr)
    gt = allel.GenotypeArray(gt_zarr)

    # GET ONLY VARIANTS BELLOW 0.05
    gt_variant_selection = gt.compress(loc_variant_selection, axis=0)

    position_selection = allel.SortedIndex(callset["variants/POS"])
    position_selection = position_selection.compress(loc_variant_selection, axis=0)

    return gt_variant_selection, position_selection
