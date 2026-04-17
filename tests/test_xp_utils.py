import numpy as np

from exp_heatmap.xp_utils import filter_by_AF


def _build_fake_callset(include_af):
    callset = {
        "calldata/GT": np.array(
            [
                [[0, 0], [0, 1]],
                [[0, 1], [1, 1]],
                [[0, 0], [0, 0]],
            ],
            dtype="i1",
        ),
        "variants/POS": np.array([10, 20, 30]),
    }
    if include_af:
        callset["variants/AF"] = np.array([[0.25], [0.75], [0.0]])
    return callset


def test_filter_by_af_uses_precomputed_values_when_available():
    _, positions = filter_by_AF(_build_fake_callset(include_af=True), af_threshold=0.2, chunked=False)
    assert list(positions) == [10, 20]


def test_filter_by_af_falls_back_to_genotype_counts_when_af_missing():
    _, positions = filter_by_AF(_build_fake_callset(include_af=False), af_threshold=0.2, chunked=False)
    assert list(positions) == [10, 20]
