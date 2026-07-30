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


def test_filter_by_af_reports_kept_and_removed_counts(caplog):
    with caplog.at_level("INFO", logger="exp_heatmap.xp_utils"):
        filter_by_AF(_build_fake_callset(include_af=True), af_threshold=0.2, chunked=False)

    assert "kept 2 of 3 variants (1 removed)" in caplog.text
    assert "total alternate-allele frequency > 0.2" in caplog.text


def test_filter_by_af_reports_the_frequency_source(caplog):
    with caplog.at_level("INFO", logger="exp_heatmap.xp_utils"):
        filter_by_AF(_build_fake_callset(include_af=True), af_threshold=0.2, chunked=False)

    assert "precomputed alternate allele frequencies from variants/AF" in caplog.text


def test_filter_by_af_warns_when_falling_back_to_genotypes(caplog):
    with caplog.at_level("WARNING", logger="exp_heatmap.xp_utils"):
        filter_by_AF(_build_fake_callset(include_af=False), af_threshold=0.2, chunked=False)

    assert "variants/AF is missing" in caplog.text
