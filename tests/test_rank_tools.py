import numpy as np
import pytest

from exp_heatmap.rank_tools import compute_rank_scores_for_series


def test_compute_rank_scores_for_series_preserves_nan_and_competition_ranks():
    values = [5.0, np.nan, 5.0, 1.0]

    descending = compute_rank_scores_for_series(values, ascending=False)
    ascending = compute_rank_scores_for_series(values, ascending=True)

    assert np.isnan(descending.iloc[1])
    assert np.isnan(ascending.iloc[1])

    assert descending.iloc[0] == pytest.approx(0.477, abs=1e-3)
    assert descending.iloc[2] == pytest.approx(0.477, abs=1e-3)
    assert descending.iloc[3] == pytest.approx(0.0, abs=1e-9)

    assert ascending.iloc[3] == pytest.approx(0.477, abs=1e-3)
    assert ascending.iloc[0] == pytest.approx(0.176, abs=1e-3)
    assert ascending.iloc[2] == pytest.approx(0.176, abs=1e-3)


def test_compute_rank_scores_for_series_returns_all_nan_for_empty_input():
    scores = compute_rank_scores_for_series([np.nan, np.nan], ascending=True)
    assert scores.isna().all()
