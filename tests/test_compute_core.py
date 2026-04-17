import numpy as np
import pytest
import pandas as pd

from exp_heatmap.compute_core import _build_output_dataframe


def test_build_output_dataframe_preserves_pair_specific_nan_values():
    base_df = pd.DataFrame(
        {
            "name": ["chr1", "chr1", "chr1"],
            "variant_pos": [10, 20, 30],
        }
    )

    result_df = _build_output_dataframe(base_df, "xpehh", [0.2, np.nan, 1.5])

    assert result_df["xpehh"].isna().sum() == 1
    assert np.isnan(result_df.loc[1, "-log10_p_value_ascending"])
    assert np.isnan(result_df.loc[1, "-log10_p_value_descending"])
    assert result_df.loc[2, "-log10_p_value_descending"] > result_df.loc[0, "-log10_p_value_descending"]
