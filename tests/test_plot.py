import pandas as pd
import pytest

from exp_heatmap.plot import (
    create_plot_input,
    downsample_heatmap_columns,
    normalize_rank_score_mode,
    plot,
)


def test_downsample_heatmap_columns_uses_midpoint_positions_and_max_reducer():
    input_df = pd.DataFrame(
        [[1, 3, 2, 4], [0, 5, 1, 2]],
        index=["A_B", "B_A"],
        columns=[100, 200, 300, 400],
    )

    downsampled, changed = downsample_heatmap_columns(input_df, max_columns=2, aggregation="max")

    assert changed is True
    assert list(downsampled.columns) == [200, 400]
    assert downsampled.loc["A_B", 200] == 3
    assert downsampled.loc["B_A", 200] == 5
    assert downsampled.loc["A_B", 400] == 4


def test_downsample_heatmap_columns_supports_mean_reducer():
    input_df = pd.DataFrame(
        [[1.0, 3.0, 2.0, 4.0]],
        index=["A_B"],
        columns=[100, 200, 300, 400],
    )

    downsampled, changed = downsample_heatmap_columns(input_df, max_columns=2, aggregation="mean")

    assert changed is True
    assert downsampled.loc["A_B", 200] == pytest.approx(2.0)
    assert downsampled.loc["A_B", 400] == pytest.approx(3.0)


def test_normalize_rank_score_mode_accepts_legacy_alias():
    assert normalize_rank_score_mode("2-tailed") == "directional"
    assert normalize_rank_score_mode("directional") == "directional"


def test_create_plot_input_infers_custom_populations_from_files(tmp_path):
    records = [
        ("FULA_GWD.tsv", "chr21", 100, 1.0, 2.0),
        ("FULA_JOLA.tsv", "chr21", 100, 1.5, 2.5),
        ("GWD_JOLA.tsv", "chr21", 100, 2.0, 3.0),
    ]
    for filename, name, pos, asc, desc in records:
        pd.DataFrame(
            [
                {
                    "name": name,
                    "variant_pos": pos,
                    "-log10_p_value_ascending": asc,
                    "-log10_p_value_descending": desc,
                }
            ]
        ).to_csv(tmp_path / filename, sep="\t", index=False)

    plot_input = create_plot_input(str(tmp_path), start=100, end=100)

    assert plot_input.attrs["population_mode"] == ("FULA", "GWD", "JOLA")
    assert list(plot_input.index) == [
        "FULA_GWD",
        "FULA_JOLA",
        "GWD_FULA",
        "GWD_JOLA",
        "JOLA_FULA",
        "JOLA_GWD",
    ]


def test_plot_passes_inferred_custom_populations(monkeypatch, tmp_path):
    for filename in ["FULA_GWD.tsv", "FULA_JOLA.tsv", "GWD_JOLA.tsv"]:
        pd.DataFrame(
            [
                {
                    "name": "chr21",
                    "variant_pos": 100,
                    "-log10_p_value_ascending": 1.0,
                    "-log10_p_value_descending": 2.0,
                }
            ]
        ).to_csv(tmp_path / filename, sep="\t", index=False)

    captured = {}

    def fake_plot_exp_heatmap(input_df, start, end, title, cmap, output, populations, max_columns, column_aggregation, dpi):
        captured["populations"] = populations
        captured["start"] = start
        captured["end"] = end
        return "ok"

    monkeypatch.setattr("exp_heatmap.plot.plot_exp_heatmap", fake_plot_exp_heatmap)

    result = plot(str(tmp_path), start=100, end=100, title="GGVP")

    assert result == "ok"
    assert captured["populations"] == ("FULA", "GWD", "JOLA")
    assert captured["start"] == 100
    assert captured["end"] == 100
