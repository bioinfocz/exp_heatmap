import pandas as pd
import pytest

import numpy as np

from exp_heatmap.plot import (
    CBAR_VMAX_1000G,
    create_plot_input,
    downsample_heatmap_columns,
    normalize_rank_score_mode,
    plot,
    plot_exp_heatmap,
    populations_1000genomes,
    resolve_population_configuration,
    warn_if_values_are_clipped,
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


def test_clipping_warning_reports_cells_above_the_color_scale_maximum(caplog):
    frame = pd.DataFrame([[1.5, 6.5], [4.0, 5.2]], index=["A_B", "B_A"], columns=[10, 20])

    with caplog.at_level("WARNING", logger="exp_heatmap.plot"):
        warn_if_values_are_clipped(frame, CBAR_VMAX_1000G)

    assert "2 of 4 rendered cells exceed" in caplog.text
    assert "6.500" in caplog.text


def test_clipping_warning_is_silent_when_data_fits_the_scale(caplog):
    frame = pd.DataFrame([[1.5, 2.0], [3.0, 4.8]], index=["A_B", "B_A"], columns=[10, 20])

    with caplog.at_level("WARNING", logger="exp_heatmap.plot"):
        warn_if_values_are_clipped(frame, CBAR_VMAX_1000G)

    assert caplog.text == ""


def test_clipping_warning_is_silent_when_the_scale_is_autoscaled(caplog):
    # Custom panels pass vmax=None, so nothing can be clipped.
    frame = pd.DataFrame([[1.5, 99.0]], index=["A_B"], columns=[10, 20])

    with caplog.at_level("WARNING", logger="exp_heatmap.plot"):
        warn_if_values_are_clipped(frame, None)

    assert caplog.text == ""


def test_clipping_warning_ignores_nan_only_input(caplog):
    frame = pd.DataFrame([[np.nan, np.nan]], index=["A_B"], columns=[10, 20])

    with caplog.at_level("WARNING", logger="exp_heatmap.plot"):
        warn_if_values_are_clipped(frame, CBAR_VMAX_1000G)

    assert caplog.text == ""


def test_clipping_warning_does_not_count_nan_cells(caplog):
    frame = pd.DataFrame([[np.nan, 6.5]], index=["A_B"], columns=[10, 20])

    with caplog.at_level("WARNING", logger="exp_heatmap.plot"):
        warn_if_values_are_clipped(frame, CBAR_VMAX_1000G)

    assert "1 of 2 rendered cells exceed" in caplog.text


def _write_pair_files(directory, pairs):
    for pop1, pop2 in pairs:
        pd.DataFrame(
            [
                {
                    "name": "chr21",
                    "variant_pos": 100,
                    "-log10_p_value_ascending": 1.0,
                    "-log10_p_value_descending": 2.0,
                }
            ]
        ).to_csv(directory / f"{pop1}_{pop2}.tsv", sep="\t", index=False)


def test_resolve_population_configuration_maps_declared_canonical_panel_to_1000g():
    assert resolve_population_configuration(populations_1000genomes, []) == "1000Genomes"


def test_resolve_population_configuration_respects_a_declared_custom_order():
    reordered = tuple(reversed(populations_1000genomes))
    assert resolve_population_configuration(reordered, []) == reordered
    assert resolve_population_configuration(("GWD", "MSL", "ESN"), []) == ("GWD", "MSL", "ESN")


def test_declared_populations_reject_a_compute_directory_missing_a_population(tmp_path):
    # A directory that is complete for three populations but was meant to hold four.
    _write_pair_files(tmp_path, [("ESN", "GWD"), ("ESN", "MSL"), ("GWD", "MSL")])
    declared = ("ESN", "FULA", "GWD", "MSL")

    with pytest.raises(ValueError) as excinfo:
        create_plot_input(str(tmp_path), start=100, end=100, populations=declared)

    message = str(excinfo.value)
    assert "3 of 6 population pairs are missing" in message
    assert "no output at all for FULA" in message


def test_declared_populations_report_a_single_missing_pair(tmp_path):
    # Both populations are present elsewhere, so only the pair itself is missing.
    _write_pair_files(tmp_path, [("ESN", "MSL"), ("GWD", "MSL")])

    with pytest.raises(ValueError) as excinfo:
        create_plot_input(
            str(tmp_path), start=100, end=100, populations=("ESN", "GWD", "MSL")
        )

    message = str(excinfo.value)
    assert "Missing pairs: ESN_GWD" in message
    assert "no output at all" not in message


def test_declared_complete_panel_passes_the_completeness_check(tmp_path):
    _write_pair_files(tmp_path, [("ESN", "GWD"), ("ESN", "MSL"), ("GWD", "MSL")])

    plot_input = create_plot_input(
        str(tmp_path), start=100, end=100, populations=("ESN", "GWD", "MSL")
    )

    assert plot_input.shape[0] == 6
    assert list(plot_input.index)[:2] == ["ESN_GWD", "ESN_MSL"]


def test_inferred_populations_silently_accept_the_same_incomplete_directory(tmp_path):
    # Without a declared panel the same directory is a self-consistent 3-population run,
    # which is why --populations exists.
    _write_pair_files(tmp_path, [("ESN", "GWD"), ("ESN", "MSL"), ("GWD", "MSL")])

    plot_input = create_plot_input(str(tmp_path), start=100, end=100)

    assert plot_input.attrs["population_mode"] == ("ESN", "GWD", "MSL")
    assert plot_input.shape[0] == 6


def test_create_plot_input_returns_a_numeric_column_index(tmp_path):
    # The two temporary sort-key columns used for row ordering promote the position index
    # to object dtype; create_plot_input must restore it so regions can be sliced.
    _write_pair_files(tmp_path, [("ESN", "GWD"), ("ESN", "MSL"), ("GWD", "MSL")])

    plot_input = create_plot_input(str(tmp_path), start=100, end=100)

    assert pd.api.types.is_numeric_dtype(plot_input.columns)


def test_plot_exp_heatmap_crops_to_a_region_that_is_not_an_exact_column(tmp_path):
    frame = pd.DataFrame(
        [[1.5, 2.0, 2.5], [1.6, 2.1, 2.6]],
        index=["A_B", "B_A"],
        columns=[100, 200, 300],
    )

    ax = plot_exp_heatmap(
        frame, start=150, end=250, output=str(tmp_path / "cropped"), populations=("A", "B")
    )

    # 200 is the only column inside 150-250.
    assert ax.collections[0].get_array().shape[0] == 2


def test_plot_exp_heatmap_falls_back_when_the_column_index_is_not_numeric(tmp_path):
    # Frames built by hand rather than by create_plot_input may carry an object index;
    # slicing one raises TypeError, which must fall back rather than propagate.
    frame = pd.DataFrame(
        [[1.5, 2.0, 2.5], [1.6, 2.1, 2.6]],
        index=["A_B", "B_A"],
        columns=pd.Index([100, 200, 300], dtype=object),
    )

    ax = plot_exp_heatmap(
        frame, start=150, end=250, output=str(tmp_path / "fallback"), populations=("A", "B")
    )

    assert ax.collections[0].get_array().shape[0] == 2


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
