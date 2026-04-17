from exp_heatmap.benchmark_utils import (
    build_comparison_regions,
    centered_window,
    normalize_population_counts,
    resolve_population_selection,
)


def test_centered_window_clips_to_available_range():
    assert centered_window(100, 500, 1000) == (100, 500)
    assert centered_window(100, 500, 200) == (200, 400)


def test_normalize_population_counts_defaults_to_small_and_full_runs():
    assert normalize_population_counts(("A", "B")) == [2]
    assert normalize_population_counts(("A", "B", "C", "D", "E")) == [3, 5]


def test_resolve_population_selection_preserves_full_1000g_mode():
    populations = ("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "BEB", "GIH", "ITU", "PJL", "STU", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "CLM", "MXL", "PEL", "PUR")
    assert resolve_population_selection(populations, "1000Genomes", 26) == "1000Genomes"
    assert resolve_population_selection(populations, "1000Genomes", 3) == ("ACB", "ASW", "ESN")


def test_build_comparison_regions_returns_two_nonempty_subregions():
    region1, region2 = build_comparison_regions(100, 1000)
    assert region1[0] < region1[1]
    assert region2[0] < region2[1]
    assert region1[1] <= region2[0]
