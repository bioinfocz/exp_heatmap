from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from exp_heatmap import cli as cli_module
from exp_heatmap import interactive as interactive_module
from exp_heatmap import plot as plot_module


class _DummyLogger:
    def info(self, *args, **kwargs):
        pass


def _patch_logging(monkeypatch):
    monkeypatch.setattr(cli_module, "setup_logging", lambda *args, **kwargs: None)
    monkeypatch.setattr(cli_module, "_finalize_log_file", lambda logger: None)
    monkeypatch.setattr(cli_module, "get_logger", lambda name: _DummyLogger())


def test_focus_command_passes_rank_scores_and_max_columns(monkeypatch, tmp_path):
    _patch_logging(monkeypatch)
    captured = {}

    monkeypatch.setattr(
        plot_module,
        "create_plot_input",
        lambda input_dir, start, end, rank_scores: pd.DataFrame([[1.0, 2.0]], index=["CEU_YRI"], columns=[start, end]),
    )

    def fake_focus(input_df, focus_population, start, end, title, output, colorscale, max_columns):
        captured["focus_population"] = focus_population
        captured["start"] = start
        captured["end"] = end
        captured["colorscale"] = colorscale
        captured["max_columns"] = max_columns

    monkeypatch.setattr(interactive_module, "create_population_focus_view", fake_focus)

    result = CliRunner().invoke(
        cli_module.cli,
        [
            "focus",
            str(tmp_path),
            "--start",
            "10",
            "--end",
            "20",
            "--population",
            "CEU",
            "--rank-scores",
            "ascending",
            "--max-columns",
            "40",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["focus_population"] == "CEU"
    assert captured["start"] == 10
    assert captured["end"] == 20
    assert captured["colorscale"] == "Blues"
    assert captured["max_columns"] == 40


def test_compare_command_loads_union_region_and_passes_max_columns(monkeypatch, tmp_path):
    _patch_logging(monkeypatch)
    captured = {}

    def fake_create_plot_input(input_dir, start, end, rank_scores):
        captured["load_start"] = start
        captured["load_end"] = end
        captured["rank_scores"] = rank_scores
        return pd.DataFrame([[1.0, 2.0, 3.0]], index=["CEU_YRI"], columns=[100, 200, 300])

    def fake_compare(input_df, region1, region2, title, output, colorscale, max_columns):
        captured["region1"] = region1
        captured["region2"] = region2
        captured["max_columns"] = max_columns

    monkeypatch.setattr(plot_module, "create_plot_input", fake_create_plot_input)
    monkeypatch.setattr(interactive_module, "create_comparison_view", fake_compare)

    result = CliRunner().invoke(
        cli_module.cli,
        [
            "compare",
            str(tmp_path),
            "--start-1",
            "120",
            "--end-1",
            "180",
            "--start-2",
            "220",
            "--end-2",
            "260",
            "--rank-scores",
            "descending",
            "--max-columns",
            "25",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["load_start"] == 120
    assert captured["load_end"] == 260
    assert captured["rank_scores"] == "descending"
    assert captured["region1"] == (120, 180)
    assert captured["region2"] == (220, 260)
    assert captured["max_columns"] == 25


def test_summary_command_writes_tsv_and_plots(monkeypatch, tmp_path):
    _patch_logging(monkeypatch)
    runner = CliRunner()
    summary_output = tmp_path / "summary"
    captured = {}

    monkeypatch.setattr(
        plot_module,
        "create_plot_input",
        lambda input_dir, start, end, rank_scores: pd.DataFrame([[1.0, 2.0]], index=["CEU_YRI"], columns=[start, end]),
    )

    summary_df = pd.DataFrame([[1.0, 2.0]], index=["EUR_AFR"], columns=[10, 20])
    monkeypatch.setattr(plot_module, "summarize_by_superpopulation", lambda *args, **kwargs: summary_df)

    def fake_plot(input_df, start, end, title, output, colorscale, max_columns, show_superpop_annotations):
        captured["start"] = start
        captured["end"] = end
        captured["colorscale"] = colorscale
        captured["max_columns"] = max_columns

    monkeypatch.setattr(interactive_module, "plot_interactive_heatmap", fake_plot)

    result = runner.invoke(
        cli_module.cli,
        [
            "summary",
            str(tmp_path),
            "--start",
            "10",
            "--end",
            "20",
            "--out",
            str(summary_output),
            "--max-columns",
            "15",
        ],
    )

    assert result.exit_code == 0, result.output
    assert Path(f"{summary_output}.tsv").exists()
    assert captured["start"] == 10
    assert captured["end"] == 20
    assert captured["max_columns"] == 15


def test_regions_command_writes_output(monkeypatch, tmp_path):
    _patch_logging(monkeypatch)
    runner = CliRunner()
    output_path = tmp_path / "regions.tsv"

    monkeypatch.setattr(
        plot_module,
        "create_plot_input",
        lambda input_dir, start, end, rank_scores: pd.DataFrame([[1.0, 2.0]], index=["CEU_YRI"], columns=[start, end]),
    )
    monkeypatch.setattr(
        plot_module,
        "extract_top_regions",
        lambda *args, **kwargs: pd.DataFrame([{"start": 10, "end": 20, "center": 15, "mean_score": 1.2, "max_score": 2.0, "top_population_pair": "CEU_YRI", "n_variants": 2}]),
    )

    result = runner.invoke(
        cli_module.cli,
        [
            "regions",
            str(tmp_path),
            "--start",
            "10",
            "--end",
            "20",
            "--out",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_path.exists()


def test_filter_vcf_command_writes_filtered_output(monkeypatch, tmp_path):
    _patch_logging(monkeypatch)
    runner = CliRunner()
    input_path = tmp_path / "input.vcf"
    output_path = tmp_path / "output.vcf"
    input_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\n"
        "21\t110\trs2\tAT\tG\t.\tPASS\t.\tGT\t0|1\n"
    )

    result = runner.invoke(
        cli_module.cli,
        [
            "filter-vcf",
            str(input_path),
            "--out",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_path.read_text().splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1",
    ]


def test_filter_vcf_command_supports_region_selection(monkeypatch, tmp_path):
    _patch_logging(monkeypatch)
    runner = CliRunner()
    input_path = tmp_path / "input.vcf"
    output_path = tmp_path / "output.vcf"
    input_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr2\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\n"
        "chr2\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1\n"
        "chr2\t300\trs3\tG\tA\t.\tPASS\t.\tGT\t0|1\n"
    )

    result = runner.invoke(
        cli_module.cli,
        [
            "filter-vcf",
            str(input_path),
            "--out",
            str(output_path),
            "--chrom",
            "2",
            "--start",
            "150",
            "--end",
            "250",
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_path.read_text().splitlines() == [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "chr2\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1",
    ]
