import glob
import os

import pandas as pd

from exp_heatmap.plot import populations_1000genomes, resolve_population_configuration


def get_segment_files(input_dir):
    return sorted(glob.glob(os.path.join(input_dir, "*.tsv")))


def read_variant_range(input_dir):
    segment_files = get_segment_files(input_dir)
    if not segment_files:
        raise ValueError(f"No TSV files found in '{input_dir}'")

    first_df = pd.read_csv(segment_files[0], sep="\t", usecols=["variant_pos"])
    return int(first_df["variant_pos"].min()), int(first_df["variant_pos"].max())


def infer_available_populations(input_dir):
    segment_files = get_segment_files(input_dir)
    population_mode = resolve_population_configuration("1000Genomes", segment_files)
    if population_mode == "1000Genomes":
        return tuple(populations_1000genomes), "1000Genomes"
    return tuple(population_mode), tuple(population_mode)


def centered_window(range_start, range_end, width):
    width = int(width)
    if width <= 0:
        raise ValueError("Window width must be positive")

    available_width = range_end - range_start
    if width >= available_width:
        return int(range_start), int(range_end)

    center = (range_start + range_end) // 2
    start = center - width // 2
    end = start + width

    if start < range_start:
        start = range_start
        end = min(range_end, start + width)
    if end > range_end:
        end = range_end
        start = max(range_start, end - width)

    return int(start), int(end)


def normalize_population_counts(available_populations, requested_counts=None):
    max_count = len(available_populations)
    if max_count < 2:
        raise ValueError("At least two populations are required for benchmarking")

    if requested_counts:
        counts = sorted(
            {
                int(count)
                for count in requested_counts
                if 2 <= int(count) <= max_count
            }
        )
        if not counts:
            raise ValueError("No valid population counts remained after filtering")
        return counts

    if max_count <= 3:
        return [max_count]
    return sorted({min(3, max_count), max_count})


def resolve_population_selection(available_populations, full_population_mode, pop_count):
    selection = tuple(available_populations[:pop_count])
    if (
        full_population_mode == "1000Genomes"
        and pop_count == len(populations_1000genomes)
        and selection == tuple(populations_1000genomes)
    ):
        return "1000Genomes"
    return selection


def build_comparison_regions(start, end):
    span = end - start
    if span <= 2:
        raise ValueError("Region is too narrow to construct a comparison view")

    sub_width = max(1, span // 3)
    region1 = (start, min(end, start + sub_width))
    region2 = (max(start, end - sub_width), end)
    if region1[0] >= region1[1] or region2[0] >= region2[1]:
        raise ValueError("Failed to construct non-empty comparison regions")
    return region1, region2
