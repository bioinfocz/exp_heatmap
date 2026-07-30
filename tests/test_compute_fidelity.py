"""
Fidelity of the compute wrapper against the underlying scikit-allel implementations.

Reviewer feedback asked for a comparison of the computed statistics against established
tools. ExP Heatmap does not implement any statistic itself: it selects haplotypes or allele
counts per population pair and delegates to scikit-allel. Comparing against a separate
implementation such as selscan would therefore compare scikit-allel with selscan rather
than test this package. What is testable here, and what these checks assert, is that the
wrapper hands the right data to the library and stores the returned values unaltered.
"""

import numpy as np
import pandas as pd
import pytest

allel = pytest.importorskip("allel")

from exp_heatmap import compute_core, xp_utils

# Relative tolerance absorbing TSV serialization loss only, roughly 10,000x the ULP.
TOLERANCE = 1e-12


def _phased_callset(n_variants=60, n_samples=12, seed=0):
    """A small phased biallelic callset with an INFO/AF field, as a dict-backed callset."""
    rng = np.random.default_rng(seed)
    genotypes = rng.integers(0, 2, size=(n_variants, n_samples, 2)).astype("i1")
    # Keep every variant clear of the 0.05 alternate-allele-frequency filter.
    genotypes[:, 0, 0] = 1
    alt_frequency = genotypes.reshape(n_variants, -1).mean(axis=1)
    return {
        "calldata/GT": genotypes,
        "variants/POS": np.arange(1, n_variants + 1, dtype="i4") * 100,
        "variants/AF": alt_frequency.reshape(-1, 1),
        "samples": np.array([f"S{i:02d}" for i in range(n_samples)]),
    }


def _panel(n_samples=12, n_populations=3):
    populations = [f"P{i}" for i in range(n_populations)]
    return pd.DataFrame(
        {
            "sample": [f"S{i:02d}" for i in range(n_samples)],
            "pop": [populations[i % n_populations] for i in range(n_samples)],
            "super_pop": ["SUP"] * n_samples,
        }
    )


def test_compute_xpehh_matches_a_direct_scikit_allel_call(tmp_path, monkeypatch):
    """
    The xpehh values written per pair must equal a direct allel.xpehh call on the same
    haplotypes, with the same parameters, element for element including NaNs.
    """
    callset = _phased_callset()
    panel = _panel()
    panel_path = tmp_path / "panel.tsv"
    panel.to_csv(panel_path, sep="\t", index=False)

    monkeypatch.setattr(compute_core.zarr, "open_group", lambda *a, **k: callset)

    output_dir = tmp_path / "out"
    compute_core.run(
        zarr_dir=str(tmp_path / "unused.zarr"),
        panel_file=str(panel_path),
        output_dir=str(output_dir),
        test="xpehh",
    )

    gt, positions = xp_utils.filter_by_AF(callset, 0.05, False)
    pairs = xp_utils.create_pop_pairs(panel)
    assert pairs, "expected at least one population pair"

    for pop1, pop2 in pairs:
        expected = allel.xpehh(
            h1=xp_utils.get_haplotypes(gt, panel, pop1),
            h2=xp_utils.get_haplotypes(gt, panel, pop2),
            pos=positions,
            map_pos=None,
            min_ehh=0.05,
            include_edges=False,
            gap_scale=20000,
            max_gap=200000,
            is_accessible=None,
            use_threads=True,
        )
        written = pd.read_csv(output_dir / f"{pop1}_{pop2}.tsv", sep="\t")["xpehh"].to_numpy()

        assert written.shape == expected.shape
        np.testing.assert_array_equal(np.isnan(written), np.isnan(expected))
        finite = ~np.isnan(expected)
        # The values pass through a float -> text -> float round trip when the TSV is
        # written, which costs up to an ULP. TOLERANCE is far below any difference a real
        # wrapper fault would cause: passing the wrong haplotypes, parameters or pair
        # ordering changes these values in the first decimal places, not the twelfth.
        np.testing.assert_allclose(
            written[finite], expected[finite], rtol=TOLERANCE, atol=0,
            err_msg=f"{pop1}_{pop2}",
        )


def test_compute_hudson_fst_matches_a_direct_scikit_allel_call(tmp_path, monkeypatch):
    """Same fidelity check for an unphased, allele-count-based statistic."""
    callset = _phased_callset(seed=7)
    panel = _panel()
    panel_path = tmp_path / "panel.tsv"
    panel.to_csv(panel_path, sep="\t", index=False)

    monkeypatch.setattr(compute_core.zarr, "open_group", lambda *a, **k: callset)

    output_dir = tmp_path / "out"
    compute_core.run(
        zarr_dir=str(tmp_path / "unused.zarr"),
        panel_file=str(panel_path),
        output_dir=str(output_dir),
        test="hudson_fst",
    )

    gt, _ = xp_utils.filter_by_AF(callset, 0.05, False)
    for pop1, pop2 in xp_utils.create_pop_pairs(panel):
        num, den = allel.hudson_fst(
            ac1=xp_utils.get_pop_allele_counts(gt, panel, pop1),
            ac2=xp_utils.get_pop_allele_counts(gt, panel, pop2),
        )
        expected = num / den
        written = pd.read_csv(
            output_dir / f"{pop1}_{pop2}.tsv", sep="\t"
        )["hudson_fst"].to_numpy()

        np.testing.assert_array_equal(np.isnan(written), np.isnan(expected))
        finite = ~np.isnan(expected)
        np.testing.assert_allclose(written[finite], expected[finite], rtol=TOLERANCE, atol=0)


def test_compute_is_deterministic_across_runs(tmp_path, monkeypatch):
    """Two runs on identical input must produce byte-identical per-pair outputs."""
    callset = _phased_callset(seed=3)
    panel = _panel()
    panel_path = tmp_path / "panel.tsv"
    panel.to_csv(panel_path, sep="\t", index=False)

    monkeypatch.setattr(compute_core.zarr, "open_group", lambda *a, **k: callset)

    first, second = tmp_path / "run1", tmp_path / "run2"
    for output_dir in (first, second):
        compute_core.run(
            zarr_dir=str(tmp_path / "unused.zarr"),
            panel_file=str(panel_path),
            output_dir=str(output_dir),
            test="xpehh",
        )

    names = sorted(path.name for path in first.glob("*.tsv"))
    assert names, "expected per-pair outputs"
    for name in names:
        assert (first / name).read_bytes() == (second / name).read_bytes(), name
