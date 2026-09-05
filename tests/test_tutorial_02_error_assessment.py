"""Regression test for tutorial/02_chemgen_error_assessment.

Reproduces `chemgen.py FFCM2_model . --custom-test custom_test.py --n-points-test N --compile`,
runs the resulting binary, and turns the README's prose tolerance ("mean
relative error < 0.01%, only a few outliers above 10%") into real assertions
against l2_norm_results.csv.
"""
import pytest

from conftest import run_binary, run_chemgen
from utils import read_l2_norm_csv

MEAN_L2_TOLERANCE = 1e-3  # README: "mean around a relative error of less than 0.01 pct"
OUTLIER_L2_TOLERANCE = 0.10  # README: "only a few outliers ever making it above a 10% difference"
MAX_OUTLIER_FRACTION = 0.05  # "only a few" -- cap it at 5% of sampled states


@pytest.mark.parametrize("tutorial_dir", ["02_chemgen_error_assessment"], indirect=True)
def test_ffcm2_source_term_error_distribution(tutorial_dir, cxx_compiler):
    n_points = 200  # reduced from the tutorial's documented 1000 for CI speed; see test_full below
    generation = run_chemgen(
        tutorial_dir, "FFCM2_model", "--custom-test", "custom_test.py",
        "--n-points-test", str(n_points), "--compile",
    )
    assert generation.returncode == 0, generation.stdout + generation.stderr

    run = run_binary(tutorial_dir)
    assert run.returncode == 0, run.stdout + run.stderr

    csv_path = tutorial_dir / "l2_norm_results.csv"
    assert csv_path.exists(), "chemgen binary did not produce l2_norm_results.csv"
    rows = read_l2_norm_csv(csv_path)
    assert len(rows) == n_points

    l2_values = [l2 for _temperature, l2 in rows]
    mean_l2 = sum(l2_values) / len(l2_values)
    outliers = [v for v in l2_values if v > OUTLIER_L2_TOLERANCE]

    assert mean_l2 < MEAN_L2_TOLERANCE, f"mean L2 error {mean_l2:.3g} exceeds {MEAN_L2_TOLERANCE:.3g}"
    outlier_fraction = len(outliers) / len(l2_values)
    assert outlier_fraction <= MAX_OUTLIER_FRACTION, (
        f"{len(outliers)}/{len(l2_values)} states ({outlier_fraction:.1%}) exceeded "
        f"{OUTLIER_L2_TOLERANCE:.0%} relative error, expected <= {MAX_OUTLIER_FRACTION:.0%}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("tutorial_dir", ["02_chemgen_error_assessment"], indirect=True)
def test_ffcm2_source_term_error_distribution_full_tutorial_size(tutorial_dir, cxx_compiler):
    """Same check at the tutorial's documented n_points=1000; run explicitly (`-m slow`)."""
    n_points = 1000
    generation = run_chemgen(
        tutorial_dir, "FFCM2_model", "--custom-test", "custom_test.py",
        "--n-points-test", str(n_points), "--compile",
    )
    assert generation.returncode == 0, generation.stdout + generation.stderr

    run = run_binary(tutorial_dir)
    assert run.returncode == 0, run.stdout + run.stderr

    rows = read_l2_norm_csv(tutorial_dir / "l2_norm_results.csv")
    assert len(rows) == n_points

    l2_values = [l2 for _temperature, l2 in rows]
    mean_l2 = sum(l2_values) / len(l2_values)
    outliers = [v for v in l2_values if v > OUTLIER_L2_TOLERANCE]

    assert mean_l2 < MEAN_L2_TOLERANCE
    assert len(outliers) / len(l2_values) <= MAX_OUTLIER_FRACTION
