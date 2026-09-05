"""Regression test for tutorial/01_mechanism_creation.

Reproduces `chemgen.py simple_mech . --compile --run-tests` and checks that
every ChemGen-vs-Cantera value pair the default test prints (source terms,
internal energy, species cp/h/s/g) still agrees within tolerance.
"""
import pytest

from conftest import run_binary, run_chemgen
from utils import max_relative_error, parse_comparison_pairs

# default_test.py compares ChemGen's generated thermo/kinetics fits against
# cantera computed live from the same mechanism; a couple of quantities carry
# a documented, small systematic difference from the polynomial fits used to
# build the monomial basis. See tutorial 02's README for the analogous
# source-term error discussion.
TOLERANCES = {
    "Source test result": 5e-2,
    "ChemGen internal energy": 1e-3,
    "Chemgen species cps": 5e-3,
    "Chemgen species enthalpies": 5e-3,
    "Chemgen species internal energies": 5e-3,
    "Chemgen species internal entropies": 5e-3,
    "Chemgen species gibbs energy": 5e-3,
    "Chemgen species gibbs reactions": 5e-2,
}


@pytest.mark.parametrize("tutorial_dir", ["01_mechanism_creation"], indirect=True)
def test_simple_mech_default_test(tutorial_dir, cxx_compiler):
    generation = run_chemgen(tutorial_dir, "simple_mech", "--compile", "--run-tests")
    assert generation.returncode == 0, generation.stdout + generation.stderr

    pairs = parse_comparison_pairs(generation.stdout)
    assert pairs, "no ChemGen/Cantera comparison pairs found in chemgen output"

    checked = set()
    for label, chemgen_values, cantera_values in pairs:
        tolerance = next((tol for prefix, tol in TOLERANCES.items() if label.startswith(prefix)), None)
        if tolerance is None:
            continue  # not a quantity we assert on (e.g. equilibrium constants)
        error = max_relative_error(chemgen_values, cantera_values)
        assert error < tolerance, f"{label}: max relative error {error:.3g} exceeds tolerance {tolerance:.3g}"
        checked.add(label)

    # every configured quantity's prefix must have matched at least one printed label
    for prefix in TOLERANCES:
        assert any(c.startswith(prefix) for c in checked), f"expected output containing '{prefix}' was not printed"


@pytest.mark.parametrize("tutorial_dir", ["01_mechanism_creation"], indirect=True)
def test_simple_mech_newton_iteration_converges(tutorial_dir, cxx_compiler):
    """The README's sample run shows temperature-from-energy Newton solve converging to ~1e-16 within 5 iterations."""
    generation = run_chemgen(tutorial_dir, "simple_mech", "--compile", "--run-tests")
    assert generation.returncode == 0, generation.stdout + generation.stderr

    residuals = {}
    for line in generation.stdout.splitlines():
        if "temperature_ for" in line and "iterations:" in line:
            it = int(line.split("for")[1].split("iterations")[0])
            residual = float(line.rsplit(":", 1)[1])
            residuals[it] = residual

    assert residuals, "no Newton iteration residual lines found"
    assert abs(residuals[max(residuals)]) < 1e-10, "Newton solve for temperature-from-energy did not converge"
