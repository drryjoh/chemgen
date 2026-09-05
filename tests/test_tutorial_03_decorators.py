"""Regression test for tutorial/03_decorators.

Note: the tutorial's README documents `chemgen.py one_reaction . --custom-test
custom_test.py --compile --run`, but this tutorial directory ships no
custom_test.py and chemgen.py has no `--run` flag (only `--run-tests`) --
that command is stale. This test instead drives the decorator configs
(double/float/point) through the default test writer, which is what the
mechanism + test_configuration.yaml here actually support, and checks both
that the decorator substitution produced the right C++ (float vs double,
pass-by-value vs const-ref) and that the numeric results still agree with
Cantera.
"""
import shutil

import pytest

from conftest import run_binary, run_chemgen
from utils import max_relative_error, parse_comparison_pairs

CORE_TOLERANCES = {
    "Source test result": 5e-2,
    "ChemGen internal energy": 1e-3,
    "Chemgen species cps": 5e-3,
    "Chemgen species enthalpies": 5e-3,
    "Chemgen species internal energies": 5e-3,
    "Chemgen species internal entropies": 5e-3,
    "Chemgen species gibbs energy": 5e-3,
}


@pytest.mark.parametrize("tutorial_dir", ["03_decorators"], indirect=True)
@pytest.mark.parametrize("config_name,expected_scalar,expected_by_value", [
    ("configuration_double.yaml", "double", False),
    ("configuration_float.yaml", "float", False),
    ("configuration_point.yaml", "double", True),
])
def test_decorator_substitution(tutorial_dir, cxx_compiler, config_name, expected_scalar, expected_by_value):
    shutil.copy(tutorial_dir / config_name, tutorial_dir / "configuration.yaml")

    generation = run_chemgen(tutorial_dir, "one_reaction", "--compile", "--run-tests")
    assert generation.returncode == 0, generation.stdout + generation.stderr

    reactions_h = (tutorial_dir / "src" / "reactions.h").read_text()
    assert f"{expected_scalar} call_forward_reaction_0(" in reactions_h

    by_value_signature = f"({expected_scalar} temperature, {expected_scalar} log_temperature)"
    by_ref_signature = f"(const {expected_scalar}& temperature, const {expected_scalar}& log_temperature)"
    if expected_by_value:
        assert by_value_signature in reactions_h
        assert by_ref_signature not in reactions_h
    else:
        assert by_ref_signature in reactions_h
        assert by_value_signature not in reactions_h

    pairs = parse_comparison_pairs(generation.stdout)
    assert pairs
    tolerance_scale = 5.0 if expected_scalar == "float" else 1.0  # single precision needs more slack
    checked_prefixes = set()
    for label, chemgen_values, cantera_values in pairs:
        tolerance = next((tol for prefix, tol in CORE_TOLERANCES.items() if label.startswith(prefix)), None)
        if tolerance is None:
            continue
        error = max_relative_error(chemgen_values, cantera_values)
        assert error < tolerance * tolerance_scale, (
            f"[{config_name}] {label}: max relative error {error:.3g} exceeds tolerance {tolerance * tolerance_scale:.3g}"
        )
        checked_prefixes.add(label)
    for prefix in CORE_TOLERANCES:
        assert any(c.startswith(prefix) for c in checked_prefixes), f"expected output containing '{prefix}'"
