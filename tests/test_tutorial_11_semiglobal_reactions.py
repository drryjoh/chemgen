"""Regression test for tutorial/11_semiglobal_reactions.

Reproduces `chemgen.py semiglobal_butadiene.yaml . --custom-test
custom_test.py --compile`, then checks:

1. The default-test source term (species production) matches Cantera for
   every species, including ones whose production comes entirely from
   reactions where they aren't part of the rate law's `orders` (CO2 is a
   product of reaction 1 but not in its `orders`; H2O only ever modifies
   reaction 2's rate without being consumed/produced by it). This is the
   case that needs true reactant/product stoichiometry and the semi-global
   rate law's `orders` tracked separately -- see process_reactions.py's
   get_rate_law_reactants().
2. RK4 (small explicit steps) and SDIRK2 (100x larger implicit steps, which
   exercises source_jacobian's Newton solve on this same semi-global
   mechanism) reach the same final time and agree with each other and with
   an independently-run Cantera IdealGasReactor.
"""
import cantera as ct
import pytest

from conftest import run_binary, run_chemgen
from utils import max_relative_error, parse_comparison_pairs

STATE = dict(temperature=1400, pressure=101325.0, X="C4H6:0.1, O2:0.5, CO:0.1, H2O:0.1, CO2:0.2")

CORE_TOLERANCES = {
    "Source test result": 5e-2,
}


@pytest.mark.parametrize("tutorial_dir", ["11_semiglobal_reactions"], indirect=True)
def test_semiglobal_source_term_matches_cantera(tutorial_dir, cxx_compiler):
    generation = run_chemgen(tutorial_dir, "semiglobal_butadiene.yaml", "--compile", "--run-tests")
    assert generation.returncode == 0, generation.stdout + generation.stderr

    pairs = parse_comparison_pairs(generation.stdout)
    assert pairs

    checked = set()
    for label, chemgen_values, cantera_values in pairs:
        tolerance = next((tol for prefix, tol in CORE_TOLERANCES.items() if label.startswith(prefix)), None)
        if tolerance is None:
            continue
        error = max_relative_error(chemgen_values, cantera_values)
        assert error < tolerance, f"{label}: max relative error {error:.3g} exceeds tolerance {tolerance:.3g}"
        checked.add(label)
    for prefix in CORE_TOLERANCES:
        assert any(c.startswith(prefix) for c in checked), f"expected output containing '{prefix}'"


@pytest.mark.parametrize("tutorial_dir", ["11_semiglobal_reactions"], indirect=True)
def test_rk4_and_sdirk2_agree_and_match_cantera(tutorial_dir, cxx_compiler):
    generation = run_chemgen(tutorial_dir, "semiglobal_butadiene.yaml", "--custom-test", "custom_test.py", "--compile")
    assert generation.returncode == 0, generation.stdout + generation.stderr

    run = run_binary(tutorial_dir)
    assert run.returncode == 0, run.stdout + run.stderr

    rk4_rows, sdirk2_rows = [], []
    for line in run.stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "RK4":
            rk4_rows.append([float(v) for v in parts[1:]])
        elif parts[0] == "SDIRK2":
            sdirk2_rows.append([float(v) for v in parts[1:]])

    assert len(rk4_rows) == 3001  # 3000 steps + initial state
    assert len(sdirk2_rows) == 31  # 30 steps + initial state

    t_rk4, T_rk4 = rk4_rows[-1][0], rk4_rows[-1][1]
    t_sdirk2, T_sdirk2 = sdirk2_rows[-1][0], sdirk2_rows[-1][1]
    assert t_rk4 == pytest.approx(t_sdirk2, rel=1e-9)

    rel_error_solvers = abs(T_rk4 - T_sdirk2) / T_sdirk2
    assert rel_error_solvers < 1e-3, (
        f"RK4 (T={T_rk4:.2f}K) and SDIRK2 (T={T_sdirk2:.2f}K) disagree by {rel_error_solvers:.3%} "
        f"at t={t_rk4:.3e}s despite integrating the same semi-global mechanism from the same state"
    )

    gas = ct.Solution(str(tutorial_dir / "semiglobal_butadiene.yaml"))
    gas.TPX = STATE["temperature"], STATE["pressure"], STATE["X"]
    reactor = ct.IdealGasReactor(gas)
    network = ct.ReactorNet([reactor])
    network.advance(t_rk4)

    for name, T_chemgen in (("RK4", T_rk4), ("SDIRK2", T_sdirk2)):
        rel_error = abs(T_chemgen - reactor.T) / reactor.T
        assert rel_error < 0.01, (
            f"{name} final temperature mismatch at t={t_rk4:.3e}s: "
            f"chemgen={T_chemgen:.2f}K cantera={reactor.T:.2f}K rel_error={rel_error:.3%}"
        )
