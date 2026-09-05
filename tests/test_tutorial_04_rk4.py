"""Regression test for tutorial/04_rk4.

Reproduces `chemgen.py ffcm2_h2.yaml . --custom-test custom_test.py --compile
--run-tests`, then independently integrates the same initial condition with
Cantera's IdealGasReactor (the same cross-check tutorial 04's post_ct.py
does by eye, as an assertion on the final temperature instead of a plot).

Note: this tutorial's custom_test.py reads test_configuration.yaml, but no
such file shipped in tutorial/04_rk4 -- post_ct.py's hardcoded initial
condition (T=1800K, P=101325Pa, O2:0.2 N2:0.6 H2:0.2) reveals what it should
contain, so a matching test_configuration.yaml was added alongside this test.
"""
import cantera as ct
import pytest

from conftest import run_binary, run_chemgen

INITIAL_STATE = dict(temperature=1800, pressure=101325.0, X="O2:0.2 N2:0.6 H2:0.2")


@pytest.mark.parametrize("tutorial_dir", ["04_rk4"], indirect=True)
def test_rk4_homogeneous_reactor_matches_cantera(tutorial_dir, cxx_compiler):
    generation = run_chemgen(tutorial_dir, "ffcm2_h2.yaml", "--custom-test", "custom_test.py", "--compile")
    assert generation.returncode == 0, generation.stdout + generation.stderr

    run = run_binary(tutorial_dir)
    assert run.returncode == 0, run.stdout + run.stderr

    trajectory = []
    for line in run.stdout.splitlines():
        parts = line.split()
        try:
            values = [float(p) for p in parts]
        except ValueError:
            continue
        if len(values) >= 2:
            trajectory.append(values)

    assert len(trajectory) == 40001, f"expected 40000 RK4 steps + initial state, got {len(trajectory)} rows"
    t0, T0 = trajectory[0][0], trajectory[0][1]
    t_final, T_final = trajectory[-1][0], trajectory[-1][1]
    assert t0 == pytest.approx(0.0, abs=1e-12)
    assert T0 == pytest.approx(INITIAL_STATE["temperature"], rel=1e-6)

    gas = ct.Solution(str(tutorial_dir / "ffcm2_h2.yaml"))
    gas.TPX = INITIAL_STATE["temperature"], INITIAL_STATE["pressure"], INITIAL_STATE["X"]
    reactor = ct.IdealGasReactor(gas)
    network = ct.ReactorNet([reactor])
    network.advance(t_final)

    rel_error = abs(T_final - reactor.T) / reactor.T
    assert rel_error < 0.02, (
        f"final temperature mismatch at t={t_final:.3e}s: "
        f"chemgen={T_final:.2f}K cantera={reactor.T:.2f}K rel_error={rel_error:.3%}"
    )
