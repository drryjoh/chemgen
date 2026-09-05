"""Placeholder for tutorial/06_implicit_solves.

This tutorial exercises SDIRK-2/4, Rosenbrock, YASS and Backward Euler with
GMRES/direct linear solvers, an order-of-accuracy convergence study, and a
timing comparison table. It's the most involved tutorial in the repo and its
`make.sh` hardcodes an absolute path to the author's local Cantera C++
install, so it isn't runnable as-is on another machine.

Left unimplemented deliberately rather than guessing at solver tolerances --
follow the pattern in test_tutorial_04_rk4.py (generate/compile/run, then
cross-check the final state against an independently-run Cantera reactor)
once someone can validate expected tolerances per solver against a real run.
"""
import pytest


@pytest.mark.skip(reason="not yet implemented -- see module docstring")
def test_placeholder():
    pass
