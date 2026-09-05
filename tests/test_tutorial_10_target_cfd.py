"""Placeholder for tutorial/10_target_cfd.

Requires a local OpenFOAM 2412 installation and its `wmake` build system --
decoupled from ChemGen's own --compile/--cmake paths entirely. Skipped by
default; set CHEMGEN_TEST_OPENFOAM=1 on a machine with OpenFOAM 2412 sourced.
"""
import os

import pytest

pytestmark = pytest.mark.requires_openfoam


@pytest.mark.skipif(
    os.environ.get("CHEMGEN_TEST_OPENFOAM") != "1",
    reason="set CHEMGEN_TEST_OPENFOAM=1 on a machine with OpenFOAM 2412 installed",
)
def test_placeholder():
    pytest.skip("not yet implemented -- see module docstring")
