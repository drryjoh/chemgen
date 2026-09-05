"""Placeholder for tutorial/07_kokkos_parallelism.

Requires building the third_party/kokkos submodule out-of-band (CMake,
Kokkos_ENABLE_THREADS=ON) before test_kokkos.cpp can even compile. Skipped
by default; set CHEMGEN_TEST_KOKKOS=1 once the submodule build step is
scripted (this repo checks out `.gitmodules` but the submodule isn't
initialized by default -- `git submodule update --init third_party/kokkos`
first).
"""
import os

import pytest

pytestmark = pytest.mark.requires_kokkos


@pytest.mark.skipif(
    os.environ.get("CHEMGEN_TEST_KOKKOS") != "1",
    reason="set CHEMGEN_TEST_KOKKOS=1 on a machine with third_party/kokkos built",
)
def test_placeholder():
    pytest.skip("not yet implemented -- see module docstring")
