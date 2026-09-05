"""Shared fixtures for the ChemGen tutorial regression suite.

These tests treat each `tutorial/NN_*` directory as a living example: they
reproduce the exact commands documented in that tutorial's README inside an
isolated temp copy, compile the generated C++, run it, and assert on the
numeric results instead of asking a human to eyeball console output.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

try:
    import cantera as ct
except ImportError:  # pragma: no cover - guarded by requires_cantera below
    ct = None

REPO_ROOT = Path(__file__).resolve().parent.parent
CHEMGEN_PY = REPO_ROOT / "bin" / "chemgen.py"
TUTORIAL_DIR = REPO_ROOT / "tutorial"


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: exercises the full tutorial workload (longer runtime)")
    config.addinivalue_line("markers", "requires_pybind: needs pybind11 + a working extension build")
    config.addinivalue_line("markers", "requires_kokkos: needs the third_party/kokkos submodule built")
    config.addinivalue_line("markers", "requires_openfoam: needs a local OpenFOAM installation")


@pytest.fixture(scope="session", autouse=True)
def require_cantera_v3():
    if ct is None:
        pytest.skip("cantera is not installed (pip install cantera)")
    if ct.__version__.split(".")[0] != "3":
        pytest.skip(f"chemgen requires cantera >=3, found {ct.__version__}")


@pytest.fixture()
def cxx_compiler():
    """Return the compiler chemgen's own --compile path uses (see compile_and_run.py)."""
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        pytest.skip("no C++17 compiler (clang++/g++) found on PATH")
    return compiler


def copy_tutorial(tutorial_name, tmp_path):
    """Copy a tutorial directory into an isolated tmp dir, mirroring `cd tutorial/NN && ...`.

    chemgen reads `configuration.yaml` / `test_configuration.yaml` relative to the
    process cwd, so tests must run with cwd set to a directory that has these
    tutorial-local files alongside it -- exactly like a person following the README.
    """
    src = TUTORIAL_DIR / tutorial_name
    dst = tmp_path / tutorial_name
    shutil.copytree(src, dst)
    return dst


@pytest.fixture()
def tutorial_dir(request, tmp_path):
    """Parametrize with the tutorial directory name, e.g. `@pytest.mark.parametrize('tutorial_dir', ['01_mechanism_creation'], indirect=True)`."""
    return copy_tutorial(request.param, tmp_path)


def run_chemgen(cwd, mechanism, *extra_args, timeout=180):
    """Invoke bin/chemgen.py exactly as the tutorials do, with cwd = a tutorial copy."""
    cmd = [sys.executable, str(CHEMGEN_PY), mechanism, "."] + list(extra_args)
    result = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout,
    )
    return result


def run_binary(cwd, timeout=180):
    binary = cwd / "bin" / "chemgen"
    assert binary.exists(), f"expected compiled binary at {binary}, chemgen build must have failed"
    result = subprocess.run([str(binary)], cwd=cwd, capture_output=True, text=True, timeout=timeout)
    return result
