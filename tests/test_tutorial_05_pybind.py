"""Regression test for tutorial/05_pybind.

Reproduces `chemgen.py ffcm2_h2.yaml . --pybind`, builds the generated
pybind11 extension exactly as the README instructs
(`python3 ./src/setup_chemgen.py build_ext --inplace`), then calls
`chemgen.source(...)` from Python and checks it against Cantera's
net_production_rates -- the same comparison tutorial 01/02 do for the C++
entry point, but through the Python binding instead.

Note: this tutorial directory also ships a `setup.py` and `test_pybind.py`
that reference a `chemwrapper` module built from `chem.cpp`/`bindings.cpp` --
neither file exists. Those two files are stale leftovers unrelated to the
actual generation path (`create_pybind.py` writes `setup_chemgen.py`, and the
real module is named `chemgen`, per the README's own build command); this
test does not use them.
"""
import subprocess
import sys
import textwrap

import cantera as ct
import pytest

from conftest import run_chemgen

pytest.importorskip("pybind11", reason="pybind11 not installed (pip install pybind11 setuptools)")

STATE = dict(temperature=1600, pressure=101325.0, X="H2:0.2, O2:0.1, N2:0.7")


@pytest.mark.requires_pybind
@pytest.mark.parametrize("tutorial_dir", ["05_pybind"], indirect=True)
def test_pybind_source_matches_cantera(tutorial_dir):
    generation = run_chemgen(tutorial_dir, "ffcm2_h2.yaml", "--pybind")
    assert generation.returncode == 0, generation.stdout + generation.stderr

    setup_file = tutorial_dir / "src" / "setup_chemgen.py"
    assert setup_file.exists(), "chemgen --pybind did not generate src/setup_chemgen.py"

    build = subprocess.run(
        [sys.executable, str(setup_file), "build_ext", "--inplace"],
        cwd=tutorial_dir, capture_output=True, text=True, timeout=180,
    )
    assert build.returncode == 0, build.stdout + build.stderr

    gas = ct.Solution(str(tutorial_dir / "ffcm2_h2.yaml"))
    gas.TPX = STATE["temperature"], STATE["pressure"], STATE["X"]
    cantera_source = [float(v) for v in gas.net_production_rates]
    concentrations = [float(v) for v in gas.concentrations]

    check_script = textwrap.dedent(f"""
        import sys
        sys.path.insert(0, ".")
        import chemgen as cg
        concentrations = {concentrations!r}
        temperature = {STATE["temperature"]!r}
        result = cg.source(concentrations, temperature)
        print(",".join(repr(v) for v in result))
    """)
    check = subprocess.run(
        [sys.executable, "-c", check_script],
        cwd=tutorial_dir, capture_output=True, text=True, timeout=60,
    )
    assert check.returncode == 0, check.stdout + check.stderr

    chemgen_source = [float(v) for v in check.stdout.strip().splitlines()[-1].split(",")]
    # cg.source() binds the ChemicalState-returning overload: index 0 is the
    # energy/temperature source term, followed by the n_species production
    # rates (in the same order as Cantera's net_production_rates).
    assert len(chemgen_source) == len(cantera_source) + 1
    chemgen_source = chemgen_source[1:]
    errors = [
        abs(cg_v - ct_v) / (abs(ct_v) if abs(ct_v) > 1e-8 else 1.0)
        for cg_v, ct_v in zip(chemgen_source, cantera_source)
    ]
    assert max(errors) < 5e-2, f"pybind source() vs Cantera net_production_rates: max relative error {max(errors):.3g}"
