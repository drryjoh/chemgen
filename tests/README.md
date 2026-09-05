# ChemGen tutorial regression tests

These tests turn the walkthroughs in [`tutorial/`](../tutorial/) into automated
regression checks: for each tutorial, a test runs the exact `chemgen.py`
command the README documents inside an isolated temp copy of that tutorial
directory, compiles the generated C++, runs it, and asserts on the numeric
result instead of asking a human to eyeball console output or a plot.

If a future change to the generator, a `.h.in` template, or a `configuration.yaml`
silently changes generated numerics, one of these should fail.

## Prerequisites

- Python 3.10+
- [Cantera](https://cantera.org/) >= 3.0 (chemgen hard-requires major version 3)
- A C++17 compiler on `PATH` (`clang++` or `g++` -- matches what
  `bin/modules/compile_and_run.py`'s `--compile` path uses)
- `pytest`

Optional, only needed for the tests that require them (see Markers below):
- `pybind11` + `setuptools` (tutorial 05)
- the `third_party/kokkos` submodule, built (tutorial 07)
- a local OpenFOAM 2412 installation (tutorial 10)

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate       # .venv\Scripts\activate on Windows
pip install -r requirements-dev.txt
```

`requirements-dev.txt` pulls in everything in `requirements.txt` (cantera,
pyyaml, matplotlib) plus the test-only extras (pytest, pybind11, seaborn).

## Running

```bash
# fast default suite: skips slow variants and anything needing kokkos/pybind/openfoam
pytest -m "not slow and not requires_kokkos and not requires_openfoam and not requires_pybind"

# everything that's runnable without special hardware/software
pytest -m "not requires_kokkos and not requires_openfoam"

# a single tutorial
pytest tests/test_tutorial_01_mechanism_creation.py -v

# the tutorial-02 test at its documented n_points=1000 (marked slow; ~20s either way,
# dominated by compile time, but this is the one to run before a release)
pytest -m slow
```

On this repo's reference machine (M-series Mac, clang++ from Xcode CLT) the
default fast suite (tutorials 01-05, `n_points=200` for 02) takes well under a
minute -- almost all of it is compiling the generated FFCM2 C++ once per test.

## Markers

Registered in `tests/conftest.py`:

| Marker | Meaning |
|---|---|
| `slow` | Runs the tutorial at its full documented workload (e.g. `n_points=1000` instead of 200). Same assertions, just closer to what the README shows. |
| `requires_pybind` | Needs `pybind11`/`setuptools` to build a Python extension (tutorial 05). |
| `requires_kokkos` | Needs `third_party/kokkos` built (tutorial 07). Currently a skipped placeholder -- see that file. |
| `requires_openfoam` | Needs a local OpenFOAM 2412 install (tutorial 10). Currently a skipped placeholder -- see that file. |

## What's covered, and what isn't yet

| Tutorial | Status |
|---|---|
| 01 mechanism creation | Full: source terms, internal energy, species cp/h/s/g vs Cantera; Newton-iteration convergence for temperature-from-energy. |
| 02 error assessment | Full: L2-norm error distribution (mean + outlier fraction) from `l2_norm_results.csv`, codifying the README's prose tolerance into assertions. |
| 03 decorators | Full: double/float/point configs all compile, produce the expected C++ signatures (scalar type, by-value vs const-ref), and still agree with Cantera numerically. |
| 04 rk4 | Full: RK4 homogeneous reactor cross-checked against an independently-run `cantera.IdealGasReactor` at the same final time. |
| 05 pybind | Full: builds the generated pybind11 extension and checks `chemgen.source(...)` against Cantera in-process. |
| 06 implicit solvers | Not implemented -- most involved tutorial (5 solvers x 2 linear solvers), and its `make.sh` hardcodes an absolute path to the author's machine. See the placeholder file for what's needed. |
| 07 kokkos | Skipped placeholder -- needs the `third_party/kokkos` submodule built first. |
| 08 derivatives, 09 large data | Not covered -- neither tutorial has a README; their `custom_test_*.py` scripts would need to be read and their intended usage documented before a runner can be written with confidence. |
| 10 target_cfd (OpenFOAM) | Skipped placeholder -- needs a local OpenFOAM 2412 install. |

## Bugs this test-writing pass found and fixed

Writing these tests surfaced two real issues in the tutorials themselves
(not in the tests) -- worth knowing about since they explain a couple of
non-obvious lines in the test code:

- `tutorial/03_decorators/configuration_point.yaml` had `jacobian_end: "n_species"`,
  which `chemgen.py`'s own validation rejects (must be `n_variables` or
  `n_species + 1`) -- this config could not generate code at all. Fixed to
  `n_variables`, matching its sibling configs.
- `tutorial/04_rk4/` was missing `test_configuration.yaml`, which its
  `custom_test.py` requires -- `post_ct.py`'s hardcoded initial condition
  (T=1800K, P=101325Pa, O2:0.2/N2:0.6/H2:0.2) revealed what it should
  contain, so that file was added.
- `tutorial/03_decorators/README.md` documented `chemgen.py one_reaction .
  --custom-test custom_test.py --compile --run`, but this tutorial ships no
  `custom_test.py` and `chemgen.py` has no `--run` flag (only `--run-tests`).
  Rewrote the "Precision" section to use the default test writer + the
  `test_configuration.yaml` already in that directory, which is what
  `test_tutorial_03_decorators.py` actually exercises.
- `tutorial/05_pybind/setup.py` and `test_pybind.py` referenced a
  `chemwrapper` module built from `chem.cpp`/`bindings.cpp`, none of which
  exist -- orphaned from before `create_pybind.py` settled on writing
  `src/setup_chemgen.py` and naming the real module `chemgen` (which is what
  that tutorial's README, and `test_tutorial_05_pybind.py`, actually use).
  Removed both files.

## Adding a new tutorial's test

1. Copy the pattern in `test_tutorial_04_rk4.py` (custom `--custom-test`
   script) or `test_tutorial_01_mechanism_creation.py` (default test writer).
2. Use `tutorial_dir` (indirect-parametrized with the tutorial's directory
   name) to get an isolated temp copy, and `run_chemgen`/`run_binary` from
   `conftest.py` to drive it -- both mirror the exact `cd tutorial/NN && chemgen.py ...`
   invocation the README shows, since chemgen reads `configuration.yaml` /
   `test_configuration.yaml` relative to the process's cwd.
3. Prefer asserting on a machine-readable artifact (a printed
   ChemGen-vs-Cantera pair, a CSV) over parsing a plot.
4. If the tutorial doesn't have a natural reference value, compute one
   in-test with Cantera (see `test_tutorial_04_rk4.py`'s `IdealGasReactor`
   cross-check) rather than hardcoding numbers from a README's example run.

## CI

See `.github/workflows/tests.yml` for a GitHub Actions workflow that runs the
fast default suite on every push/PR.
