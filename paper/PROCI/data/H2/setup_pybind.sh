#!/usr/bin/env zsh
set -euo pipefail

# load aliases/functions defined in ~/.zshrc
[[ -f ~/.zshrc ]] && source ~/.zshrc
setopt aliases   # allow alias expansion in scripts

cg   # now this alias/function exists
export PATH="$(cd ~/chemgen/bin/ && pwd):$PATH"
chemgen.py ffcm2_h2.yaml chemgen_conservative --pybind
chemgen.py ffcm2_h2.yaml chemgen_temperature --pybind --temperature-equation
chemgen.py ffcm2_h2.yaml chemgen_dtemperature_ignore --pybind --ignore-temp-dependence
chemgen.py ffcm2_h2.yaml chemgen_conservative_ignore_species --pybind --ignore-other-species 
chemgen.py ffcm2_h2.yaml chemgen_temperature_ignore_species --pybind --temperature-equation --ignore-other-species 
python3 ./setup.py build_ext --inplace
