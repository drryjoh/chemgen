#!/usr/bin/env zsh
set -euo pipefail

# load aliases/functions defined in ~/.zshrc
[[ -f ~/.zshrc ]] && source ~/.zshrc
setopt aliases   # allow alias expansion in scripts

cg   # now this alias/function exists
export PATH="$(cd ~/chemgen/bin/ && pwd):$PATH"
chemgen.py FFCM2_model.yaml chemgen_conservative --pybind
chemgen.py FFCM2_model.yaml chemgen_temperature --pybind --temperature-equation
python3 ./setup.py build_ext --inplace
