#!/bin/bash
## o3-mini
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
export PATH="$(${SCRIPT_DIR}/../../bin && pwd):$PATH"
# export PATH="$(cd ../../bin && pwd):$PATH"
export DYLD_LIBRARY_PATH=/Users/jsampa/miniconda3/envs/chemgen/lib:$DYLD_LIBRARY_PATH
chemgen.py FFCM2_model.yaml . --custom-test custom_test.py --compile
# /Users/jsampa/Code/chemgen/tutorial/06_implicit_solves/bin/chemgen