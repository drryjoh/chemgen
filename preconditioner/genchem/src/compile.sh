#!/bin/bash
# generate and compile
# SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# export PATH="$(${SCRIPT_DIR}/../../bin && pwd):$PATH"
# export DYLD_LIBRARY_PATH=/Users/jsampa/miniconda3/envs/chemgen/lib:$DYLD_LIBRARY_PATH
# chemgen.py FFCM2_model.yaml . --custom-test custom_test.py --compile

# DO NOT TOUCH AT ALL
# exec 2>/dev/null
rm -rf test
export DYLD_LIBRARY_PATH=/Users/jsampa/miniconda3/envs/chemgen/lib:$DYLD_LIBRARY_PATH
clang++ -O3 -march=native -I/Users/jsampa/miniconda3/envs/chemgen/include -L/Users/jsampa/miniconda3/envs/chemgen/lib -lyaml-cpp  chemgen.cpp -o test
./test