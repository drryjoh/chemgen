#!/bin/bash

#JAYS VERSION
python3 $HOME/code/chemgen/bin/chemgen.py FFCM2_model.yaml bin/ --custom-test custom_test.py --temperature-equation --ignore-other-species \
    --get-sparsity --save-sparsity \
    --compile \
    # --skip \
    # > error.txt 2>&1
    # --get-sparsity \
    # --print-sparsity \

# #PURE COMPILING
# clang++ -std=c++23 -O2 -I/Users/jsampa/homebrew/opt/yaml-cpp/include -L/Users/jsampa/homebrew/opt/yaml-cpp/lib -lyaml-cpp -I$CONDA_PREFIX/include/eigen3 -o bin/bin/chemgen bin/src/chemgen.cpp
# ./bin/bin/chemgen