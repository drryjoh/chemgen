#!/bin/bash

#JAYS VERSION
python3 /Users/jsampa/Code/chemgen/bin/chemgen.py FFCM2_model.yaml bin/ --custom-test custom_test.py --temperature-equation --ignore-other-species \
    --compile \
    --skip \
    # --get-sparsity \

# #PURE COMPILING
# clang++ -std=c++23 -O2 -I/Users/jsampa/homebrew/opt/yaml-cpp/include -L/Users/jsampa/homebrew/opt/yaml-cpp/lib -lyaml-cpp -I$CONDA_PREFIX/include/eigen3 -o bin/bin/chemgen bin/src/chemgen.cpp
# ./bin/bin/chemgen