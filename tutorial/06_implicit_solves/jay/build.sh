#!/bin/bash

#JAYS VERSION
python3 /Users/jsampa/Code/chemgen/bin/chemgen.py FFCM2_model.yaml bin/ --get-sparsity --plot-sparsity --skip-tests

# #PURE COMPILING
# BREW_PREFIX=$(brew --prefix yaml-cpp)
# clang++ -std=c++23 -O2 -I"${BREW_PREFIX}/include" -L"${BREW_PREFIX}/lib" -lyaml-cpp -I$CONDA_PREFIX/include/eigen3 chemgen.cpp -o main