#!/bin/bash

BREW_PREFIX=$(brew --prefix yaml-cpp)
clang++ -std=c++23 -O2 -I"${BREW_PREFIX}/include" -L"${BREW_PREFIX}/lib" -lyaml-cpp -I$CONDA_PREFIX/include/eigen3 chemgen.cpp -o main