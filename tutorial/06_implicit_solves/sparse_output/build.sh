#!/bin/bash
MECH=isobutene-ic4h8-493-2716.yaml

#FIRST
python3 $HOME/code/chemgen/bin/chemgen.py ./mechanisms/$MECH . --custom-test custom_test.py --temperature-equation --ignore-other-species \
    --get-sparsity \
    --save-sparsity

#SECOND
python3 $HOME/code/chemgen/bin/chemgen.py ./mechanisms/$MECH . --custom-test custom_test.py --temperature-equation --ignore-other-species \
    --generate-permutation

#THIRD
python3 $HOME/code/chemgen/bin/chemgen.py ./mechanisms/$MECH . --custom-test custom_test.py --temperature-equation --ignore-other-species \
    --compile

rm -rf __pycache__

# #PURE COMPILING
# clang++ -std=c++23 -O2 -I/Users/jsampa/homebrew/opt/yaml-cpp/include -L/Users/jsampa/homebrew/opt/yaml-cpp/lib -lyaml-cpp -I$CONDA_PREFIX/include/eigen3 -o bin/bin/chemgen bin/src/chemgen.cpp
# ./bin/bin/chemgen