#!/bin/bash
MECH=jetA-detailed-NOx-203-1589.yaml

#FIRST
python3 $HOME/code/chemgen/bin/chemgen.py ./mechanisms/$MECH . --custom-test custom_test.py --temperature-equation --ignore-other-species \
    --get-sparsity \
    --save-sparsity \
    --skip-tests

#SECOND
python3 $HOME/code/chemgen/bin/chemgen.py ./mechanisms/$MECH . --custom-test custom_test.py --temperature-equation --ignore-other-species \
    --generate-permutation \
    --save-permutation \
    --save-output-sparsity \
    --save-input-sparsity \
    --compile

# #OTHER
# python3 $HOME/code/chemgen/bin/chemgen.py ./mechanisms/$MECH . --custom-test custom_test.py --temperature-equation --ignore-other-species \
#     --skip \
#     --compile

rm -rf __pycache__

# #PURE COMPILING
# clang++ -std=c++23 -O2 -I/Users/jsampa/homebrew/opt/yaml-cpp/include -L/Users/jsampa/homebrew/opt/yaml-cpp/lib -lyaml-cpp -I$CONDA_PREFIX/include/eigen3 -o bin/bin/chemgen bin/src/chemgen.cpp
# ./bin/bin/chemgen