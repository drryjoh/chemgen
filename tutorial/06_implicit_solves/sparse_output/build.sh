#!/bin/bash
# MECH=hydrogen-10-28.yaml
# MECH=jetA-detailed-NOx-203-1589.yaml
MECH=ic8-874-6864.yaml
# MECH=md-nc7-3787-10264.yaml
# MECH=mmc5-7171-38324.yaml

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

# #compile only
# g++ -std=c++23 -O3 -ffast-math -march=native -mtune=native \
#   -funroll-loops -fvectorize -fslp-vectorize \
#   -I$CONDA_PREFIX/include \
#   -L$CONDA_PREFIX/lib -lyaml-cpp \
#   -I$CONDA_PREFIX/include/eigen3 \
#   -Wl,-rpath,$CONDA_PREFIX/lib \
#   -L./src -lchemgen_library \
#   -framework Accelerate \
#   -Wno-deprecated-declarations \
#   -Wno-nan-infinity-disabled \
#   -o ./bin/chemgen ./src/chemgen.cpp
# ./bin/chemgen