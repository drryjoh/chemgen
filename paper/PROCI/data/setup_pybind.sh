#!/bin/zsh
cg
chemgen.py ffcm2_h2.yaml chemgen_conservative --pybind 
chemgen.py ffcm2_h2.yaml chemgen_temperature --pybind  
python3 ./setup.py build_ext --inplace
