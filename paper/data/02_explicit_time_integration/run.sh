#!/bin/bash
cg
# run RK test case 
python3 ~/chemgen/bin/chemgen.py ffcm2_h2.yaml . --compile --custom-test custom_test.py
./bin/chemgen
./post_ct.py
