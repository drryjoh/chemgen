pip3 install pybind11 setuptools


chemgen.py ffcm2_h2.yaml . --pybind --jacobian-temperature --force; python3 ./src/setup_chemgen.py build_ext --inplace 2>error
