pip3 install pybind11 setuptools



chemgen.py ffcm2_h2.yaml . --pybind; python3 ./src/setup_chemgen.py build_ext --inplace 2>error

chemgen.py ffcm2_h2.yaml . --pybind --ignore-temp-dependence; python3 ./src/setup_chemgen.py build_ext --inplace 2>error
