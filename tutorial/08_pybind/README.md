pip3 install pybind11 setuptools
chemgen.py andrea.yaml . --pybind --remove_reactions
python3 ./src/setup_chemgen.py build_ext --inplace 2>error
