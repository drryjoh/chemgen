g++ -pg -g -O2 \
  -I/Users/jsampa/homebrew/opt/yaml-cpp/include \
  -L/Users/jsampa/homebrew/opt/yaml-cpp/lib -lyaml-cpp \
  -I$CONDA_PREFIX/include/eigen3 \
  bin/src/chemgen.cpp -o profile
