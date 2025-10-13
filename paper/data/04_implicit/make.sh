g++ src/chemgen.cpp -o test_cantera \
  -I/Users/ryanjohnson/cantera-cpp/include \
  -L/Users/ryanjohnson/cantera-cpp/lib \
  -lcantera -lyaml-cpp -std=c++17 \
  -framework Accelerate -O3
