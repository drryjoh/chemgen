#!/bin/bash

clang++ -std=c++17 -I/Users/eching/homebrew/opt/yaml-cpp/include -L/Users/eching/homebrew/opt/yaml-cpp/lib -lyaml-cpp -I/Users/eching/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -c src/chemgen_library.cpp -o src/chemgen_library.o
echo "Done compiling mylibrary"

ar rcs src/libchemgen_library.a src/chemgen_library.o
echo "Done creating static library"

clang++ -std=c++17 -I/Users/eching/homebrew/opt/yaml-cpp/include -L/Users/eching/homebrew/opt/yaml-cpp/lib -lyaml-cpp -I/Users/eching/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -L./src -lchemgen_library -o ./bin/chemgen src/chemgen.cpp
echo "Done compiling executable"
