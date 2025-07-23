#!/bin/bash

clang++ -c mylibrary.cpp -o mylibrary.o #2>&1 >/dev/null
# clang++ -c mylibrary_cpp.h -o mylibrary.o # Warning
echo "Done compiling mylibrary"
ar rcs libmylibrary.a mylibrary.o
echo "Done creating static library"
clang++ main.cpp -L. -lmylibrary -o myprogram
echo "Done compiling executable"
./myprogram
