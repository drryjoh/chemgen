#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <iomanip> 
#include "SMALLEST.hpp" 
#include <chrono>

using namespace std;
using Scalar = double;

int main() {
    
    ifstream file("../data_200/A_199.csv");
    if (!file) {
        std::cerr << "Error opening file A.csv" << std::endl;
        return 1;
    }

    // load A into 2d array matrix and flatten
    std::array<Scalar, 96*96> A;
    for (int i = 0; i < 96; ++i) {
        for (int j = 0; j < 96; ++j) {
            file >> A[(i*96)+j];
            if (file.peek() == ',') file.ignore(); // Skip comma
        }
    }

    // Pass the input to your generated CNN function
    // auto start = std::chrono::high_resolution_clock::now();
    auto output = SMALLEST<Scalar>(A);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = end - start;
    // std::cout << "NN time = " << elapsed.count() << " seconds" << std::endl; 

    // Print the results with high precision
    std::cout << std::scientific << std::setprecision(15);  // Set precision and scientific notation
    std::cout << "Output:\n";  // Print each value on a new line
    for(int row = 0; row < 96; ++row) {
        for(int col = 0; col < 96; ++col) {
            std::cout << output[(row*96)+col] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << std::endl;

    return 0;
}

/*
Compile and run:
clang++ -std=c++23 -Wall -O3 -march=native -o cpp_outputs cpp_outputs.cpp
./test
*/