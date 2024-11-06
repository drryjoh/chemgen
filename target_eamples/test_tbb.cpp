#include <iostream>
#include <vector>
#include <chrono>
#include <tbb/tbb.h>

// Function to simulate work (replace these with your actual functions)
void function_work(int id) {
    // Simulate some computation
    for (volatile int i = 0; i < 1000000 * (id + 1); ++i);
}

void run_in_serial(const std::vector<void(*)(int)>& functions) {
    for (int i = 0; i < functions.size(); ++i) {
        functions[i](i);
    }
}

void run_in_parallel(const std::vector<void(*)(int)>& functions) {
    tbb::parallel_for(0, static_cast<int>(functions.size()), [&](int i) {
        functions[i](i);
    });
}

int main() {
    const int num_functions = 50;

    // Create a vector of 50 function pointers (could be different functions if needed)
    std::vector<void(*)(int)> functions(num_functions, function_work);

    // Measure serial execution time
    auto start_serial = std::chrono::high_resolution_clock::now();
    run_in_serial(functions);
    auto end_serial = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> serial_time = end_serial - start_serial;

    // Measure parallel execution time
    auto start_parallel = std::chrono::high_resolution_clock::now();
    run_in_parallel(functions);
    auto end_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_time = end_parallel - start_parallel;

    // Output results
    std::cout << "Serial execution time: " << serial_time.count() << " seconds\n";
    std::cout << "Parallel execution time: " << parallel_time.count() << " seconds\n";

    return 0;
}
