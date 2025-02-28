#include <Kokkos_Core.hpp>
#include <iostream>

// Expensive function to be run on the GPU
KOKKOS_INLINE_FUNCTION
double expensive_function(int i) {
    double result = 0.0;
    for (int j = 0; j < 1000; ++j) {
        result += 1.0 / (i + j + 1.0);
    }
    return result;
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv); // Initialize Kokkos runtime
    {
        int N = 1000000; // Number of elements to process

        // Allocate a result array on the host and device
        Kokkos::View<double*> results("results", N);

        // Launch parallel kernel on the GPU
        Kokkos::parallel_for("ExpensiveFunctionKernel", N, KOKKOS_LAMBDA(int i) {
            results(i) = expensive_function(i); // Run the function in parallel
        });

        // Optional: Transfer data back to CPU and summarize results
        double total_sum = 0.0;
        Kokkos::parallel_reduce("SumResults", N, KOKKOS_LAMBDA(int i, double& sum) {
            sum += results(i); // Sum results on the GPU
        }, total_sum);

        // Output the result on the CPU
        std::cout << "Total sum of results: " << total_sum << std::endl;
    }
    Kokkos::finalize(); // Finalize Kokkos runtime
    return 0;
}

