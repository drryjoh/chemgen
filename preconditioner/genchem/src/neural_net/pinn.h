#include <array>
#include <random>
#include <iostream>
#include <cmath>

// train (solve) pinn
auto train_pinn(const SpeciesJacobian &A, const Species &b){

    Species x;
    int N = 96;

    // intialize weights
    std::mt19937_64 rng(1234);
    std::uniform_real_distribution<double> unif(-0.1,0.1);
    for(int i=0;i<N;++i)
        x[i] = unif(rng);

    // hyperparams
    const double lr     = 1e-3;
    const int    epochs = 5000;

    // buffers
    Species r;      
    Species grad; 

    for(int e = 0; e < epochs; ++e){

        // compute residual
        for(int i = 0; i < N; ++i){
            double s = 0.0;
            for(int j = 0; j < N; ++j)
                s += A[i][j] * x[j];
            r[i] = s - b[i];
        }

        // monitor loss (optional)
        if((e+1) % 500 == 0){
            double loss = 0;
            for(int i=0;i<N;++i) loss += r[i]*r[i];
            loss /= N;
            std::cout << "Epoch " << (e+1)
                      << " | loss=" << loss << "\n";
        }

        // back prop
        for(int k = 0; k < N; ++k){
            double g = 0.0;
            for(int i = 0; i < N; ++i)
                g += r[i] * A[i][k];
            grad[k] = 2.0 * g;
        }

        // update
        for(int k = 0; k < N; ++k)
            x[k] -= lr * grad[k];
    }

    return x;
}

