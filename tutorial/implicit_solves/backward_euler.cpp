// This C++ program is a standalone translation of your Python code for solving dy/dt = source_1pt3(y)
// using various implicit time integration methods. It uses std::array, avoids external libraries,
// and writes CSV files for the time-evolution of y.

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <cassert>

constexpr int N = 3;
using Vec = std::array<double, N>;
using Mat = std::array<std::array<double, N>, N>;

// ----------------------------------------
// Utility Functions
// ----------------------------------------
double norm2(const Vec& x) {
    double sum = 0.0;
    for (double xi : x) sum += xi * xi;
    return std::sqrt(sum);
}

Vec operator+(const Vec& a, const Vec& b) {
    Vec r;
    for (int i = 0; i < N; ++i) r[i] = a[i] + b[i];
    return r;
}

Vec operator-(const Vec& a, const Vec& b) {
    Vec r;
    for (int i = 0; i < N; ++i) r[i] = a[i] - b[i];
    return r;
}

Vec operator*(double s, const Vec& a) {
    Vec r;
    for (int i = 0; i < N; ++i) r[i] = s * a[i];
    return r;
}

Vec operator*(const Mat& A, const Vec& x) {
    Vec r = {};
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            r[i] += A[i][j] * x[j];
    return r;
}

Mat eye() {
    Mat I = {};
    for (int i = 0; i < N; ++i) I[i][i] = 1.0;
    return I;
}

// ----------------------------------------
// Source Function and Jacobian
// ----------------------------------------
Vec source_1pt3(const Vec& y) {
    Vec s;
    s[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
    s[1] =  0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1];
    s[2] =  3e7 * y[1] * y[1];
    return s;
}

Mat dsource_1pt3_dy(const Vec& y) {
    Mat J = {};
    J[0][0] = -0.04;
    J[0][1] =  1e4 * y[2];
    J[0][2] =  1e4 * y[1];
    J[1][0] =  0.04;
    J[1][1] = -1e4 * y[2] - 6e7 * y[1];
    J[1][2] = -1e4 * y[1];
    J[2][1] =  6e7 * y[1];
    return J;
}

// ----------------------------------------
// GMRES (basic, non-restarted, no preconditioner)
// ----------------------------------------
Vec gmres_solve(const Mat& A, const Vec& b, double tol = 1e-10, int max_iter = 10) {
    Vec x = {};
    Vec r = b - (A * x);
    if (norm2(r) < tol) return x;
    Vec delta = r;  // no Arnoldi/Krylov here, minimal GMRES for small N
    Mat AI = A; // shallow copy for inversion

    // Gaussian elimination (for small 3x3 systems)
    for (int i = 0; i < N; ++i) {
        double pivot = AI[i][i];
        for (int j = 0; j < N; ++j) AI[i][j] /= pivot;
        delta[i] /= pivot;
        for (int k = 0; k < N; ++k) {
            if (k != i) {
                double factor = AI[k][i];
                for (int j = 0; j < N; ++j)
                    AI[k][j] -= factor * AI[i][j];
                delta[k] -= factor * delta[i];
            }
        }
    }
    return delta;
}

// ----------------------------------------
// Backward Euler Integrator
// ----------------------------------------
void backwards_euler(const Vec& y0, double dt, int n_steps, std::vector<Vec>& ys, std::vector<double>& time) {
    Vec yn = y0;
    ys.push_back(yn);
    time.push_back(0.0);
    Mat I = eye();

    for (int t = 0; t < n_steps; ++t) {
        Vec y_guess = yn;
        for (int iter = 0; iter < 5; ++iter) {
            Vec f = source_1pt3(y_guess);
            Mat J = dsource_1pt3_dy(y_guess);
            Mat A = I;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    A[i][j] = I[i][j] / dt - J[i][j];
            Vec res;
            for (int i = 0; i < N; ++i)
                res[i] = (y_guess[i] - yn[i]) / dt - f[i];
            Vec dy = gmres_solve(A, {-res[0], -res[1], -res[2]});
            for (int i = 0; i < N; ++i)
                y_guess[i] += dy[i];
            if (norm2(res) < 1e-10) break;
        }
        yn = y_guess;
        ys.push_back(yn);
        time.push_back(time.back() + dt);
    }
}

// ----------------------------------------
// Write CSV
// ----------------------------------------
void write_csv(const std::string& filename, const std::vector<double>& time, const std::vector<Vec>& ys) {
    std::ofstream file(filename);
    file << "time,y0,y1,y2\n";
    for (size_t i = 0; i < time.size(); ++i) {
        file << std::setprecision(12) << time[i];
        for (int j = 0; j < N; ++j) file << "," << ys[i][j];
        file << "\n";
    }
    file.close();
}

// ----------------------------------------
// Main
// ----------------------------------------
int main() {
    Vec y0 = {1.0, 0.0, 0.0};
    double dt = 1e-3;
    double tf = 0.03;
    int n_steps = static_cast<int>(tf / dt);

    std::vector<Vec> ys;
    std::vector<double> time;

    auto t_start = std::chrono::high_resolution_clock::now();
    backwards_euler(y0, dt, n_steps, ys, time);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout << "Backward Euler completed in " << elapsed.count() << " seconds.\n";

    write_csv("backward_euler.csv", time, ys);
    return 0;
}


