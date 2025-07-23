#pragma once
#include <iostream>
#include <array>
#include <random>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <algorithm> 
#include <cstddef>
#include <Eigen/Dense>

// template<typename Scalar>
// using activationFunction = void(*)(Scalar&, Scalar, Scalar);


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template<typename Scalar, int input_size, int output_size, typename ActFun>
inline void Dense_MLP_LU(Scalar* __restrict outputs, const Scalar* __restrict inputs, const Scalar * __restrict weights, const Scalar * __restrict biases, ActFun activation_function, Scalar alpha) noexcept 
{
    // Fixed-size Eigen types
    using InputVec = Eigen::Matrix<Scalar, input_size, 1>;
    using OutputVec = Eigen::Matrix<Scalar, output_size, 1>;
    using WeightMat = Eigen::Matrix<Scalar, input_size, output_size, Eigen::RowMajor>;
    // Map inputs, weights, biases, outputs
    Eigen::Map<const InputVec> input_vec(inputs);
    Eigen::Map<const WeightMat> weight_matrix(weights);
    Eigen::Map<const OutputVec> bias_vec(biases);
    Eigen::Map<OutputVec> output_vec(outputs);
    // Perform fast multiply-add with no temporaries
    output_vec.noalias() = weight_matrix.transpose() * input_vec + bias_vec;
    // Apply activation with SIMD
    #pragma omp simd
    for(int i = 0; i < output_size; ++i) {
        activation_function(outputs[i], outputs[i], alpha);
    }
}

//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template <typename Scalar = double>
inline auto MLP_LU(const std::array<Scalar, 3256>& initial_input) {

    // Eigen type aliases for better performance
    using VectorXs = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MatrixXs = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using Vector3256s = Eigen::Matrix<Scalar, 3256, 1>;
    using Vector9409s = Eigen::Matrix<Scalar, 9409, 1>;

    // Dense layer 1
    constexpr std::array<Scalar, 3256> weights_1 = {};
    constexpr std::array<Scalar, 1> biases_1 = {};

    // Dense layer 3
    constexpr std::array<Scalar, 1> weights_3 = {};
    constexpr std::array<Scalar, 1> biases_3 = {};

    // Dense layer 5
    constexpr std::array<Scalar, 9409> weights_5 = {};
    constexpr std::array<Scalar, 9409> biases_5 = {};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 

    constexpr std::array<Scalar, 3256> input_norm_std = {};

    constexpr std::array<Scalar, 3256> input_min_mean = {};

    // Final output
    constexpr static std::array<Scalar, 9409> output_norm_std = {};
    constexpr static std::array<Scalar, 9409> output_min_mean = {};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept 
    {
        output = input;
    };

    static constexpr Scalar kC0 = 0.044715;
    static constexpr Scalar kSqrt2PiInv = Scalar(0.7978845608028654);
    auto gelu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept 
    {
        Scalar x3 = input * input * input;
        Scalar y  = kSqrt2PiInv * (input + kC0 * x3);
        output     = Scalar(0.5) * input * (Scalar(1) + std::tanh(y));
    };

    auto nonzero_diag_activation = +[](Scalar& output, Scalar input, int index) noexcept
    {
        constexpr int M = 97;
        constexpr int FLAT_DIM = M * M;
        constexpr Scalar EPS = Scalar(1e-16);
        
        bool is_diagonal = (index % (M + 1)) == 0;
        
        if (is_diagonal) {
            Scalar abs_x = std::abs(input);
            Scalar sign_x = (input >= Scalar(0)) ? Scalar(1) : Scalar(-1);
            if (input == Scalar(0)) {
                sign_x = Scalar(1);
            }
            output = sign_x * std::max(abs_x, EPS);
        } else {
            output = input;
        }
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 

    
    // model input and flattened
    constexpr int flat_size = 3256; 
    std::array<Scalar, flat_size> model_input;

    // normalize input using Eigen operations
    Eigen::Map<const Vector3256s> input_vec(initial_input.data());
    Eigen::Map<const Vector3256s> input_min_vec(input_min_mean.data());
    Eigen::Map<const Vector3256s> input_std_vec(input_norm_std.data());
    Eigen::Map<Vector3256s> model_input_vec(model_input.data());
    
    model_input_vec = (input_vec - input_min_vec).cwiseQuotient(input_std_vec); 

    if (model_input.size() != 3256) { throw std::invalid_argument("Invalid input size. Expected size: 3256"); }

    // Dense, layer 1
    static std::array<Scalar, 1> layer_1_output;
    Dense_MLP_LU<Scalar, 3256, 1>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        linear, 0.0);

    // Activation, layer 2
    static std::array<Scalar, 1> layer_2_output;
    for (int i = 0; i < 1; ++i) {
        gelu(layer_2_output[i], layer_1_output[i], 0.0);
    }

    // Dense, layer 3
    static std::array<Scalar, 1> layer_3_output;
    Dense_MLP_LU<Scalar, 1, 1>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        linear, 0.0);

    // Activation, layer 4
    static std::array<Scalar, 1> layer_4_output;
    for (int i = 0; i < 1; ++i) {
        gelu(layer_4_output[i], layer_3_output[i], 0.0);
    }

    // Dense, layer 5
    static std::array<Scalar, 9409> layer_5_output;
    Dense_MLP_LU<Scalar, 1, 9409>(
        layer_5_output.data(), layer_4_output.data(),
        weights_5.data(), biases_5.data(),
        linear, 0.0);

    // Activation, layer 6
    static std::array<Scalar, 9409> layer_6_output;
    // Apply nonzero_diag_activation - this needs to be done element-wise due to index dependency
    for (int i = 0; i < 9409; ++i) {
        nonzero_diag_activation(layer_6_output[i], layer_5_output[i], i);
    }


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    static std::array<Scalar, 9409> model_output;

    // denormalize output using Eigen operations
    Eigen::Map<const Vector9409s> layer_6_vec(layer_6_output.data());
    Eigen::Map<const Vector9409s> output_std_vec(output_norm_std.data());
    Eigen::Map<const Vector9409s> output_mean_vec(output_min_mean.data());
    Eigen::Map<Vector9409s> model_output_vec(model_output.data());
    
    model_output_vec = layer_6_vec.cwiseProduct(output_std_vec) + output_mean_vec;

    return model_output;

}