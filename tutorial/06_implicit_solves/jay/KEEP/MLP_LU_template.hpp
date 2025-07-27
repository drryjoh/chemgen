#pragma once
#include <iostream>
#include <array>
#include <random>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <algorithm> 
#include <cstddef> 

// template<typename Scalar>
// using activationFunction = void(*)(Scalar&, Scalar, Scalar);


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template<typename Scalar, int output_size, typename ActFun>
inline void Dense_MLP_LU(Scalar* __restrict outputs, const Scalar* __restrict inputs, const Scalar * __restrict weights, const Scalar * __restrict biases, int input_size, ActFun activation_function, Scalar alpha) noexcept 
{
    for(int i = 0; i < output_size; ++i){
        Scalar sum = 0;
        
        for(int j = 0; j < input_size; ++j){
            sum += inputs[j] * weights[j * output_size + i];
        }
        sum += biases[i];
        activation_function(outputs[i], sum, alpha);
    }
}

//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template <typename Scalar = double>
inline auto MLP_LU(const std::array<Scalar, 3256>& initial_input) {

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

    // normalize input
    for (int i = 0; i < 3256; i++) { model_input[i] = (initial_input[i] - input_min_mean[i]) / (input_norm_std[i]); } 

    if (model_input.size() != 3256) { throw std::invalid_argument("Invalid input size. Expected size: 3256"); }

    // Dense, layer 1
    static std::array<Scalar, 1> layer_1_output;
    Dense_MLP_LU<Scalar, 1>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        3256, linear, 0.0);

    // Activation, layer 2
    static std::array<Scalar, 1> layer_2_output;
    for (int i = 0; i < 1; ++i) {
        gelu(layer_2_output[i], layer_1_output[i], 0.0);
    }

    // Dense, layer 3
    static std::array<Scalar, 1> layer_3_output;
    Dense_MLP_LU<Scalar, 1>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        1, linear, 0.0);

    // Activation, layer 4
    static std::array<Scalar, 1> layer_4_output;
    for (int i = 0; i < 1; ++i) {
        gelu(layer_4_output[i], layer_3_output[i], 0.0);
    }

    // Dense, layer 5
    static std::array<Scalar, 9409> layer_5_output;
    Dense_MLP_LU<Scalar, 9409>(
        layer_5_output.data(), layer_4_output.data(),
        weights_5.data(), biases_5.data(),
        1, linear, 0.0);

    // Activation, layer 6
    static std::array<Scalar, 9409> layer_6_output;
    for (int i = 0; i < 9409; ++i) {
        nonzero_diag_activation(layer_6_output[i], layer_5_output[i], i);
    }


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    static std::array<Scalar, 9409> model_output;

    for (int i = 0; i < 9409; i++) { model_output[i] = (layer_6_output[i] * output_norm_std[i]) + output_min_mean[i]; }

    return model_output;

}