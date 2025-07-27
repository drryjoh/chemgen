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
inline void Dense(Scalar* __restrict outputs, const Scalar* __restrict inputs, const Scalar * __restrict weights, const Scalar * __restrict biases, int input_size, ActFun activation_function, Scalar alpha) noexcept {
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
auto mlp_zigzag_4_leakyrelu(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

    constexpr int flat_size = 9216; 
    std::array<Scalar, flat_size> model_input;
    for (int i0 = 0; i0 < 96; i0++) {
      for (int i1 = 0; i1 < 96; i1++) {
            int flatIndex = i0 * 96 + i1 * 1;
            model_input[flatIndex] = initial_input[i0][i1];
        }
    }
    if (model_input.size() != 9216) { throw std::invalid_argument("Invalid input size. Expected size: 9216"); }


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // Dense layer 1
    constexpr std::array<Scalar, 4> weights_1 = {-1.009672046e+00, -7.186513543e-01, 7.725011110e-01, -5.321492553e-01};
    constexpr std::array<Scalar, 4> biases_1 = {-1.060881452e-10, -1.524557425e-10, 2.793653503e-11, 3.147091962e-11};

    // Dense layer 2
    constexpr std::array<Scalar, 8> weights_2 = {-1.827411652e-01, 7.571220398e-01, 1.405467987e-01, 6.466276646e-01, -8.225846291e-01, 1.982920170e-01, 4.706001282e-01, -6.819891930e-01};
    constexpr std::array<Scalar, 2> biases_2 = {-2.118784037e-10, -1.372242736e-10};

    // Dense layer 3
    constexpr std::array<Scalar, 8> weights_3 = {8.605480194e-02, 9.250085354e-01, -3.097934723e-01, 3.228013515e-01, -7.704234123e-02, 4.965837002e-01, 2.882075310e-01, -7.727777958e-01};
    constexpr std::array<Scalar, 4> biases_3 = {3.388196407e-12, -2.629516049e-10, -2.916725672e-12, -4.753422482e-12};

    // Dense layer 4
    constexpr std::array<Scalar, 4> weights_4 = {3.729720116e-01, -4.945461750e-01, -3.505933285e-01, -2.743428946e-01};
    constexpr std::array<Scalar, 1> biases_4 = {2.277270386e-11};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto leakyrelu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input > 0 ? input : alpha * input;
    };

    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input;
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // Reshape, layer 1
    static std::array<Scalar, 4> layer_1_output;
    Dense<Scalar, 4>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        9216, leakyrelu, 0.1);

    // Reshape, layer 2
    static std::array<Scalar, 2> layer_2_output;
    Dense<Scalar, 2>(
        layer_2_output.data(), layer_1_output.data(),
        weights_2.data(), biases_2.data(),
        4, leakyrelu, 0.1);

    // Reshape, layer 3
    static std::array<Scalar, 4> layer_3_output;
    Dense<Scalar, 4>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        2, leakyrelu, 0.1);

    // Reshape, layer 4
    static std::array<Scalar, 1> layer_4_output;
    Dense<Scalar, 1>(
        layer_4_output.data(), layer_3_output.data(),
        weights_4.data(), biases_4.data(),
        4, linear, 0.0);

    // Final output
    static std::array<std::array<Scalar, 96>, 96> model_output;
    for(int i = 0; i < 96; i++) {
        for(int j = 0; j < 96; j++) {
            model_output[i][j] = layer_4_output[i * 96 + j];
        }
    }

    return model_output;
}
