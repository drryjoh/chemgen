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
auto mlp_zigzag_4_relu(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

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
    constexpr std::array<Scalar, 4> weights_1 = {6.957931519e-01, -1.929575205e-01, -8.212412000e-01, -6.948803067e-01};
    constexpr std::array<Scalar, 4> biases_1 = {-1.892003264e-12, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00};

    // Dense layer 2
    constexpr std::array<Scalar, 8> weights_2 = {-3.758630753e-01, 4.819409847e-01, -8.407347202e-01, -3.895025253e-01, -7.700242996e-01, 5.328996181e-01, 5.211508274e-01, -5.647726059e-01};
    constexpr std::array<Scalar, 2> biases_2 = {0.000000000e+00, -3.925798055e-12};

    // Dense layer 3
    constexpr std::array<Scalar, 8> weights_3 = {8.484883308e-01, -5.977959633e-01, 3.989272118e-01, -4.616789818e-01, 3.107583523e-01, -1.188626289e-01, -3.185095787e-01, 4.404437542e-01};
    constexpr std::array<Scalar, 4> biases_3 = {-2.764495056e-11, 0.000000000e+00, 0.000000000e+00, 2.014398170e-14};

    // Dense layer 4
    constexpr std::array<Scalar, 4> weights_4 = {8.107534647e-01, 4.818141460e-02, 2.358092070e-01, -2.624163628e-01};
    constexpr std::array<Scalar, 1> biases_4 = {2.564331343e-12};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        Scalar tmp = input;
        output = tmp;
    };

    auto relu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        Scalar tmp = 0.;
        if (tmp > 0.) tmp = input;
        output = tmp;
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // Reshape, layer 1
    static std::array<Scalar, 4> layer_1_output;
    Dense<Scalar, 4>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        9216, relu, 0.0);

    // Reshape, layer 2
    static std::array<Scalar, 2> layer_2_output;
    Dense<Scalar, 2>(
        layer_2_output.data(), layer_1_output.data(),
        weights_2.data(), biases_2.data(),
        4, relu, 0.0);

    // Reshape, layer 3
    static std::array<Scalar, 4> layer_3_output;
    Dense<Scalar, 4>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        2, relu, 0.0);

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
