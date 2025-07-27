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
auto mlp_funnel_4_relu(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

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
    constexpr std::array<Scalar, 8> weights_1 = {4.581825733e-01, 1.147647500e-01, 2.543224096e-01, -5.062156916e-01, -2.333761454e-01, 3.467770815e-01, 6.691775322e-01, 3.592488766e-01};
    constexpr std::array<Scalar, 8> biases_1 = {0.000000000e+00, 0.000000000e+00, 0.000000000e+00, -1.422822693e-10, -9.067024909e-11, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00};

    // Dense layer 2
    constexpr std::array<Scalar, 32> weights_2 = {-1.108323932e-01, 2.073547244e-01, 6.473726630e-01, -5.648123026e-01, -3.118313849e-01, 6.487177014e-01, 5.751097798e-01, -6.788904667e-01, -3.965077698e-01, 1.365320086e-01, 6.770048738e-01, -4.617431760e-01, 6.596452594e-01, -3.759397566e-01, -1.526710391e-01, -1.269068718e-01, 3.460392356e-01, -4.784220159e-01, 4.841033816e-01, 4.719334245e-01, -3.924274445e-02, 4.610794783e-02, 1.914350390e-01, -4.267891049e-01, -1.051690578e-01, -6.944962740e-01, 5.139730573e-01, 7.065581679e-01, -5.332332850e-02, -6.555446386e-01, 3.001442552e-01, -4.096553922e-01};
    constexpr std::array<Scalar, 4> biases_2 = {-2.171427205e-10, 0.000000000e+00, 1.432230307e-10, -1.798242666e-10};

    // Dense layer 3
    constexpr std::array<Scalar, 8> weights_3 = {4.368276596e-01, 5.032410622e-01, -4.164409637e-02, -2.390978336e-01, -5.462944508e-01, -2.447993755e-01, -2.546441555e-01, 6.132521629e-01};
    constexpr std::array<Scalar, 2> biases_3 = {-1.126263180e-10, -3.337257415e-10};

    // Dense layer 4
    constexpr std::array<Scalar, 2> weights_4 = {3.770524263e-01, 1.089878678e+00};
    constexpr std::array<Scalar, 1> biases_4 = {-7.370681322e-12};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input;
    };

    auto relu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input > 0 ? input : 0;
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // Reshape, layer 1
    static std::array<Scalar, 8> layer_1_output;
    Dense<Scalar, 8>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        9216, relu, 0.0);

    // Reshape, layer 2
    static std::array<Scalar, 4> layer_2_output;
    Dense<Scalar, 4>(
        layer_2_output.data(), layer_1_output.data(),
        weights_2.data(), biases_2.data(),
        8, relu, 0.0);

    // Reshape, layer 3
    static std::array<Scalar, 2> layer_3_output;
    Dense<Scalar, 2>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        4, relu, 0.0);

    // Reshape, layer 4
    static std::array<Scalar, 1> layer_4_output;
    Dense<Scalar, 1>(
        layer_4_output.data(), layer_3_output.data(),
        weights_4.data(), biases_4.data(),
        2, linear, 0.0);

    // Final output
    static std::array<std::array<Scalar, 96>, 96> model_output;
    for(int i = 0; i < 96; i++) {
        for(int j = 0; j < 96; j++) {
            model_output[i][j] = layer_4_output[i * 96 + j];
        }
    }

    return model_output;
}
