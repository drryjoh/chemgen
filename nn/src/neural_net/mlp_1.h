#pragma once
#include <iostream>
#include <array>
#include <random>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <algorithm> 
#include <cstddef> 


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template<typename Scalar, int output_size, typename ActFun>
// __attribute__((hot))
inline void Dense(Scalar* __restrict outputs, const Scalar* __restrict inputs, const Scalar * __restrict weights, const Scalar * __restrict biases, int input_size, ActFun activation_function, Scalar alpha) noexcept {
    // #pragma clang loop vectorize(enable)
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
auto mlp_1ds_3hl(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

    constexpr int flat_size = 9216; 
    std::array<Scalar, flat_size> model_input;
    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < 96; ++i) {
        for (int j = 0; j < 96; ++j) {
            model_input[i * 96 + j] = initial_input[i][j];
        }
    }
    if (model_input.size() != 9216) { throw std::invalid_argument("Invalid input size. Expected size: 9216"); }


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // Dense layer 1
    constexpr std::array<Scalar, 8> weights_1 = {5.826035142e-01, -5.843822956e-01, 7.196350694e-01, -9.835799038e-02, 1.361447573e-01, -4.910477102e-01, 5.014869571e-01, 2.078979611e-01};
    constexpr std::array<Scalar, 8> biases_1 = {7.625377923e-02, 4.501570016e-02, -2.721223794e-02, 4.279249907e-02, -7.704723626e-02, -3.496622294e-02, -6.103385240e-02, 7.810089737e-02};

    // Dense layer 2
    constexpr std::array<Scalar, 32> weights_2 = {5.093327165e-01, -7.448778749e-01, -5.386130810e-01, 8.573333919e-02, -2.039985210e-01, 5.436694026e-01, -5.871812701e-01, -2.313625515e-01, 2.533885837e-01, 4.147711992e-01, -4.128219783e-01, 3.149982989e-01, -4.475685954e-01, 1.457353383e-01, -3.727089465e-01, 4.530324042e-01, 2.475436032e-01, 2.065295130e-01, -3.352370262e-01, -6.049255133e-01, 3.324220479e-01, -6.761999726e-01, 5.370734334e-01, -1.645411849e-01, 1.833023578e-01, 6.276357770e-01, -2.989415228e-01, -1.060793027e-01, -6.232430339e-01, 5.162175894e-01, -4.950504899e-01, 6.595013738e-01};
    constexpr std::array<Scalar, 4> biases_2 = {-5.624772981e-02, -5.844761804e-02, 0.000000000e+00, 9.023229033e-02};

    // Dense layer 3
    constexpr std::array<Scalar, 4> weights_3 = {-9.214295447e-02, -1.323880255e-01, 1.103854179e-02, 3.071466386e-01};
    constexpr std::array<Scalar, 1> biases_3 = {6.458479911e-02};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto relu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input > 0 ? input : 0;
    };

    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input;
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
    static std::array<Scalar, 1> layer_3_output;
    Dense<Scalar, 1>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        4, linear, 0.0);

    // Final output
    static std::array<std::array<Scalar, 96>, 96> model_output;
    // #pragma omp parallel for collapse(2)
    for(int i = 0; i < 96; i++) {
        for(int j = 0; j < 96; j++) {
            model_output[i][j] = layer_3_output[i * 96 + j];
        }
    }

    return model_output;
}