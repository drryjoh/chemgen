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
auto mlp_reverse_funnel_4_relu(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

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
    constexpr std::array<Scalar, 2> weights_1 = {5.070616007e-01, -3.809258938e-01};
    constexpr std::array<Scalar, 2> biases_1 = {-1.206519016e-13, -1.570957491e-15};

    // Dense layer 2
    constexpr std::array<Scalar, 8> weights_2 = {-7.410044670e-01, -5.292685032e-01, 6.322827339e-01, 6.360354424e-01, -6.227693558e-01, -5.413010120e-01, -9.173018932e-01, 1.217648983e-01};
    constexpr std::array<Scalar, 4> biases_2 = {0.000000000e+00, 0.000000000e+00, 1.854988433e-14, -2.197437517e-13};

    // Dense layer 3
    constexpr std::array<Scalar, 32> weights_3 = {2.064623833e-01, -6.813353300e-01, -4.773889184e-01, -6.112479568e-01, 5.372260213e-01, 6.477665305e-01, 2.865326405e-02, 3.524845243e-01, 1.863268018e-01, 6.646514535e-01, 2.931807637e-01, 5.539940000e-01, -6.111695766e-01, 6.651338935e-01, 1.265243292e-01, -5.008494854e-02, -2.453386784e-02, 9.981775284e-02, 6.236757636e-01, -4.730135798e-01, -3.025567234e-01, 5.157981515e-01, 3.540926576e-01, -3.576254547e-01, 1.824523807e-01, 3.758251667e-02, 4.725467563e-01, 6.434758306e-01, 5.957369208e-01, -5.777087212e-01, -2.194060385e-01, -4.132058322e-01};
    constexpr std::array<Scalar, 8> biases_3 = {6.399141640e-15, 6.375550078e-15, 4.847229389e-14, 1.830257950e-14, -1.097512872e-15, -5.560599182e-14, 1.184296368e-12, 0.000000000e+00};

    // Dense layer 4
    constexpr std::array<Scalar, 8> weights_4 = {-6.025623679e-01, -5.888378620e-01, 3.610516787e-01, 3.359398842e-01, -4.721916616e-01, -7.588384748e-01, -5.117872953e-01, -8.102403879e-01};
    constexpr std::array<Scalar, 1> biases_4 = {1.342532043e-13};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input;
    };

    auto relu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input > 0 ? input : 0;
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // Reshape, layer 1
    static std::array<Scalar, 2> layer_1_output;
    Dense<Scalar, 2>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        9216, relu, 0.0);

    // Reshape, layer 2
    static std::array<Scalar, 4> layer_2_output;
    Dense<Scalar, 4>(
        layer_2_output.data(), layer_1_output.data(),
        weights_2.data(), biases_2.data(),
        2, relu, 0.0);

    // Reshape, layer 3
    static std::array<Scalar, 8> layer_3_output;
    Dense<Scalar, 8>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        4, relu, 0.0);

    // Reshape, layer 4
    static std::array<Scalar, 1> layer_4_output;
    Dense<Scalar, 1>(
        layer_4_output.data(), layer_3_output.data(),
        weights_4.data(), biases_4.data(),
        8, linear, 0.0);

    // Final output
    static std::array<std::array<Scalar, 96>, 96> model_output;
    for(int i = 0; i < 96; i++) {
        for(int j = 0; j < 96; j++) {
            model_output[i][j] = layer_4_output[i * 96 + j];
        }
    }

    return model_output;
}
