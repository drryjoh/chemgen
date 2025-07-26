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
auto mlp_funnel_4_leakyrelu(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

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
    constexpr std::array<Scalar, 8> weights_1 = {4.791432619e-01, -3.779920340e-01, -4.071167111e-01, -6.327797174e-01, -9.583812952e-02, -7.205038667e-01, -5.699862242e-01, 2.873727083e-01};
    constexpr std::array<Scalar, 8> biases_1 = {3.737172518e-13, -7.559097098e-11, -4.500406992e-11, 1.315032198e-11, 6.426670394e-12, 5.108701309e-11, -3.174542573e-11, 6.115419044e-13};

    // Dense layer 2
    constexpr std::array<Scalar, 32> weights_2 = {-6.060683727e-02, -2.492296398e-01, 1.901065707e-01, 3.044866920e-01, 2.210976481e-01, -5.655177832e-01, -6.954991817e-01, 5.437760949e-01, -2.016283274e-01, -2.971828282e-01, -6.348375678e-01, 4.710010886e-01, 3.621390462e-01, -1.770614386e-01, 4.461417794e-01, -2.638998032e-01, 7.059398293e-01, 1.162777543e-01, 8.230674267e-02, -5.312535763e-01, -2.253878415e-01, -2.381853163e-01, 4.208460450e-01, -3.674092591e-01, 3.606432080e-01, -2.438240647e-01, 4.652912021e-01, 3.291004300e-01, -6.215798855e-02, -1.170448065e-01, 5.281059146e-01, -2.521307170e-01};
    constexpr std::array<Scalar, 4> biases_2 = {-5.231554773e-11, 3.712291325e-12, 2.637159241e-11, -6.967432881e-11};

    // Dense layer 3
    constexpr std::array<Scalar, 8> weights_3 = {6.621413231e-01, -1.506648064e-01, -5.669374466e-01, 5.879852772e-01, -3.482966423e-01, -2.783846855e-02, 8.501431942e-01, 4.535365105e-02};
    constexpr std::array<Scalar, 2> biases_3 = {-7.565900684e-11, -7.116003750e-13};

    // Dense layer 4
    constexpr std::array<Scalar, 2> weights_4 = {7.071503401e-01, 1.374034762e+00};
    constexpr std::array<Scalar, 1> biases_4 = {-1.323620141e-12};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto leakyrelu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input > 0 ? input : alpha * input;
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
        9216, leakyrelu, 0.1);

    // Reshape, layer 2
    static std::array<Scalar, 4> layer_2_output;
    Dense<Scalar, 4>(
        layer_2_output.data(), layer_1_output.data(),
        weights_2.data(), biases_2.data(),
        8, leakyrelu, 0.1);

    // Reshape, layer 3
    static std::array<Scalar, 2> layer_3_output;
    Dense<Scalar, 2>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        4, leakyrelu, 0.1);

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
