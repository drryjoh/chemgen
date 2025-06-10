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
auto mlp_reverse_funnel_4_leakyrelu(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

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
    constexpr std::array<Scalar, 2> weights_1 = {-1.227754593e+00, -5.355906487e-01};
    constexpr std::array<Scalar, 2> biases_1 = {-5.592578067e-10, 1.886276829e-10};

    // Dense layer 2
    constexpr std::array<Scalar, 8> weights_2 = {-7.662572861e-01, 6.441817284e-01, -5.046100616e-01, 7.503793240e-01, -7.177255154e-01, -5.204215050e-01, -4.627714157e-01, -1.049389839e-01};
    constexpr std::array<Scalar, 4> biases_2 = {8.350065639e-10, -3.365774326e-10, 3.439332152e-10, -4.146113453e-10};

    // Dense layer 3
    constexpr std::array<Scalar, 32> weights_3 = {5.803613067e-01, 2.953531146e-01, -1.568781734e-01, 5.663098693e-01, 8.138102293e-02, 6.921345592e-01, 5.415384173e-01, 3.804787993e-01, -6.116075516e-01, 4.076409340e-02, 3.153359294e-01, 5.366163850e-01, -4.414782822e-01, 2.649651170e-01, -6.568056941e-01, -3.187835813e-01, 2.389287949e-01, 1.631993651e-01, -4.697889090e-02, -2.737228572e-01, -2.337543368e-01, -1.098692417e-03, -4.880079031e-01, 5.414258838e-01, -6.170127988e-01, 1.994828582e-01, 5.929701924e-01, 9.833890200e-02, -6.440626383e-01, -2.055050135e-01, -3.579226732e-01, -2.299014330e-01};
    constexpr std::array<Scalar, 8> biases_3 = {8.754606068e-11, -1.538119805e-11, -5.096980060e-10, -1.090088436e-11, 3.033679974e-10, 4.930426623e-10, -2.055307313e-10, 5.073695353e-10};

    // Dense layer 4
    constexpr std::array<Scalar, 8> weights_4 = {-1.252528429e-01, -5.262034535e-01, -6.961978674e-01, -3.928350806e-01, -7.678035498e-01, -7.241420746e-01, 2.886551619e-01, -7.355484366e-01};
    constexpr std::array<Scalar, 1> biases_4 = {3.466904541e-11};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto leakyrelu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input > 0 ? input : alpha * input;
    };

    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input;
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // Reshape, layer 1
    static std::array<Scalar, 2> layer_1_output;
    Dense<Scalar, 2>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        9216, leakyrelu, 0.1);

    // Reshape, layer 2
    static std::array<Scalar, 4> layer_2_output;
    Dense<Scalar, 4>(
        layer_2_output.data(), layer_1_output.data(),
        weights_2.data(), biases_2.data(),
        2, leakyrelu, 0.1);

    // Reshape, layer 3
    static std::array<Scalar, 8> layer_3_output;
    Dense<Scalar, 8>(
        layer_3_output.data(), layer_2_output.data(),
        weights_3.data(), biases_3.data(),
        4, leakyrelu, 0.1);

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
