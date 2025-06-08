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
void Dense(Scalar* __restrict outputs, const Scalar* __restrict inputs, const Scalar * __restrict weights, const Scalar * __restrict biases, int input_size, ActFun activation_function, Scalar alpha) noexcept {
    for(int i = 0; i < output_size; ++i){
        Scalar sum = 0;
        // #pragma unroll
        for(int j = 0; j < input_size; ++j){
            sum += inputs[j] * weights[j * output_size + i];
        }
        sum += biases[i];
        activation_function(outputs[i], sum, alpha);
    }
}

template<typename Scalar, int N>
void Reshape(Scalar * __restrict outputs, const Scalar * __restrict inputs) noexcept {
    // #pragma unroll
    for (int i = 0; i < N; ++i) {
        outputs[i] = inputs[i];
    }
}   

//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template <typename Scalar = double>
auto mlp_all(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

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
    constexpr std::array<Scalar, 8> weights_1 = {7.980300188e-01, -2.156733871e-01, -7.709990144e-01, -3.354940712e-01, 8.094003201e-01, -5.145971775e-01, 5.809540749e-01, 5.354760885e-01};
    constexpr std::array<Scalar, 8> biases_1 = {-5.688558305e-11, -3.555638839e-10, -2.385202036e-10, -1.518967591e-10, 8.322533634e-11, -5.707728387e-11, -8.950304387e-11, 9.897933168e-11};

    // Dense layer 2
    constexpr std::array<Scalar, 32> weights_2 = {5.689938664e-01, -5.217683315e-01, -3.381918669e-01, -6.110509038e-01, -7.100427151e-02, 3.391248584e-01, 4.269918799e-01, 5.334430337e-01, 6.967098117e-01, 1.909294128e-01, 4.356622100e-01, 3.399886489e-01, 3.173649311e-03, 3.980503678e-01, 6.048424840e-01, -2.168468833e-01, -4.400248826e-01, -2.241908908e-01, 5.009821057e-01, 6.661260724e-01, 4.335193038e-01, -1.101539731e-01, 5.663999915e-01, 1.991868019e-03, -2.109974623e-02, 7.060433030e-01, -5.688214302e-01, -6.780164838e-01, -6.951868534e-02, -3.152322769e-01, 6.143129468e-01, 3.821434379e-01};
    constexpr std::array<Scalar, 4> biases_2 = {-1.819190144e-13, -3.113339031e-10, -5.079194444e-12, -3.302975116e-10};

    // Dense layer 3
    constexpr std::array<Scalar, 4> weights_3 = {-1.115548611e-02, 4.986207485e-01, 2.756124735e-01, 5.288119316e-01};
    constexpr std::array<Scalar, 1> biases_3 = {-1.842875332e-11};


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

    // Reshape, layer 4
    static std::array<Scalar, 9216> layer_4_output;
    Reshape<Scalar, 9216>(
        layer_4_output.data(), layer_3_output.data());

    // Final output
    static std::array<std::array<Scalar, 96>, 96> model_output;
    for(int i = 0; i < 96; i++) {
        for(int j = 0; j < 96; j++) {
            model_output[i][j] = layer_4_output[i * 96 + j];
        }
    }

    return model_output;
}
