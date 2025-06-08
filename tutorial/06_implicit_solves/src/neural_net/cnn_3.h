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


template<typename Scalar, int N>
void Reshape(Scalar * __restrict outputs, const Scalar * __restrict inputs) noexcept {
    #pragma unroll
    for (int i = 0; i < N; ++i) {
        outputs[i] = inputs[i];
    }
}   

// template <typename Scalar>
// void DepthwiseForsSeparableConv2D(Scalar *__restrict outputs, const Scalar *__restrict inputs, const Scalar *__restrict weights, const Scalar *__restrict biases,
//                                   int out_height, int out_width,
//                                   int in_channels, int in_height, int in_width,
//                                   int kernel_height, int kernel_width, int stride_height, int stride_width,
//                                   int padding_height, int padding_width) noexcept
// {
//     for (int c = 0; c < in_channels; ++c)
//     {
//         for (int oh = 0; oh < out_height; ++oh)
//         {
//             for (int ow = 0; ow < out_width; ++ow)
//             {
//                 Scalar sum = 0;
//                 for (int kh = 0; kh < kernel_height; ++kh)
//                 {
//                     #pragma unroll
//                     for (int kw = 0; kw < kernel_width; ++kw)
//                     {
//                         int in_h = oh * stride_height - padding_height + kh;
//                         int in_w = ow * stride_width - padding_width + kw;
//                         if (in_h >= 0 && in_h < in_height && in_w >= 0 && in_w < in_width)
//                         {
//                             int input_index = (in_h * in_width * in_channels) + (in_w * in_channels) + c;
//                             int weight_index = (kh * kernel_width + kw) * in_channels + c;
//                             sum += inputs[input_index] * weights[weight_index];
//                         }
//                     }
//                 }
//                 sum += biases[c];
//                 int output_index = ((oh * out_width + ow) * in_channels) + c;
//                 outputs[output_index] = sum;
//             }
//         }
//     }
// }
template <typename Scalar>
void DepthwiseForsSeparableConv2D(Scalar *__restrict outputs, const Scalar *__restrict inputs, 
                                 const Scalar *__restrict weights, const Scalar *__restrict biases,
                                 int out_height, int out_width,
                                 int in_channels, int in_height, int in_width,
                                 int kernel_height, int kernel_width, int stride_height, int stride_width,
                                 int padding_height, int padding_width) noexcept
{
    // Pre-calculate strides
    const int in_stride = in_width * in_channels;
    const int out_stride = out_width * in_channels;
    
    #pragma omp parallel for collapse(3)
    for (int c = 0; c < in_channels; ++c) {
        for (int oh = 0; oh < out_height; ++oh) {
            for (int ow = 0; ow < out_width; ++ow) {
                Scalar sum = 0;
                const int base_h = oh * stride_height - padding_height;
                const int base_w = ow * stride_width - padding_width;
                
                // Kernel loop optimization
                for (int kh = 0; kh < kernel_height; ++kh) {
                    const int in_h = base_h + kh;
                    if (in_h >= 0 && in_h < in_height) {
                        const int in_h_offset = in_h * in_stride;
                        const int kh_offset = kh * kernel_width * in_channels;
                        
                        #pragma unroll
                        for (int kw = 0; kw < kernel_width; ++kw) {
                            const int in_w = base_w + kw;
                            if (in_w >= 0 && in_w < in_width) {
                                sum += inputs[in_h_offset + in_w * in_channels + c] * 
                                      weights[kh_offset + kw * in_channels + c];
                            }
                        }
                    }
                }
                outputs[oh * out_stride + ow * in_channels + c] = sum + biases[c];
            }
        }
    }
}

// template <typename Scalar, int out_channels, int out_height, int out_width, typename ActFun>
// void SeparableConv2D(
//     Scalar *__restrict outputs,
//     const Scalar *__restrict inputs,
//     const Scalar *__restrict depthwise_weights,
//     const Scalar *__restrict pointwise_weights,
//     const Scalar *__restrict biases,
//     int in_channels, int in_height, int in_width,
//     int kernel_height, int kernel_width,
//     int stride_height, int stride_width,
//     int padding_height, int padding_width,
//     ActFun activation_function, Scalar alpha) noexcept
// {
//     std::vector<Scalar> depthwise_output(out_height * out_width * in_channels, 0);
//     std::vector<Scalar> zero_bias(in_channels, 0);
//     DepthwiseForsSeparableConv2D(
//         depthwise_output.data(), inputs, depthwise_weights, zero_bias.data(), out_height, out_width,
//         in_channels, in_height, in_width,
//         kernel_height, kernel_width,
//         stride_height, stride_width,
//         padding_height, padding_width);
//     for (int oc = 0; oc < out_channels; ++oc)
//     {
//         for (int i = 0; i < out_height * out_width; ++i)
//         {
//             Scalar sum = 0;
//             #pragma unroll
//             for (int ic = 0; ic < in_channels; ++ic)
//             {
//                 int index = i * in_channels + ic;
//                 int weight_index = ic * out_channels + oc;
//                 sum += depthwise_output[index] * pointwise_weights[weight_index];
//             }
//             sum += biases[oc];
//             int output_index = i * out_channels + oc;
//             activation_function(outputs[output_index], sum, alpha);
//         }
//     }
// }
// Replace the problematic SeparableConv2D implementation with this fixed version:
template <typename Scalar, int out_channels, int out_height, int out_width, typename ActFun>
void SeparableConv2D(
    Scalar *__restrict outputs,
    const Scalar *__restrict inputs,
    const Scalar *__restrict depthwise_weights,
    const Scalar *__restrict pointwise_weights,
    const Scalar *__restrict biases,
    int in_channels, int in_height, int in_width,
    int kernel_height, int kernel_width,
    int stride_height, int stride_width,
    int padding_height, int padding_width,
    ActFun activation_function, Scalar alpha) noexcept
{
    // Use dynamic allocation instead of VLA or constexpr
    std::vector<Scalar> depthwise_output(out_height * out_width * in_channels);
    std::vector<Scalar> zero_bias(in_channels, 0);

    // Process depthwise convolution
    DepthwiseForsSeparableConv2D(
        depthwise_output.data(), inputs, depthwise_weights, zero_bias.data(),
        out_height, out_width,
        in_channels, in_height, in_width,
        kernel_height, kernel_width,
        stride_height, stride_width,
        padding_height, padding_width);

    // Process pointwise convolution with tiling
    constexpr int TILE_SIZE = 32;
    
    #pragma omp parallel for collapse(2)
    for (int oh = 0; oh < out_height; oh += TILE_SIZE) {
        for (int ow = 0; ow < out_width; ow += TILE_SIZE) {
            const int tile_h = std::min(TILE_SIZE, out_height - oh);
            const int tile_w = std::min(TILE_SIZE, out_width - ow);
            
            // Process tiles
            for (int i = 0; i < tile_h; ++i) {
                for (int j = 0; j < tile_w; ++j) {
                    for (int oc = 0; oc < out_channels; ++oc) {
                        Scalar sum = 0;
                        const int base_idx = ((oh + i) * out_width + (ow + j)) * in_channels;
                        
                        #pragma unroll
                        for (int ic = 0; ic < in_channels; ++ic) {
                            sum += depthwise_output[base_idx + ic] * 
                                  pointwise_weights[ic * out_channels + oc];
                        }
                        
                        const int out_idx = ((oh + i) * out_width + (ow + j)) * out_channels + oc;
                        activation_function(outputs[out_idx], sum + biases[oc], alpha);
                    }
                }
            }
        }
    }
}


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template <typename Scalar = double>
auto cnn_3(const std::array<std::array<Scalar, 96>, 96>& initial_input) {

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


    // Layer 1: SeparableConv2D
    constexpr std::array<Scalar, 9> sepDepthwise_1 = {-2.662068605e-01, -3.992708623e-01, 3.533266783e-01, -2.172652818e-02, -5.384864807e-01, 3.982106745e-01, -4.763884544e-01, -2.648780048e-01, -3.948049843e-01};
    constexpr std::array<Scalar, 8> sepPointwise_1 = {2.596438229e-01, -7.146233320e-02, -5.260425806e-01, 1.440117508e-01, -1.602113829e-03, -8.784379065e-02, -4.713116288e-01, -3.961979151e-01};
    constexpr std::array<Scalar, 8> sepPointwiseBias_1 = {1.082712710e-01, -6.966398563e-03, 4.841496237e-03, -9.824381769e-02, -3.581824154e-02, 1.053147614e-01, 9.868926555e-02, 1.382113397e-01};

    // Layer 2: SeparableConv2D
    constexpr std::array<Scalar, 72> sepDepthwise_2 = {-1.691844463e-01, -1.150052249e-01, -1.945297420e-01, -1.032571271e-01, -1.217945293e-01, 1.981890500e-01, 8.515224606e-02, -2.259479612e-01, 2.969646268e-02, -4.509807378e-02, -1.006831750e-01, 9.832426906e-02, -2.678280510e-02, 2.382079214e-01, 1.986369342e-01, -2.925040424e-01, 9.041867405e-02, -7.917849720e-02, 7.069707662e-02, -2.233926207e-01, 9.597239643e-02, 6.091910228e-02, -7.717417181e-02, 1.568229944e-01, 1.057828963e-01, -1.118271425e-01, -2.478459328e-01, 1.676519364e-01, -5.265040696e-02, 1.547919214e-01, 2.623814940e-01, -8.592943847e-02, 6.613541394e-03, 1.847235411e-01, 1.045466065e-01, -1.618703753e-01, 2.462237328e-01, -1.286749542e-01, 3.425717819e-03, -6.354456395e-02, -2.066545337e-01, 2.847931087e-01, 1.093171015e-01, -1.804537475e-01, -2.286223322e-01, 1.875100434e-01, 3.046096265e-01, -8.686222881e-02, 2.906547673e-02, -2.615128458e-02, 4.804091901e-02, -2.438906953e-02, 2.187240124e-01, -5.347016267e-03, 1.884182841e-01, -1.023382172e-01, -1.637305617e-01, -2.382215261e-01, 9.406100214e-02, -1.719868183e-01, -6.886783987e-02, 3.731203033e-03, -2.149452083e-02, 1.958181709e-01, -1.508332938e-01, 2.302172594e-02, 8.988729119e-02, 1.224084198e-01, -4.493437707e-03, 3.027247190e-01, 2.251225412e-01, -6.424529850e-02};
    constexpr std::array<Scalar, 32> sepPointwise_2 = {4.126499593e-01, 4.727780521e-01, 7.008569241e-01, -5.160987377e-01, -2.012285590e-01, 4.255146720e-03, 7.953302562e-02, 2.430222631e-01, 1.294487268e-01, -1.229873076e-01, 4.912725091e-01, -6.647866368e-01, -6.360049248e-01, -5.391793251e-01, 1.526619643e-01, 7.263239622e-01, 1.686060429e-01, 1.000533551e-01, 2.452695668e-01, 6.689151525e-01, 5.101766586e-01, -2.025424540e-01, -5.015365481e-01, -6.801244020e-01, -5.730250478e-02, -1.696929634e-01, -2.336276919e-01, -5.202180743e-01, -9.609027207e-02, 2.176066935e-01, 2.314888984e-01, -1.023961045e-02};
    constexpr std::array<Scalar, 4> sepPointwiseBias_2 = {-9.114835411e-02, -7.669080794e-02, -1.037833020e-01, -6.347419322e-02};

    // Layer 3: SeparableConv2D
    constexpr std::array<Scalar, 4> sepDepthwise_3 = {2.016551197e-01, 8.109325767e-01, 2.098383904e-01, 7.711197734e-01};
    constexpr std::array<Scalar, 4> sepPointwise_3 = {4.865278304e-01, 8.910318017e-01, 9.092211723e-01, 4.985754192e-01};
    constexpr std::array<Scalar, 1> sepPointwiseBias_3 = {-9.440560639e-02};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input;
    };

    auto relu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept {
        output = input > 0 ? input : 0;
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    // SeparableConv2D, layer 1
    static std::array<Scalar, (96 * 96 * 8)> layer_1_output;
    SeparableConv2D<Scalar, 8, 96, 96>(
        layer_1_output.data(), model_input.data(),
        sepDepthwise_1.data(), sepPointwise_1.data(), sepPointwiseBias_1.data(),
        1, 96, 96,
        3, 3, 1, 1, 1, 1,
        relu, 0.0);

    // SeparableConv2D, layer 2
    static std::array<Scalar, (96 * 96 * 4)> layer_2_output;
    SeparableConv2D<Scalar, 4, 96, 96>(
        layer_2_output.data(), layer_1_output.data(),
        sepDepthwise_2.data(), sepPointwise_2.data(), sepPointwiseBias_2.data(),
        8, 96, 96,
        3, 3, 1, 1, 1, 1,
        relu, 0.0);

    // SeparableConv2D, layer 3
    static std::array<Scalar, (96 * 96 * 1)> layer_3_output;
    SeparableConv2D<Scalar, 1, 96, 96>(
        layer_3_output.data(), layer_2_output.data(),
        sepDepthwise_3.data(), sepPointwise_3.data(), sepPointwiseBias_3.data(),
        4, 96, 96,
        1, 1, 1, 1, 0, 0,
        linear, 0.0);

    // SeparableConv2D, layer 4
    static std::array<Scalar, 9216> layer_4_output;
    Reshape<Scalar, 9216>(
        layer_4_output.data(), layer_3_output.data());

    // Final output
    static std::array< std::array<Scalar, 96>, 96> model_output;
    for(int i = 0; i < 96; i++) {
        for(int j = 0; j < 96; j++) {
            model_output[i][j] = layer_4_output[i * 96 + j];
        }
    }

    return model_output;
}
