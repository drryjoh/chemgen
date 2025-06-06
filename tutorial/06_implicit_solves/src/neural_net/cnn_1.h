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

template <typename Scalar>
__attribute__((hot))
void DepthwiseForsSeparableConv2D(Scalar *__restrict outputs, const Scalar *__restrict inputs, const Scalar *__restrict weights, const Scalar *__restrict biases,
                                  int out_height, int out_width,
                                  int in_channels, int in_height, int in_width,
                                  int kernel_height, int kernel_width, int stride_height, int stride_width,
                                  int padding_height, int padding_width) noexcept
{
    #pragma omp parallel for collapse(2)
    for (int c = 0; c < in_channels; ++c)
    {
        for (int oh = 0; oh < out_height; ++oh)
        {
            for (int ow = 0; ow < out_width; ++ow)
            {
                Scalar sum = 0;
                for (int kh = 0; kh < kernel_height; ++kh)
                {
                    #pragma clang loop vectorize(enable)
                    for (int kw = 0; kw < kernel_width; ++kw)
                    {
                        int in_h = oh * stride_height - padding_height + kh;
                        int in_w = ow * stride_width - padding_width + kw;
                        if (in_h >= 0 && in_h < in_height && in_w >= 0 && in_w < in_width)
                        {
                            int input_index = (in_h * in_width * in_channels) + (in_w * in_channels) + c;
                            int weight_index = (kh * kernel_width + kw) * in_channels + c;
                            sum += inputs[input_index] * weights[weight_index];
                        }
                    }
                }
                sum += biases[c];
                int output_index = ((oh * out_width + ow) * in_channels) + c;
                outputs[output_index] = sum;
            }
        }
    }
}

template <typename Scalar, int out_channels, int out_height, int out_width, typename ActFun>
__attribute__((hot))
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
    std::vector<Scalar> depthwise_output(out_height * out_width * in_channels, 0);
    std::vector<Scalar> zero_bias(in_channels, 0);
    DepthwiseForsSeparableConv2D(
        depthwise_output.data(), inputs, depthwise_weights, zero_bias.data(), out_height, out_width,
        in_channels, in_height, in_width,
        kernel_height, kernel_width,
        stride_height, stride_width,
        padding_height, padding_width);

    #pragma omp parallel for collapse(2)
    for (int oc = 0; oc < out_channels; ++oc)
    {
        for (int i = 0; i < out_height * out_width; ++i)
        {
            Scalar sum = 0;
            #pragma clang loop vectorize(enable)
            for (int ic = 0; ic < in_channels; ++ic)
            {
                int index = i * in_channels + ic;
                int weight_index = ic * out_channels + oc;
                sum += depthwise_output[index] * pointwise_weights[weight_index];
            }
            sum += biases[oc];
            int output_index = i * out_channels + oc;
            activation_function(outputs[output_index], sum, alpha);
        }
    }
}

template <typename Scalar, int N>
__attribute__((hot))
void Reshape(Scalar *__restrict outputs, const Scalar *__restrict inputs) noexcept
{
    #pragma clang loop vectorize(enable)
    for (int i = 0; i < N; ++i)
    {
        outputs[i] = inputs[i];
    }
}

//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//

template <typename Scalar = double>
auto cnn_1(const std::array<std::array<Scalar, 96>, 96> &initial_input)
{

    constexpr int flat_size = 9216;
    std::array<Scalar, flat_size> model_input;
    #pragma omp parallel for collapse(2)
    for (int i0 = 0; i0 < 96; i0++)
    {
        for (int i1 = 0; i1 < 96; i1++)
        {
            int flatIndex = i0 * 96 + i1 * 1;
            model_input[flatIndex] = initial_input[i0][i1];
        }
    }
    if (model_input.size() != 9216)
    {
        throw std::invalid_argument("Invalid input size. Expected size: 9216");
    }

    //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//

    // Layer 1: SeparableConv2D
    constexpr std::array<Scalar, 9> sepDepthwise_1 = {-3.194657266e-01, -3.036856651e-04, -1.974890232e-01, -4.558111429e-01, -1.342993975e-02, 1.759145856e-01, 2.385321856e-01, -2.408222556e-01, -5.319681168e-01};
    constexpr std::array<Scalar, 16> sepPointwise_1 = {-5.669549704e-01, -4.273068309e-01, -4.470078945e-01, -2.155958712e-01, 3.856537938e-01, -5.370312929e-01, -3.147925436e-01, 1.794828773e-01, 1.297131181e-01, -5.101598501e-01, 3.920166492e-01, -4.006438851e-01, 2.864408493e-02, -3.715263009e-01, -5.849651098e-01, 5.528271198e-03};
    constexpr std::array<Scalar, 16> sepPointwiseBias_1 = {4.836391983e-23, -3.281982716e-22, 6.662691438e-23, -5.070745045e-24, 8.159794819e-23, -4.681225820e-22, -7.427730006e-23, -1.689771487e-22, 1.571534743e-22, -3.971002103e-22, -1.985183614e-22, 2.014975167e-22, -1.198653824e-22, 6.742931213e-23, -1.326279276e-22, -8.016485298e-23};

    // Layer 2: SeparableConv2D
    constexpr std::array<Scalar, 144> sepDepthwise_2 = {-1.828941405e-01, 6.074102223e-02, 1.933371872e-01, -2.597969770e-02, 7.282619178e-02, -1.613444388e-01, -9.405256808e-02, 8.745475113e-02, 8.284918964e-02, 2.267056704e-02, -1.709357053e-01, -7.463126630e-02, 1.576061696e-01, -4.400908947e-03, -2.857579291e-02, -5.873835087e-02, 1.098805219e-01, -1.654600203e-01, 4.089291394e-02, 1.573240310e-01, -4.101255536e-02, -1.866163313e-02, -4.173950851e-02, -4.556100070e-02, -3.978040814e-02, 1.784670800e-01, 1.179122478e-01, 7.313089073e-02, 8.014659584e-02, 1.234695315e-03, -1.935534775e-01, -1.961134672e-01, 8.437083662e-02, 1.122142375e-02, -1.026864797e-01, 3.609621525e-02, 2.144399285e-02, 8.648455143e-03, 1.201379746e-01, 4.138506949e-02, 1.906290203e-01, -1.150963306e-01, -1.474794596e-01, -2.167363465e-02, 6.626665592e-03, 7.789351046e-02, 1.533574313e-01, 1.265193075e-01, -5.947305262e-02, -1.042188108e-01, -1.660682261e-01, -1.484819949e-01, -1.468001902e-01, 4.909846187e-02, 9.458924830e-02, 1.841193289e-01, -7.517271489e-02, 1.304035038e-01, 1.847884208e-01, 1.170968264e-01, -3.098058701e-02, -1.772590131e-01, 4.248675704e-02, -7.829166949e-02, -3.335078061e-02, 5.505631864e-02, 1.610680968e-01, 1.464160532e-01, 7.315893471e-02, -1.926940978e-01, -9.565208107e-02, 9.662790596e-02, -4.860086739e-02, 1.119725406e-02, 1.170700043e-01, 1.776761860e-01, -1.720257849e-01, 1.873487681e-01, -7.154479623e-02, 3.033876419e-02, -1.939541399e-01, -1.916230023e-01, -2.495869994e-03, 1.663246602e-01, 2.873386443e-02, -7.222227752e-02, 5.027459562e-02, -1.669072658e-01, -1.291565597e-02, -6.603766978e-02, -2.786126733e-02, 9.814269841e-02, 7.986001670e-02, -1.584332734e-01, 8.983056247e-02, 6.545995176e-02, -3.482927382e-02, 7.653404772e-02, 9.219615161e-02, 1.747177392e-01, -1.394797862e-02, -1.242075264e-01, 1.762511283e-01, -8.932900429e-02, -1.742781252e-01, -8.857420087e-02, -5.421304703e-02, 1.511689872e-01, 1.154673249e-01, 7.219819725e-02, 1.943603158e-02, -1.631765068e-02, -4.684677720e-02, -4.366305470e-02, 1.067570299e-01, 9.916618466e-03, -7.543583959e-02, 3.798741102e-02, -1.582898796e-01, -4.754564166e-02, 7.114674151e-02, 1.598043591e-01, -1.289439797e-01, -4.086939991e-02, 1.621159911e-02, 7.174097002e-02, 1.288489252e-01, 1.210371107e-01, 3.289863467e-03, 1.675782651e-01, -1.329993755e-01, -1.789908707e-01, -1.427309215e-02, -1.745583862e-01, 4.053738713e-02, -2.595651150e-02, -4.005095363e-02, 6.250752509e-02, -3.605528176e-02, 8.090643585e-02, 1.892099530e-01, -1.603445411e-01, -1.521949470e-01, 1.701466590e-01};
    constexpr std::array<Scalar, 128> sepPointwise_2 = {3.058035374e-01, -3.844397068e-01, -3.206810951e-01, 1.922999620e-01, 2.460091114e-01, -2.789487839e-01, -2.750343084e-01, -2.153762579e-01, -4.055135250e-01, -4.913132191e-01, -1.189451218e-01, -3.082761765e-01, 4.704407454e-01, -6.588852406e-02, 4.367129803e-01, -1.349384785e-01, -5.025994778e-02, 3.092460632e-01, 1.350045204e-02, 1.216026545e-01, 1.170319319e-01, -2.823212147e-01, 3.864884377e-03, 2.308752537e-01, -4.959213734e-01, 3.462629318e-01, 2.385395765e-01, 1.402449608e-02, -8.114266396e-02, 2.273044586e-01, 3.158487082e-01, 1.997592449e-01, -1.800526381e-01, -3.414208889e-01, 3.108975887e-01, 1.381292343e-01, -3.685181141e-01, -4.838943481e-03, -2.449201345e-01, 3.929287195e-01, -4.275003672e-01, 4.184910059e-01, 3.733992577e-01, -4.954630136e-01, -2.538058758e-01, 2.868409157e-01, -1.480547190e-01, -5.164182186e-02, 3.103518486e-02, -4.691792727e-01, 1.157697439e-01, -3.708064556e-02, 3.529808521e-01, 2.013683319e-03, -4.878575802e-01, 1.964991093e-01, -3.523126841e-01, -2.281296253e-01, -2.820088863e-01, 4.920297861e-01, -2.388468981e-01, 1.775932312e-01, 2.852872610e-01, 2.049248219e-01, -1.401129961e-01, -4.080679417e-01, 1.777737141e-01, 2.996531725e-01, 1.343181133e-01, -1.567560434e-01, 1.489995718e-01, 1.473931074e-01, -1.332372427e-01, 5.633270741e-02, -4.887466431e-01, 4.767827988e-01, -1.720434427e-01, 3.315891027e-01, 4.767048359e-02, 1.355928183e-01, 4.393696785e-01, -2.406582832e-01, -6.960940361e-02, 3.713665009e-01, 3.432186842e-01, -3.061507940e-01, 2.689907551e-01, -1.721572876e-01, -8.249413967e-02, 3.567901850e-01, 1.650588512e-01, -4.692745209e-01, -4.242643118e-01, 7.983493805e-02, -2.130259275e-01, 3.689675331e-01, 2.329266071e-02, 1.359157562e-01, -1.136033535e-01, -1.020848751e-02, 4.021387100e-01, 3.721988201e-02, -3.037834167e-01, -4.086544514e-01, -1.557934284e-01, 3.202357292e-01, 4.147429466e-01, -4.493849277e-01, 2.445739508e-01, 9.506440163e-02, -4.446213245e-01, -3.484926224e-01, -4.498994350e-02, 3.079851866e-01, 3.652124405e-01, -2.785505056e-01, 1.490961313e-01, 4.402399063e-01, -2.908433676e-01, 2.035136223e-01, 4.402571917e-01, -6.340491772e-02, 1.773126125e-01, 2.107881308e-01, 2.579500675e-01, 3.628349304e-02, 2.982240915e-01, 1.573982239e-01};
    constexpr std::array<Scalar, 8> sepPointwiseBias_2 = {5.318528496e-22, -1.240257364e-21, 4.844014272e-23, -2.687248287e-21, -2.829391127e-22, 7.732962062e-24, -8.457726757e-23, 1.588597615e-22};

    // Layer 3: SeparableConv2D
    constexpr std::array<Scalar, 72> sepDepthwise_3 = {-2.479105592e-01, 2.073075473e-01, -1.569555104e-01, -1.038270593e-01, -1.516292691e-01, 2.287142277e-01, 1.349220276e-01, -1.612072587e-01, -1.596421897e-01, -1.949797571e-01, 2.834209800e-02, 2.362205982e-01, 1.695966721e-01, -6.490184367e-02, 5.759587884e-02, 6.061992049e-02, 8.630293608e-02, -2.217021585e-02, 6.048533320e-02, 2.528211474e-01, 4.221957922e-02, 1.928371191e-01, -2.510090172e-01, 8.371746540e-02, 1.000213027e-01, -9.199190140e-02, 2.490195036e-01, 7.981112599e-02, -8.136315644e-02, -5.886256695e-04, -1.967440248e-01, -1.472487897e-01, -1.688047051e-01, 2.171027958e-01, -9.040483832e-02, 2.649760842e-01, -1.804365665e-01, -1.764585227e-01, 2.453678846e-02, 1.288628578e-01, -7.248890400e-02, 1.016153097e-01, -2.570966184e-01, 2.049084604e-01, -1.688877642e-01, -3.278811276e-02, 1.630732417e-02, 1.127856076e-01, 3.329625726e-02, -1.700693369e-01, 1.009075642e-01, 1.147482693e-01, -1.407119632e-02, -3.558978438e-02, -1.354883164e-01, 4.804134369e-03, 2.242375016e-01, 2.669328451e-01, -3.041492403e-02, 1.038917005e-01, 1.269389689e-01, 2.461024523e-01, 2.594601512e-01, -1.459965557e-01, 6.662672758e-02, -2.645448744e-01, 1.875334978e-02, -1.848568916e-01, 3.343278170e-02, 4.284545779e-02, -1.903593540e-02, 1.606906056e-01};
    constexpr std::array<Scalar, 8> sepPointwise_3 = {5.052667856e-01, -4.886442125e-01, 6.013858318e-01, 5.992780924e-01, -7.385904789e-01, -2.443490028e-01, -5.580663681e-02, 7.488929033e-01};
    constexpr std::array<Scalar, 1> sepPointwiseBias_3 = {1.764519403e-22};

    //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//

    auto linear = +[](Scalar &output, Scalar input, Scalar alpha) noexcept
    {
        output = input;
    };

    auto relu = +[](Scalar &output, Scalar input, Scalar alpha) noexcept
    {
        output = input > 0 ? input : 0;
    };

    //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//

    // SeparableConv2D, layer 1
    std::array<Scalar, (96 * 96 * 16)> layer_1_output;
    SeparableConv2D<Scalar, 16, 96, 96>(
        layer_1_output.data(), model_input.data(),
        sepDepthwise_1.data(), sepPointwise_1.data(), sepPointwiseBias_1.data(),
        1, 96, 96,
        3, 3, 1, 1, 1, 1,
        relu, 0.0);

    // SeparableConv2D, layer 2
    std::array<Scalar, (96 * 96 * 8)> layer_2_output;
    SeparableConv2D<Scalar, 8, 96, 96>(
        layer_2_output.data(), layer_1_output.data(),
        sepDepthwise_2.data(), sepPointwise_2.data(), sepPointwiseBias_2.data(),
        16, 96, 96,
        3, 3, 1, 1, 1, 1,
        relu, 0.0);

    // SeparableConv2D, layer 3
    std::array<Scalar, (96 * 96 * 1)> layer_3_output;
    SeparableConv2D<Scalar, 1, 96, 96>(
        layer_3_output.data(), layer_2_output.data(),
        sepDepthwise_3.data(), sepPointwise_3.data(), sepPointwiseBias_3.data(),
        8, 96, 96,
        3, 3, 1, 1, 1, 1,
        linear, 0.0);

    // SeparableConv2D, layer 4
    std::array<Scalar, 9216> layer_4_output;
    Reshape<Scalar, 9216>(
        layer_4_output.data(), layer_3_output.data());

    // Final output
    std::array<std::array<Scalar, 96>, 96> model_output;
    for (int i = 0; i < 96; i++)
    {
        for (int j = 0; j < 96; j++)
        {
            model_output[i][j] = layer_4_output[i * 96 + j];
        }
    }

    return model_output;
}
