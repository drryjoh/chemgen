#pragma once
#include <iostream>
#include <array>
#include <random>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <algorithm> 
#include <cstddef>
#include <arm_neon.h> 

// template<typename Scalar>
// using activationFunction = void(*)(Scalar&, Scalar, Scalar);


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template <typename Scalar, int size>
inline void UnitNormalization_MLP_LU(Scalar * __restrict outputs,
                              const Scalar * __restrict inputs,
                              Scalar epsilon) noexcept
{
    if constexpr (std::is_same_v<Scalar, float>) {
        // NEON optimized version for float32
        float32x4_t sum_sq_vec = vdupq_n_f32(0.0f);
        int i = 0;
        
        // Vectorized sum of squares
        for(; i <= size - 4; i += 4) {
            float32x4_t input_vec = vld1q_f32(&inputs[i]);
            sum_sq_vec = vfmaq_f32(sum_sq_vec, input_vec, input_vec);
        }
        
        float sum_sq = vaddvq_f32(sum_sq_vec);
        
        // Handle remaining elements
        for(; i < size; ++i) {
            sum_sq += inputs[i] * inputs[i];
        }
        
        float inv_norm = 1.0f / std::sqrt(sum_sq + epsilon);
        float32x4_t inv_norm_vec = vdupq_n_f32(inv_norm);
        
        // Vectorized normalization
        i = 0;
        for(; i <= size - 4; i += 4) {
            float32x4_t input_vec = vld1q_f32(&inputs[i]);
            float32x4_t result_vec = vmulq_f32(input_vec, inv_norm_vec);
            vst1q_f32(&outputs[i], result_vec);
        }
        
        // Handle remaining elements
        for(; i < size; ++i) {
            outputs[i] = inputs[i] * inv_norm;
        }
    } else if constexpr (std::is_same_v<Scalar, double>) {
        // NEON optimized version for float64
        float64x2_t sum_sq_vec = vdupq_n_f64(0.0);
        int i = 0;
        
        // Vectorized sum of squares
        for(; i <= size - 2; i += 2) {
            float64x2_t input_vec = vld1q_f64(&inputs[i]);
            sum_sq_vec = vfmaq_f64(sum_sq_vec, input_vec, input_vec);
        }
        
        double sum_sq = vaddvq_f64(sum_sq_vec);
        
        // Handle remaining elements
        for(; i < size; ++i) {
            sum_sq += inputs[i] * inputs[i];
        }
        
        double inv_norm = 1.0 / std::sqrt(sum_sq + epsilon);
        float64x2_t inv_norm_vec = vdupq_n_f64(inv_norm);
        
        // Vectorized normalization
        i = 0;
        for(; i <= size - 2; i += 2) {
            float64x2_t input_vec = vld1q_f64(&inputs[i]);
            float64x2_t result_vec = vmulq_f64(input_vec, inv_norm_vec);
            vst1q_f64(&outputs[i], result_vec);
        }
        
        // Handle remaining elements
        for(; i < size; ++i) {
            outputs[i] = inputs[i] * inv_norm;
        }
    } else {
        // Fallback with compiler vectorization hints
        Scalar sum_sq = 0;
        
        #pragma clang loop vectorize(enable) unroll(enable)
        for (int i = 0; i < size; ++i) 
        {
            sum_sq += inputs[i] * inputs[i];
        }
        
        Scalar inv_norm = Scalar(1) / std::sqrt(sum_sq + epsilon);
        
        #pragma clang loop vectorize(enable) unroll(enable)
        for (int i = 0; i < size; ++i) 
        {
            outputs[i] = inputs[i] * inv_norm;
        }
    }
}

template<typename Scalar, int output_size, typename ActFun>
inline void Dense_MLP_LU(Scalar* __restrict outputs, const Scalar* __restrict inputs, const Scalar * __restrict weights, const Scalar * __restrict biases, int input_size, ActFun activation_function, Scalar alpha) noexcept 
{
    // Vectorized implementation for better performance on ARM NEON
    if constexpr (std::is_same_v<Scalar, float>) {
        // NEON optimized version for float32
        for(int i = 0; i < output_size; ++i){
            float32x4_t sum_vec = vdupq_n_f32(0.0f);
            int j = 0;
            
            // Process 4 elements at a time
            for(; j <= input_size - 4; j += 4) {
                float32x4_t input_vec = vld1q_f32(&inputs[j]);
                float32x4_t weight_vec = vld1q_f32(&weights[j * output_size + i]);
                sum_vec = vfmaq_f32(sum_vec, input_vec, weight_vec);
            }
            
            // Sum the vector elements
            float sum = vaddvq_f32(sum_vec);
            
            // Handle remaining elements
            for(; j < input_size; ++j) {
                sum += inputs[j] * weights[j * output_size + i];
            }
            
            sum += biases[i];
            activation_function(outputs[i], sum, alpha);
        }
    } else if constexpr (std::is_same_v<Scalar, double>) {
        // NEON optimized version for float64
        for(int i = 0; i < output_size; ++i){
            float64x2_t sum_vec = vdupq_n_f64(0.0);
            int j = 0;
            
            // Process 2 elements at a time for double precision
            for(; j <= input_size - 2; j += 2) {
                float64x2_t input_vec = vld1q_f64(&inputs[j]);
                float64x2_t weight_vec = vld1q_f64(&weights[j * output_size + i]);
                sum_vec = vfmaq_f64(sum_vec, input_vec, weight_vec);
            }
            
            // Sum the vector elements
            double sum = vaddvq_f64(sum_vec);
            
            // Handle remaining elements
            for(; j < input_size; ++j) {
                sum += inputs[j] * weights[j * output_size + i];
            }
            
            sum += biases[i];
            activation_function(outputs[i], sum, alpha);
        }
    } else {
        // Fallback for other types with compiler auto-vectorization hints
        for(int i = 0; i < output_size; ++i){
            Scalar sum = 0;
            
            // Unroll loop for better vectorization
            #pragma clang loop vectorize(enable) unroll(enable)
            for(int j = 0; j < input_size; ++j){
                sum += inputs[j] * weights[j * output_size + i];
            }
            sum += biases[i];
            activation_function(outputs[i], sum, alpha);
        }
    }
}

//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


template <typename Scalar = double>
inline auto MLP_LU(const std::array<Scalar, 85>& initial_input) {

    // Dense layer 1
    constexpr std::array<Scalar, 85> weights_1 = {2.314898942e-01, -2.360008097e-01, 2.132867318e-01, -6.455120126e-02, 1.586883213e-01, -1.360542058e-01, -1.732204110e-01, -1.759786338e-01, 2.231475962e-01, -2.129380108e-01, 1.301280961e-01, -1.580768727e-01, 1.328211401e-01, -1.754969708e-01, -1.711459106e-01, 9.202156739e-03, 5.552272755e-02, 1.863689419e-01, 1.077865183e-01, -4.995939357e-02, 4.722182535e-02, 1.587195453e-01, 1.118639994e-01, -5.845070075e-02, 4.530017281e-02, 7.340224637e-02, 3.429364944e-02, -2.029391651e-01, 2.395414525e-01, -3.660656256e-02, 1.011974795e-01, 1.238144227e-01, 1.849097818e-01, 1.059037864e-01, -1.871590837e-01, -4.491684968e-02, 1.715630475e-02, 3.172636893e-02, -9.167300981e-02, 8.671409470e-02, 2.589532899e-01, 2.550519900e-01, -2.273201791e-03, -3.697515224e-02, 1.165656594e-01, 2.549429517e-01, -1.679677033e-01, 2.526549053e-01, -9.570362503e-03, -2.286582277e-01, -2.132330966e-01, 2.606927488e-01, -2.209862662e-01, -1.184543654e-02, 1.166788954e-01, -2.485191479e-01, 1.701342637e-01, 1.070382915e-01, -2.465438932e-01, -3.504524127e-02, -9.671205561e-02, -2.597528005e-01, -7.058684764e-02, 1.481669608e-02, -1.838954712e-02, 1.004246063e-01, 2.060069023e-01, -7.475693104e-02, -9.233381554e-02, 2.366638035e-02, -1.591278912e-01, -7.698108727e-03, -2.065316863e-01, -2.324172769e-01, -8.046079263e-02, 9.713791325e-02, 1.203385602e-01, 1.028325071e-01, 1.934519967e-02, -2.098418534e-01, 1.104714443e-01, -1.626467454e-01, -1.637481850e-01, -2.430916480e-01, -1.193722821e-01};
    constexpr std::array<Scalar, 1> biases_1 = {1.060416778e-10};

    // Layer 2: Normalization
    constexpr Scalar epsilon_2 = 1.000000000e-05;

    // Dense layer 4
    constexpr std::array<Scalar, 85> weights_4 = {3.631022996e-01, 8.696174307e-05, -3.630182071e-02, -2.091381110e-02, -2.797661561e-02, 6.575272279e-01, -2.793872337e-02, -2.355301434e-02, -8.504942185e-01, 2.680931946e-02, -7.016853051e-01, 9.025993085e-01, 4.053229440e-01, 1.187000884e+00, 9.594328206e-01, -7.942843210e-01, -1.072281340e+00, 1.243301206e+00, -1.118940745e+00, -9.762162234e-01, -1.358493240e+00, -6.770197887e-01, -1.257935674e+00, -1.061007781e+00, 6.705536850e-01, 1.024484159e+00, 1.101222883e+00, 1.033919465e+00, 9.810726704e-01, 1.081941179e+00, 8.152692491e-01, -4.289987276e-02, 1.155492535e+00, -1.046623803e-01, -2.730516758e-01, -1.115444575e+00, -7.994265868e-01, 9.121985594e-01, 1.150044110e+00, 4.135038912e-01, -2.866357479e-01, 9.117669994e-01, -1.098062361e+00, -1.113841813e+00, 1.089138882e+00, -1.165388920e+00, 7.292374855e-02, -9.557927620e-01, 7.747234209e-01, 8.628410338e-01, -4.894282267e-01, 1.093913705e+00, 9.134582824e-01, -9.819466480e-01, 8.795045746e-01, -1.045818006e+00, -1.057007981e+00, 6.321521261e-01, -1.397918108e+00, -9.413889874e-01, 9.313437330e-01, 8.617879436e-01, -7.412077388e-01, 1.061172823e+00, 1.134024702e+00, -6.405270078e-01, -6.592480959e-01, 2.697742915e-01, 1.455900387e-01, 2.134715581e-01, 3.829560068e-02, -1.015404708e+00, -6.588732071e-01, -1.215095548e-01, 9.801791017e-01, -6.753806852e-01, -2.164453209e-01, 2.084377409e-01, -2.998163407e-01, -3.079006688e-01, -2.079924451e-01, 2.417955805e-01, 1.205455246e-01, -2.722151138e-01, 6.748557702e-01};
    constexpr std::array<Scalar, 85> biases_4 = {-2.402963693e-01, -1.258491762e-05, 9.071429495e-02, -1.123159873e-01, -1.413350004e-01, -4.079719855e-01, -1.418226768e-01, -1.437059997e-01, 4.596975259e-01, 1.418645873e-01, 4.276597859e-01, -4.823787875e-01, -3.189435044e-01, -5.213176973e-01, -4.599276610e-01, 4.277641588e-01, 5.233105064e-01, -6.293946503e-01, 4.600716745e-01, 4.736500604e-01, 6.322165755e-01, 4.089538327e-01, 5.797780154e-01, 4.953602486e-01, -3.821297931e-01, -4.767892673e-01, -4.795935889e-01, -4.064202943e-01, -4.374773187e-01, -4.420298543e-01, -4.529267967e-01, -1.449502838e-01, -5.513432815e-01, 4.696499982e-02, 1.438493454e-01, 4.815231282e-01, 4.077321847e-01, -4.665464583e-01, -5.564919939e-01, -3.138853363e-01, 2.834722465e-01, -3.244363617e-01, 4.447596679e-01, 4.564656240e-01, -4.711962707e-01, 4.946351700e-01, 9.205599346e-02, 4.170564597e-01, -4.147690596e-01, -3.586846012e-01, 2.542648606e-01, -4.428777290e-01, -4.873803942e-01, 4.242108851e-01, -4.366921985e-01, 4.962185360e-01, 4.286654782e-01, -3.899650192e-01, 6.556478547e-01, 4.194553740e-01, -4.669982409e-01, -4.466196885e-01, 3.094359207e-01, -4.203965425e-01, -5.254098436e-01, 1.494557326e-01, 3.908265234e-01, -1.214271466e-02, -2.612598721e-01, 6.690835314e-02, 1.254197174e-01, 4.308710354e-01, 3.282078698e-01, -1.052351770e-01, -3.882006580e-01, 4.115067135e-01, -7.991674846e-02, 8.768404418e-02, -3.569189095e-02, -1.551791215e-02, -9.013526137e-02, 7.296194440e-02, 1.190120193e-01, 2.638056924e-01, -3.856030815e-01};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 

    constexpr std::array<Scalar, 85> input_norm_std = {1.000000000e-16, 1.000000000e-16, 2.699196042e+09, 1.225573460e+09, 3.585557697e+05, 2.945345211e+11, 2.268399242e+11, 2.202964040e+11, 2.245129839e+11, 4.785749833e+11, 1.767507920e+11, 3.694415689e+13, 2.272483551e+14, 8.308113410e-01, 1.988106695e+06, 1.768274764e+06, 1.616055417e+06, 1.779017703e+02, 4.424210472e+06, 4.269940354e+05, 7.569700421e+06, 4.052924595e+06, 9.001187697e-01, 1.988108478e+06, 1.630712521e+06, 1.363674934e+06, 7.099660926e+05, 3.750293440e+06, 4.269886274e+05, 2.385414749e+07, 4.912584057e+06, 1.994398589e-01, 8.981420385e+05, 6.617463497e+05, 1.846904766e+06, 7.008693724e+05, 2.204679604e+06, 1.289861854e+05, 7.858080926e+05, 1.553042511e+06, 2.517413144e-01, 7.696207163e+01, 2.816978748e+05, 7.362852004e+05, 7.104179469e+05, 6.946001647e+05, 6.542859651e+01, 1.952537612e+07, 5.091584868e+03, 1.302251770e+00, 2.188139695e+05, 2.948529096e+05, 2.196400471e+06, 7.010350694e+05, 2.570917976e+06, 6.841591117e+05, 2.781134037e+07, 3.138461390e+08, 9.946365495e-01, 1.092438460e+06, 9.872170376e+05, 7.965732162e+05, 3.881127116e+02, 2.256592712e+06, 5.555667734e+05, 9.377280524e+06, 3.160177987e+08, 5.059059085e-03, 7.284679955e+01, 3.109994738e+03, 8.571206846e+02, 1.166401262e+04, 7.061261057e+03, 1.137319210e+02, 3.630818335e+07, 3.201108209e+08, 4.701282463e-03, 4.633427833e+00, 6.201303376e+01, 5.208729819e+01, 1.283541415e-01, 6.610950526e+03, 1.793710801e+02, 4.093875414e+06, 3.205435036e+08};

    constexpr std::array<Scalar, 85> input_min_mean = {1.000000000e+07, 1.000000000e+07, 5.193701600e+09, 2.320363673e+09, 1.079035772e+07, -5.642217886e+11, 6.602699036e+11, 1.195565605e+11, 5.172801725e+11, -3.703384732e+12, 3.451709759e+11, -9.045821782e+13, -3.194809897e+14, -1.026887167e+00, 1.377799693e+07, -2.895859241e+06, 5.292719737e+06, 2.032669319e+02, 2.959405476e+06, -8.156360346e+05, -1.828587162e+07, -9.273672690e+06, 1.651242366e+00, -3.777999288e+06, 1.376581281e+07, -6.981430077e+06, 1.616272685e+06, -4.529014644e+06, 8.156297585e+05, 5.723225405e+07, 1.155627815e+07, 5.843209498e-03, 1.721045000e+06, -2.151761742e+06, 1.824085827e+07, -1.590914028e+06, -4.233968955e+06, 2.209057257e+05, 1.586600862e+06, 3.222272971e+06, 6.285813221e-01, 9.437268313e+01, 8.447428861e+05, -1.692377185e+06, 1.161704042e+07, -1.600127699e+06, 5.806498943e+01, -4.402790013e+07, 7.263390460e+03, -2.943621933e+00, 3.357977044e+05, -1.090183711e+06, -6.113345739e+06, -1.591159043e+06, 2.623538313e+07, -1.257718969e+06, -6.221931519e+07, 4.259891371e+08, 1.671925162e+00, -2.056842762e+06, 1.559905030e+06, 1.253809284e+06, 3.229240400e+02, -8.815987349e+06, 1.103682059e+07, -1.910541493e+07, -4.337766141e+08, 6.387358409e-03, -8.274805722e+01, -3.724295439e+03, 1.716152376e+03, -2.616572849e+04, 8.032254794e+03, 1.561446281e+02, 9.802429617e+07, -4.440044982e+08, -2.041951975e-03, -1.162464279e+01, 1.604255538e+00, 8.817133054e-01, -1.824039574e-02, -6.183407352e+02, -2.176810422e+02, -4.127333019e+06, 4.562798139e+08};

    // Final output
    constexpr static std::array<Scalar, 85> output_norm_std = {2.698099241e+09, 1.000000000e-16, 1.360113272e-01, 7.290339559e-02, 2.067187864e+06, 8.415143423e+08, 1.058928678e+11, 8.043731986e+10, 6.048636129e+08, 9.618811632e+11, 5.021911140e+08, 8.817091838e+10, 6.473283789e+11, 4.983159985e-05, 1.938723542e+06, 1.708468495e+06, 1.628213640e+06, 4.166804763e+04, 4.650544217e+06, 3.972900846e+05, 1.311152054e+07, 3.940110523e+07, 4.749758191e-05, 1.205644031e-01, 9.008179545e+05, 1.239138828e+06, 6.724312708e+05, 4.015256189e+06, 2.400146283e+05, 2.744514515e+07, 3.654291018e+07, 1.267489856e-05, 5.725563465e-02, 2.075649334e-02, 1.372896094e+06, 6.081549860e+05, 1.783669229e+06, 2.349116128e+05, 7.186415731e+06, 2.108743302e+07, 2.801450101e-05, 1.017534390e-03, 2.708712737e-02, 5.817627235e-02, 5.952357928e+05, 1.436511234e+06, 8.927348095e+03, 1.792043965e+07, 1.153460472e+07, 6.679816340e-05, 7.983191316e-03, 2.615128917e-02, 1.142995679e-01, 7.075292680e-02, 1.375352150e+06, 4.854338910e+05, 3.793151204e+07, 2.468397973e+08, 4.754860844e-05, 6.345213516e-02, 3.917690573e-02, 7.482902599e-02, 1.883284556e-03, 1.353524450e-01, 2.237879655e+05, 1.486755607e+07, 2.099654335e+08, 2.941858653e-07, 2.337541398e-05, 2.677625029e-04, 1.620542530e-04, 1.007268423e-03, 4.499516690e-04, 3.705379793e-05, 3.629958072e+07, 3.199787182e+08, 1.770991055e-07, 2.344877583e-05, 1.667212969e-05, 1.824947319e-05, 2.414205974e-05, 3.345740600e-04, 3.533121121e-05, 3.286993798e-02, 2.909947461e+08};
    constexpr static std::array<Scalar, 85> output_min_mean = {5.194271710e+09, 1.000000000e+07, -7.060437424e-03, -4.345204434e-01, 3.313358789e+05, 1.196833998e+09, 1.669418380e+10, 1.341793233e+10, -1.047219969e+09, -1.566718909e+11, -7.292051093e+08, 1.813328742e+11, 7.203278109e+11, 4.369409514e-05, 1.370554761e+07, -2.823248859e+06, 5.303702824e+06, 6.323011614e+04, 2.626483312e+06, -7.714603998e+05, -2.891469499e+07, -5.282128498e+07, -7.951535687e-05, -2.509400389e-01, 1.275295588e+07, -5.736156457e+06, 1.541342586e+06, -3.844521726e+06, 5.092583172e+05, 6.493670374e+07, 5.626524899e+07, 3.144379795e-06, 1.178434785e-01, -1.345160738e-01, 1.686567884e+07, -1.389276309e+06, -4.791821253e+06, 4.034844034e+05, 1.447916948e+07, 2.313506187e+07, -3.603955638e-05, 2.408058109e-03, 6.503923786e-02, -7.761618482e-02, 1.135875945e+07, -1.488032668e+06, -4.453174487e+03, -4.012164927e+07, 1.832750808e+07, 1.466378762e-04, 1.034948481e-02, -6.804218321e-02, -3.812274185e-01, -1.585203713e-01, 2.269952403e+07, -9.490773225e+05, -8.891404780e+07, 3.361235565e+08, -7.715883217e-05, -1.330282380e-01, 7.311008808e-02, 1.333921450e-01, 1.506439074e-04, -3.247195464e-01, 1.046947332e+07, -3.873063203e+07, -2.960232188e+08, -3.512671432e-07, 2.026860565e-05, -3.239873404e-04, -2.348689318e-05, -2.257383512e-03, 1.447322702e-04, 2.469458129e-05, 9.801685728e+07, -4.437921553e+08, 7.911338681e-08, -1.069169766e-05, 8.832567378e-06, 9.113200621e-06, 1.038304821e-05, -4.785943738e-05, -2.234954521e-05, -3.585501508e-02, 4.301550828e+08};


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept 
    {
        output = input;
    };

    auto nonzero_diag = +[](Scalar& output, Scalar input, int index) noexcept
    {
        constexpr Scalar EPS = Scalar(1e-4);

        int r = get_lu_perm_row_index(index);
        int c = get_lu_perm_col_index(index);

        if (r == c) {
            Scalar abs_x  = std::abs(input);
            Scalar sign_x = (input >= Scalar(0) ? Scalar(1) : Scalar(-1));
            if (input == Scalar(0)) sign_x = Scalar(1);
            output = sign_x * std::max(abs_x, EPS);
        }
        else {
            output = input;
        }
    };

    static constexpr Scalar kC0 = 0.044715;
    static constexpr Scalar kSqrt2PiInv = Scalar(0.7978845608028654);
    auto gelu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept 
    {
        Scalar x3 = input * input * input;
        Scalar y  = kSqrt2PiInv * (input + kC0 * x3);
        output     = Scalar(0.5) * input * (Scalar(1) + std::tanh(y));
    };


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 

    
    // model input and flattened
    constexpr int flat_size = 85; 
    alignas(16) std::array<Scalar, flat_size> model_input;

    // Vectorized input normalization
    if constexpr (std::is_same_v<Scalar, float>) {
        // NEON optimized version for float32
        int i = 0;
        for(; i <= 85 - 4; i += 4) {
            float32x4_t input_vec = vld1q_f32(&initial_input[i]);
            float32x4_t mean_vec = vld1q_f32(&input_min_mean[i]);
            float32x4_t std_vec = vld1q_f32(&input_norm_std[i]);
            
            float32x4_t normalized = vdivq_f32(vsubq_f32(input_vec, mean_vec), std_vec);
            vst1q_f32(&model_input[i], normalized);
        }
        // Handle remaining elements
        for(; i < 85; i++) {
            model_input[i] = (initial_input[i] - input_min_mean[i]) / input_norm_std[i];
        }
    } else if constexpr (std::is_same_v<Scalar, double>) {
        // NEON optimized version for float64
        int i = 0;
        for(; i <= 85 - 2; i += 2) {
            float64x2_t input_vec = vld1q_f64(&initial_input[i]);
            float64x2_t mean_vec = vld1q_f64(&input_min_mean[i]);
            float64x2_t std_vec = vld1q_f64(&input_norm_std[i]);
            
            float64x2_t normalized = vdivq_f64(vsubq_f64(input_vec, mean_vec), std_vec);
            vst1q_f64(&model_input[i], normalized);
        }
        // Handle remaining elements
        for(; i < 85; i++) {
            model_input[i] = (initial_input[i] - input_min_mean[i]) / input_norm_std[i];
        }
    } else {
        // Fallback with compiler vectorization hints
        #pragma clang loop vectorize(enable) unroll(enable)
        for (int i = 0; i < 85; i++) { 
            model_input[i] = (initial_input[i] - input_min_mean[i]) / (input_norm_std[i]); 
        }
    } 

    if (model_input.size() != 85) { throw std::invalid_argument("Invalid input size. Expected size: 85"); }

    // Dense, layer 1
    alignas(16) static std::array<Scalar, 1> layer_1_output;
    Dense_MLP_LU<Scalar, 1>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        85, linear, 0.0);

    // UnitNormalization, layer 2
    alignas(16) static std::array<Scalar, 1> layer_2_output;
    UnitNormalization_MLP_LU<Scalar, 1>(
        layer_2_output.data(), layer_1_output.data(),
        epsilon_2);

    // Activation, layer 3 - Vectorized GELU application
    alignas(16) static std::array<Scalar, 1> layer_3_output;
    if constexpr (std::is_same_v<Scalar, float> && 1 >= 4) {
        // For larger arrays, use NEON vectorization
        int i = 0;
        for(; i <= 1 - 4; i += 4) {
            for(int j = 0; j < 4; ++j) {
                gelu(layer_3_output[i + j], layer_2_output[i + j], 0.0);
            }
        }
        for(; i < 1; ++i) {
            gelu(layer_3_output[i], layer_2_output[i], 0.0);
        }
    } else {
        // For small arrays or non-float types
        #pragma clang loop unroll(enable)
        for (int i = 0; i < 1; ++i) {
            gelu(layer_3_output[i], layer_2_output[i], 0.0);
        }
    }

    // Dense, layer 4
    alignas(16) static std::array<Scalar, 85> layer_4_output;
    Dense_MLP_LU<Scalar, 85>(
        layer_4_output.data(), layer_3_output.data(),
        weights_4.data(), biases_4.data(),
        1, linear, 0.0);

    // Activation, layer 5 - Vectorized nonzero_diag application
    alignas(16) static std::array<Scalar, 85> layer_5_output;
    if constexpr (std::is_same_v<Scalar, float> && 85 >= 4) {
        // For larger arrays, process with vectorization hints
        int i = 0;
        #pragma clang loop vectorize(enable) unroll(enable)
        for(; i < 85; ++i) {
            nonzero_diag(layer_5_output[i], layer_4_output[i], i);
        }
    } else {
        // For small arrays or non-float types
        #pragma clang loop unroll(enable)
        for (int i = 0; i < 85; ++i) {
            nonzero_diag(layer_5_output[i], layer_4_output[i], i);
        }
    }


//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\// 


    alignas(16) static std::array<Scalar, 85> model_output;

    // Vectorized output denormalization
    if constexpr (std::is_same_v<Scalar, float>) {
        // NEON optimized version for float32
        int i = 0;
        for(; i <= 85 - 4; i += 4) {
            float32x4_t layer_vec = vld1q_f32(&layer_5_output[i]);
            float32x4_t std_vec = vld1q_f32(&output_norm_std[i]);
            float32x4_t mean_vec = vld1q_f32(&output_min_mean[i]);
            
            float32x4_t denormalized = vfmaq_f32(mean_vec, layer_vec, std_vec);
            vst1q_f32(&model_output[i], denormalized);
        }
        // Handle remaining elements
        for(; i < 85; i++) {
            model_output[i] = (layer_5_output[i] * output_norm_std[i]) + output_min_mean[i];
        }
    } else if constexpr (std::is_same_v<Scalar, double>) {
        // NEON optimized version for float64
        int i = 0;
        for(; i <= 85 - 2; i += 2) {
            float64x2_t layer_vec = vld1q_f64(&layer_5_output[i]);
            float64x2_t std_vec = vld1q_f64(&output_norm_std[i]);
            float64x2_t mean_vec = vld1q_f64(&output_min_mean[i]);
            
            float64x2_t denormalized = vfmaq_f64(mean_vec, layer_vec, std_vec);
            vst1q_f64(&model_output[i], denormalized);
        }
        // Handle remaining elements
        for(; i < 85; i++) {
            model_output[i] = (layer_5_output[i] * output_norm_std[i]) + output_min_mean[i];
        }
    } else {
        // Fallback with compiler vectorization hints
        #pragma clang loop vectorize(enable) unroll(enable)
        for (int i = 0; i < 85; i++) { 
            model_output[i] = (layer_5_output[i] * output_norm_std[i]) + output_min_mean[i]; 
        }
    }

    return model_output;

}