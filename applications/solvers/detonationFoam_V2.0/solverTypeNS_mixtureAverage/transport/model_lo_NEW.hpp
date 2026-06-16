#pragma once
#include <iostream>
#include <array>
#include <random>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <limits>

/*\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//*/ 


template<typename Scalar, int output_size, typename ActFun>
inline void Dense_model_lo_NEW(Scalar* __restrict outputs, const Scalar* __restrict inputs, const Scalar * __restrict weights, const Scalar * __restrict biases, int input_size, ActFun activation_function, Scalar alpha) noexcept 
{
    for(int i = 0; i < output_size; ++i){
        Scalar sum = 0;
        
        for(int j = 0; j < input_size; ++j){
            sum += inputs[j] * weights[j * output_size + i];
        }
        sum += biases[i];
        activation_function(outputs[i], sum, alpha);
    }
}


/*\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//*/ 


template <typename Scalar = double>
inline auto model_lo_NEW(const std::array<Scalar, 4>& initial_input) {

    // Dense layer 1
    constexpr std::array<Scalar, 4> weights_1 = {7.00837731361389160e-01, 7.71281242370605469e-01, 2.80559480190277100e-01, 9.29196238517761230e-01};
    constexpr std::array<Scalar, 1> biases_1 = {1.31508278846740723e+00};

    // Dense layer 2
    constexpr std::array<Scalar, 1> weights_2 = {8.66168081760406494e-01};
    constexpr std::array<Scalar, 1> biases_2 = {-1.24522316455841064e+00};


/*\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//*/ 


    constexpr std::array<Scalar, 4> input_scale = {1.77591421598875004e-01, 1.98512375693752019e-01, 2.42619805685777096e-01, 7.55877352536265448e+02};

    constexpr std::array<Scalar, 4> input_shift = {2.46778118661087148e-01, 3.51785758261501325e-01, 4.01436120647625239e-01, 1.33417697899999303e+03};

    constexpr std::array<Scalar, 1> output_scale = {1.69752328136196334e-06};

    constexpr std::array<Scalar, 1> output_shift = {3.56793146875000022e-06};


/*\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//*/ 


    auto relu = +[](Scalar& output, Scalar input, Scalar alpha) noexcept 
    {
        output = input > 0 ? input : 0;
    };

    auto linear = +[](Scalar& output, Scalar input, Scalar alpha) noexcept 
    {
        output = input;
    };


/*\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//*/ 

    
    // model input and flattened
    constexpr int flat_size = 4; 
    std::array<Scalar, flat_size> model_input;

    // normalize input
    for (int i = 0; i < 4; i++) { model_input[i] = (initial_input[i] - input_shift[i]) / (input_scale[i]); } 

    if (model_input.size() != 4) { throw std::invalid_argument("Invalid input size. Expected size: 4"); }

    // Dense, layer 1
    static std::array<Scalar, 1> layer_1_output;
    Dense_model_lo_NEW<Scalar, 1>(
        layer_1_output.data(), model_input.data(),
        weights_1.data(), biases_1.data(),
        4, relu, 0.0);

    // Dense, layer 2
    static std::array<Scalar, 1> layer_2_output;
    Dense_model_lo_NEW<Scalar, 1>(
        layer_2_output.data(), layer_1_output.data(),
        weights_2.data(), biases_2.data(),
        1, linear, 0.0);


/*\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//\//*/ 


    static std::array<Scalar, 1> model_output;

    for (int i = 0; i < 1; i++) { model_output[i] = (layer_2_output[i] * output_scale[i]) + output_shift[i]; }

    return model_output;

}