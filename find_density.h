//
// Created by MacBook on 13.03.2019.
//

#ifndef PROJECT_FIND_DENSITY_H
#define PROJECT_FIND_DENSITY_H

#include <math.h>
#define VCL_FASTEXP
#include "vectorclass/vectormath_exp.h"
#include "vectorclass/vectorclass.h"

#define M_PI 3.14159265358979323846264338327950288


template<typename T>
static inline T find_density_opt_0(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
    T one_over_sqrt_two_pi = 1. / sqrt(2. * M_PI);
    T arg_minus_mean = x - mean;

    return exp(-0.5 * arg_minus_mean * arg_minus_mean
               * one_over_div * one_over_div)
           * one_over_div * one_over_sqrt_two_pi;
}


template<typename T>
static inline T find_density_opt_1(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
    T one_over_sqrt_two_pi = 1. / sqrt(2. * M_PI);
    T arg_minus_mean = x - mean;

    return - exp(-0.5 * arg_minus_mean * arg_minus_mean
                 * one_over_div * one_over_div)
           * one_over_sqrt_two_pi
           * one_over_div * one_over_div * one_over_div
           * arg_minus_mean;
}

template<typename T>
static inline T find_minus_density_opt_1(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
    T one_over_sqrt_two_pi = 1. / sqrt(2. * M_PI);
    //T arg_minus_mean = x - mean;
    T mean_minus_arg = mean - x;

    return exp(-0.5 * mean_minus_arg * mean_minus_arg
               * one_over_div * one_over_div)
           * one_over_sqrt_two_pi
           * one_over_div * one_over_div * one_over_div
           * mean_minus_arg;
}

template<typename T>
static inline T find_direct_density_opt_1(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
    return exp(-0.5 * (x - mean) * (x - mean)
               * one_over_div * one_over_div)
           * one_over_div * one_over_div * one_over_div
           * (mean - x) / sqrt(2. * M_PI);
}

template<typename T>
static inline T find_density_opt_2(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
    T one_over_sqrt_two_pi = 1. / sqrt(2. * M_PI);
    T arg_minus_mean = x - mean;


    return exp(-0.5 * arg_minus_mean * arg_minus_mean
               * one_over_div * one_over_div)
           * one_over_sqrt_two_pi
           * one_over_div * one_over_div * one_over_div
           * (1. - arg_minus_mean * arg_minus_mean
                   * one_over_div * one_over_div);

}

template<class T>
struct ContainerEnvironment {
    T* ptr_to_arg;
    T* ptr_to_mean;
    T* ptr_to_demension;
    T* ptr_to_result;
    size_t size_arg_result;
    size_t size_mean_demen;
};


template <typename T, class TFunction, class TContainer>
struct global_density {
    global_density(TFunction method, ContainerEnvironment<TContainer> &container_to_store);
};

//@todo как сделать чтобы можно было выбрать один из типов double, float...
template <class TFunction, class TContainer>
struct global_density<double, TFunction, TContainer> {
    global_density(TFunction method, ContainerEnvironment<TContainer> &container_to_store) {

        for (size_t i = 0; i != container_to_store.size_arg_result; ++i) {
            container_to_store.ptr_to_result[i] = 0;
            for (size_t j = 0; j != container_to_store.size_mean_demen; ++j) {
                container_to_store.ptr_to_result[i] +=
                        method(
                                container_to_store.ptr_to_arg[i],
                                container_to_store.ptr_to_mean[j],
                                container_to_store.ptr_to_demension[j]
                        );
            }
            container_to_store.ptr_to_result[i] /= static_cast<double>(container_to_store.size_arg_result);
        }

    }
};


//@todo как сделать чтобы можно было выбрать один из типов Vec4d, Vec8d, Vec8f...
template <class TFunction, class TContainer>
struct global_density<Vec4d, TFunction, TContainer> {
    global_density(TFunction method, ContainerEnvironment<TContainer> &container_to_store) {

        int shift = 4;
        Vec4d argVec(0), meanVec(0), demensionVec(0);
        for (size_t i = 0; i < container_to_store.size_arg_result; i += shift) {
            Vec4d ansVec(0);
            argVec.load(&container_to_store.ptr_to_arg[i]);
            for (size_t j = 0; j < container_to_store.size_mean_demen; j += shift) {
                meanVec.load(&container_to_store.ptr_to_mean[j]);
                demensionVec.load(&container_to_store.ptr_to_mean[j]);
                ansVec += method(
                        argVec,
                        meanVec,
                        demensionVec
                );
            }
            ansVec /= static_cast<double>(container_to_store.size_arg_result);
            ansVec.store(&container_to_store.ptr_to_result[i]);
        }
    }
};

#endif //PROJECT_FIND_DENSITY_H

