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
#define one_over_sqrt_two_pi 0.398942280401432702863218082712

template<typename T>
static inline T find_density_opt_0(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
    T arg_minus_mean = x - mean;

    return exp(-0.5 * arg_minus_mean * arg_minus_mean
               * one_over_div * one_over_div)
           * one_over_div * one_over_sqrt_two_pi;
}


template<typename T>
static inline T find_density_opt_1(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
    T arg_minus_mean = x - mean;

    return - exp(-0.5 * arg_minus_mean * arg_minus_mean
                 * one_over_div * one_over_div)
           * one_over_sqrt_two_pi
           * one_over_div * one_over_div * one_over_div
           * arg_minus_mean;
}

template<typename T>
static inline T find_density_opt_2(T &x, T &mean, T &deviation) {
    T one_over_div = 1. / deviation;
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
    T* ptr_to_dimension;
    T* ptr_to_result;
    size_t size_arg_result;
    size_t size_mean_dimension;
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
            for (size_t j = 0; j != container_to_store.size_mean_dimension; ++j) {
                container_to_store.ptr_to_result[i] +=
                        method(
                                container_to_store.ptr_to_arg[i],
                                container_to_store.ptr_to_mean[j],
                                container_to_store.ptr_to_dimension[j]
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

        constexpr size_t shift = 4;
        for (size_t i = 0; i < container_to_store.size_arg_result; i += shift) {
            Vec4d ansVec(0);
            Vec4d argVec;
            argVec.load(&container_to_store.ptr_to_arg[i]);
            for (size_t j = 0; j < container_to_store.size_mean_dimension; j += shift) {
                Vec4d meanVec;
                meanVec.load(&container_to_store.ptr_to_mean[j]);
                Vec4d dimensionVec;
                dimensionVec.load(&container_to_store.ptr_to_dimension[j]);
                ansVec += method(
                        argVec,
                        meanVec,
                        dimensionVec
                );
            }
            ansVec /= static_cast<double>(container_to_store.size_arg_result);
            ansVec.store(&container_to_store.ptr_to_result[i]);
        }
    }
};

#endif //PROJECT_FIND_DENSITY_H

