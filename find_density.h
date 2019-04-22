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


template<typename T>
constexpr bool is_vector_class =
    std::is_same<T, Vec4f>::value 
    || std::is_same<T, Vec2d>::value
    || std::is_same<T, Vec8f>::value
    || std::is_same<T, Vec4d>::value
//    || std::is_same<T, Vec16f>::value SPECIAL DEFINE
//    || std::is_same<T, Vec8d>::value  SPECIAL DEFINE
    ;


template <class T, class TFunction, class TContainer, class Enable = void>
struct global_density {
    global_density(TFunction method, ContainerEnvironment<TContainer> &container_to_store);
};


template <class T, class TFunction, class TContainer>
struct global_density<T, TFunction, TContainer, typename std::enable_if<std::is_floating_point<T>::value >::type> {
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
            container_to_store.ptr_to_result[i] /= static_cast<T>(container_to_store.size_arg_result);
        }

    }
};


template <class T, class TFunction, class TContainer>
struct global_density<T, TFunction, TContainer, typename std::enable_if_t<is_vector_class<T> > > {
    global_density(TFunction method, ContainerEnvironment<TContainer> &container_to_store) {

        size_t shift = T().size();
        for (size_t i = 0; i < container_to_store.size_arg_result; i += shift) {
            T ansVec(0);
            T argVec;
            argVec.load(&container_to_store.ptr_to_arg[i]);
            for (size_t j = 0; j < container_to_store.size_mean_dimension; j += shift) {
                T meanVec;
                meanVec.load(&container_to_store.ptr_to_mean[j]);
                T dimensionVec;
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

