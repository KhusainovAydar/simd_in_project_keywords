//
// Created by MacBook on 13.03.2019.
//

#ifndef PROJECT_FIND_DENSITY_H
#define PROJECT_FIND_DENSITY_H

#include <iostream>
#include <cmath>
#include <chrono>
#include "vectorclass/vectorclass.h"
#include "vectorclass/vectormath_exp.h"

template<typename T>
T find_density_0(T &x, T &mean, T &deviation) {
    T sigma = deviation;
    T ret = 0;

    T my_exp = (x - mean);
    my_exp *= -my_exp;
    my_exp /= 2 * sigma * sigma;
    my_exp = exp(my_exp);

    ret = my_exp;
    ret /= sqrt(2 * M_PI) * sigma;

    return ret;
}

template<typename T>
T find_density_1(T &x, T &mean, T &deviation) {
    T sigma = deviation;
    T ret = 0;

    T my_exp = (x - mean);
    my_exp *= -my_exp;
    my_exp /= 2 * sigma * sigma;

    my_exp = exp(my_exp);


    ret = my_exp * (x - mean);
    ret /= sqrt(2 * M_PI);

    ret /= sigma * sigma;
    ret /= sigma;

    return ret;
}

template<typename T>
T find_density_2(T &x, T &mean, T &deviation) {
    auto start = std::chrono::system_clock::now();

    T sigma = deviation;

    T my_exp = -(x - mean) * (x - mean) / 2;
    my_exp /= sigma * sigma;
    my_exp = exp(my_exp);

    T ret1 = (x - mean) * (x - mean);
    ret1 *= my_exp;

    ret1 /= sqrt(2 * M_PI);

    ret1 /= sigma * sigma;
    ret1 /= sigma * sigma;
    ret1 /= sigma;

    T ret2 = -my_exp;

    ret2 /= sqrt(2 * M_PI);


    ret2 /= sigma * sigma;
    ret2 /= sigma;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    if (std::is_same<T, double>::value) {
//        std::cerr << "DOUBLE ";
    } else {
//        std::cerr << "VECTOR " << std::fixed << elapsed_seconds.count() << '\n';
    }
//    std::cerr << std::fixed << elapsed_seconds.count() << '\n';
    return ret1 - ret2;
}


#endif //PROJECT_FIND_DENSITY_H
