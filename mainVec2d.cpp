#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#define VCL_FASTEXP
#include "vectorclass/vectormath_exp.h"
#include "vectorclass/vectorclass.h"
#include "find_density.h"



double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

constexpr const size_t shift_vec = 2;
constexpr const size_t NMAX = 25000;
constexpr const size_t NMAXVEC = NMAX / shift_vec;

double a[NMAX + 10], b[NMAX + 10], c[NMAX + 10];
double ans[NMAX + 10];
Vec2d avec[NMAXVEC], bvec[NMAXVEC], cvec[NMAXVEC], ansvec[NMAXVEC];

void init() {
    double L = 50;
    double R = 100;
    for (size_t i = 0; i != NMAX; ++i) {
        a[i] = fRand(L, R);
        b[i] = fRand(L, R);
        c[i] = fRand(L, R);
    }
    size_t index = 0;
    for (size_t j = 0; j < NMAX; j += shift_vec, ++index) {
        avec[index].load(&a[j]);
        bvec[index].load(&b[j]);
        cvec[index].load(&c[j]);
    }
}


template<class TVector, class TFunction, class TContainer>
double get_time_of_2(TFunction method_of_compute, ContainerEnvironment<TContainer> &container_to_store, size_t loops) {
    auto start = std::chrono::system_clock::now();

    for (size_t i = 0; i != loops; ++i) {
        global_density<TVector, TFunction, TContainer>(method_of_compute, container_to_store);
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    return elapsed_seconds.count() / static_cast<double>(loops);
}

int main() {
    init();

    // Uninitialized record type
    ContainerEnvironment<double> cont;
    cont.ptr_to_dimension = c;
    cont.ptr_to_arg = a;
    cont.ptr_to_mean = b;
    cont.ptr_to_result = ans;
    cont.size_arg_result = 500;
    cont.size_mean_dimension = 2500;

    size_t number_of_loops = 1;

    std::cout << "global_density_opt_0 time = "
              << get_time_of_2<double>(find_density_opt_0<double>, cont, number_of_loops) << std::endl;

    std::cout << "global_density_vec_opt_0 time = "
              << get_time_of_2<Vec2d>(find_density_opt_0<Vec2d>, cont, number_of_loops) << std::endl;

    std::cout << "global_density_opt_1 time = "
              << get_time_of_2<double>(find_density_opt_1<double>, cont, number_of_loops) << std::endl;

    std::cout << "global_density_vec_opt_1 time = "
              << get_time_of_2<Vec2d>(find_density_opt_1<Vec2d>, cont, number_of_loops) << std::endl;

    std::cout << "global_density_opt_2 time = "
              << get_time_of_2<double>(find_density_opt_2<double>, cont, number_of_loops) << std::endl;

    std::cout << "global_density_vec_opt_2 time = "
              << get_time_of_2<Vec2d>(find_density_opt_2<Vec2d>, cont, number_of_loops) << std::endl;

    return 0;
}

