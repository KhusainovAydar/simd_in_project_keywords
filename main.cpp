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

constexpr const size_t NMAX = 25000;
constexpr const size_t NMAXVEC = NMAX / 4;
constexpr const size_t LOOPS = 1000;


double a[NMAX + 10], b[NMAX + 10], c[NMAX + 10];
double ans[NMAX + 10];
Vec4d avec[NMAXVEC], bvec[NMAXVEC], cvec[NMAXVEC], ansvec[NMAXVEC];

void init() {
    double L = 50;
    double R = 100;
    for (size_t i = 0; i != NMAX; ++i) {
        a[i] = fRand(L, R);
        b[i] = fRand(L, R);
        c[i] = fRand(L, R);
    }
    size_t index = 0;
    for (size_t j = 0; j < NMAX; j += 4, ++index) {
        avec[index].load(&a[j]);
        bvec[index].load(&b[j]);
        cvec[index].load(&c[j]);
    }
}


void solve_opt_0() {
    for (size_t j = 0; j != NMAX; ++j)
        ans[j] = find_density_opt_0<double>(a[j], b[j], c[j]);
}


void solve_vec_opt_0() {
    for (size_t j = 0; j != NMAXVEC; ++j)
        ansvec[j] = find_density_opt_0<Vec4d>(avec[j], bvec[j], cvec[j]);
}


void solve_opt_1() {
    for (size_t j = 0; j != NMAX; ++j)
        ans[j] = find_density_opt_1<double>(a[j], b[j], c[j]);
}


void solve_direct_opt_1() {
    for (size_t j = 0; j != NMAX; ++j)
        ans[j] = find_direct_density_opt_1<double>(a[j], b[j], c[j]);
}


void solve_what_opt_1() {
    for (size_t j = 0; j != NMAX; ++j)
        ans[j] = find_minus_density_opt_1<double>(a[j], b[j], c[j]);
}


void solve_vec_opt_1() {
    for (size_t j = 0; j != NMAXVEC; ++j)
        ansvec[j] = find_density_opt_1<Vec4d>(avec[j], bvec[j], cvec[j]);
}


void solve_vec_opt_unroll_1() {
    constexpr const size_t block = 10;
    static_assert(NMAXVEC % block == 0, "block must divide NMAXVEC");
    for (size_t j = 0; j != NMAXVEC; j += 10) {
        ansvec[j] = find_density_opt_1<Vec4d>(avec[j], bvec[j], cvec[j]);
        ansvec[j + 1] = find_density_opt_1<Vec4d>(avec[j + 1], bvec[j + 1], cvec[j + 1]);
        ansvec[j + 2] = find_density_opt_1<Vec4d>(avec[j + 2], bvec[j + 2], cvec[j + 2]);
        ansvec[j + 3] = find_density_opt_1<Vec4d>(avec[j + 3], bvec[j + 3], cvec[j + 3]);
        ansvec[j + 4] = find_density_opt_1<Vec4d>(avec[j + 4], bvec[j + 4], cvec[j + 4]);
        ansvec[j + 5] = find_density_opt_1<Vec4d>(avec[j + 5], bvec[j + 5], cvec[j + 5]);
        ansvec[j + 6] = find_density_opt_1<Vec4d>(avec[j + 6], bvec[j + 6], cvec[j + 6]);
        ansvec[j + 7] = find_density_opt_1<Vec4d>(avec[j + 7], bvec[j + 7], cvec[j + 7]);
        ansvec[j + 8] = find_density_opt_1<Vec4d>(avec[j + 8], bvec[j + 8], cvec[j + 8]);
        ansvec[j + 9] = find_density_opt_1<Vec4d>(avec[j + 9], bvec[j + 9], cvec[j + 9]);
    }
}


void solve_opt_2() {
    for (size_t j = 0; j != NMAX; ++j)
        ans[j] = find_density_opt_2<double>(a[j], b[j], c[j]);
}


void solve_vec_opt_2() {
    for (size_t j = 0; j != NMAXVEC; ++j)
        ansvec[j] = find_density_opt_2<Vec4d>(avec[j], bvec[j], cvec[j]);
}


template<class TFunction>
double get_time_of(TFunction Function) {
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i != LOOPS; ++i)
        (Function)();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    return elapsed_seconds.count();
}


template<class TVector, class TFunction, class TContainer>
double get_time_of_2(TFunction method_of_compute, ContainerEnvironment<TContainer> &container_to_store, int loops) {
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
    cont.ptr_to_demension = c;
    cont.ptr_to_arg = a;
    cont.ptr_to_mean = b;
    cont.ptr_to_result = ans;
    cont.size_arg_result = 500;
    cont.size_mean_demen = 2500;


    std::cout << get_time_of_2<Vec4d>(find_density_opt_0<Vec4d>, cont, 1) << std::endl;
    std::cout << get_time_of_2<Vec4d>(find_density_opt_1<Vec4d>, cont, 1) << std::endl;
    std::cout << get_time_of_2<Vec4d>(find_density_opt_2<Vec4d>, cont, 1) << std::endl;

    std::cout << get_time_of_2<double>(find_density_opt_0<double>, cont, 1) << std::endl;
    std::cout << get_time_of_2<double>(find_density_opt_1<double>, cont, 1) << std::endl;
    std::cout << get_time_of_2<double>(find_density_opt_2<double>, cont, 1) << std::endl;

//    std::cout << "density_0 time = " << get_time_of(solve_0) << std::endl;
//    std::cout << "density_vec_0 time = " << get_time_of(solve_vec_0) << std::endl;
//    std::cout << "density_opt_0 time = " << get_time_of(solve_opt_0) << std::endl;
//    std::cout << "density_vec_opt_0 time = " << get_time_of(solve_vec_opt_0) << std::endl;
//    std::cout << "density_1 time = " << get_time_of(solve_1) << std::endl;
//    std::cout << "density_vec_1 time = " << get_time_of(solve_vec_1) << std::endl;
//    std::cout << "minus_density_opt_1 time = " << get_time_of(solve_what_opt_1) << std::endl;
//    std::cout << "density_direct_opt_1 time = " << get_time_of(solve_direct_opt_1) << std::endl;
//    std::cout << "density_opt_1 time = " << get_time_of(solve_opt_1) << std::endl;
//    std::cout << "density_vec_opt_1 time = " << get_time_of(solve_vec_opt_1) << std::endl;
//    std::cout << "density_vec_opt_unroll_1 time = " << get_time_of(solve_vec_opt_unroll_1) << std::endl;
//    std::cout << "density_2 time = " << get_time_of(solve_2) << std::endl;
//    std::cout << "density_vec_2 time = " << get_time_of(solve_vec_2) << std::endl;
//    std::cout << "density_opt_2 time = " << get_time_of(solve_opt_2) << std::endl;
//    std::cout << "density_vec_opt_2 time = " << get_time_of(solve_vec_opt_2) << std::endl;
    return 0;
}

