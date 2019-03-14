#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#define VCL_FASTEXP
#include "vectorclass/vectorclass.h"
#include "find_density.h"

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

const int NMAX = 25000;
const size_t LOOPS = 10000;


double a[NMAX + 10], b[NMAX + 10], c[NMAX + 10];
double ans[3][NMAX + 10];

void init() {
    double L = 50;
    double R = 100;
    for (size_t i = 0; i != NMAX; ++i) {
        a[i] = fRand(L, R);
        b[i] = fRand(L, R);
        c[i] = fRand(L, R);
    }
}

void print() {
    std::cout.precision(8);
    for (size_t i = 0; i != 3; ++i) {
        for (size_t j = 0; j != NMAX; ++j) {
            std::cout << std::fixed << ans[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

void solve() {
    for (size_t i = 0; i != 3; ++i) {
        for (size_t j = 0; j != NMAX; ++j) {
            switch (i) {
                case 0:
                    ans[i][j] = find_density_0<double>(a[j], b[j], c[j]);
                    break;
                case 1:
                    ans[i][j] = find_density_1<double>(a[j], b[j], c[j]);
                    break;
                case 2:
                    ans[i][j] = find_density_2<double>(a[j], b[j], c[j]);
                    break;
                default:
                    std::cerr << "ERROR\n";
                    break;
            }
        }
    }
}

void newsolve() {
    Vec4d avec, bvec, cvec, ansvec;
    for (size_t i = 0; i != 3; ++i) {
        for (size_t j = 0; j < NMAX; j += 4) {
            avec.load(&a[j]);
            bvec.load(&b[j]);
            cvec.load(&c[j]);
            switch (i) {
                case 0:
                    ansvec = find_density_0<Vec4d>(avec, bvec, cvec);
                    break;
                case 1:
                    ansvec = find_density_1<Vec4d>(avec, bvec, cvec);
                    break;
                case 2:
                    ansvec = find_density_2<Vec4d>(avec, bvec, cvec);
                    break;
                default:
                    std::cerr << "ERROR\n";
                    break;
            }
            ansvec.store(&ans[i][j]);
        }
    }
}

double get_time_without_simd() {
    auto start = std::chrono::system_clock::now();

    for (size_t i = 0; i != LOOPS; ++i) {
        solve();
        if (NMAX < 15 && LOOPS < 3) {
            print();
        }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    return elapsed_seconds.count();
}

double get_time_with_simd() {
    auto start = std::chrono::system_clock::now();

    for (size_t i = 0; i != LOOPS; ++i) {
        newsolve();
        if (NMAX < 15 && LOOPS < 3) {
            print();
        }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    return elapsed_seconds.count();
}


int main() {
    init();
    std::cout << std::fixed << get_time_without_simd() << '\n';
    std::cout << std::fixed << get_time_with_simd() << '\n';
    return 0;
}
