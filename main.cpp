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
double ans[3][NMAX + 10];
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


void solve_0() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_density_0<double>(a[j], b[j], c[j]);
}

void solve_opt_0() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_density_opt_0<double>(a[j], b[j], c[j]);
}


void solve_vec_0() {
  for (size_t j = 0; j != NMAXVEC; ++j)
    ansvec[j] = find_density_0<Vec4d>(avec[j], bvec[j], cvec[j]);
}

void solve_vec_opt_0() {
  for (size_t j = 0; j != NMAXVEC; ++j)
    ansvec[j] = find_density_opt_0<Vec4d>(avec[j], bvec[j], cvec[j]);
}


void solve_1() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_density_1<double>(a[j], b[j], c[j]);
}

void solve_opt_1() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_density_opt_1<double>(a[j], b[j], c[j]);
}

void solve_direct_opt_1() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_direct_density_opt_1<double>(a[j], b[j], c[j]);
}

void solve_what_opt_1() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_minus_density_opt_1<double>(a[j], b[j], c[j]);
}


void solve_vec_1() {
  for (size_t j = 0; j != NMAXVEC; ++j)
    ansvec[j] = find_density_1<Vec4d>(avec[j], bvec[j], cvec[j]);
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


void solve_2() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_density_2<double>(a[j], b[j], c[j]);
}

void solve_opt_2() {
  for (size_t j = 0; j != NMAX; ++j)
    ans[0][j] = find_density_opt_2<double>(a[j], b[j], c[j]);
}

void solve_vec_2() {
  for (size_t j = 0; j != NMAXVEC; ++j)
    ansvec[j] = find_density_2<Vec4d>(avec[j], bvec[j], cvec[j]);
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

void check() {
  double res, resopt;
  for (size_t j = 0; j != NMAX; ++j) {
    res = find_density_0<double>(a[j], b[j], c[j]);
    resopt = find_density_opt_0<double>(a[j], b[j], c[j]);
    if (abs(res - resopt) > 1e-17)
      std::cout << "density 0\n"
                << "j = " << j
                << " res = " << res
                << " resopt = " << resopt << std::endl;
  }
  for (size_t j = 0; j != NMAX; ++j) {
    res = find_density_1<double>(a[j], b[j], c[j]);
    resopt = find_density_opt_1<double>(a[j], b[j], c[j]);
    if (abs(res - resopt) > 1e-19)
      std::cout << "density 1\n"
                << "j = " << j
                << " res = " << res
                << " resopt = " << resopt << std::endl;
  }
  for (size_t j = 0; j != NMAX; ++j) {
    res = find_density_1<double>(a[j], b[j], c[j]);
    resopt = find_minus_density_opt_1<double>(a[j], b[j], c[j]);
    if (abs(res - resopt) > 1e-19)
      std::cout << "density what 1\n"
                << "j = " << j
                << " res = " << res
                << " resopt = " << resopt << std::endl;
  }
  for (size_t j = 0; j != NMAX; ++j) {
    res = find_density_2<double>(a[j], b[j], c[j]);
    resopt = find_density_opt_2<double>(a[j], b[j], c[j]);
    if (abs(res - resopt) > 1e-20)
      std::cout << "density 2\n"
                << "j = " << j
                << " res = " << res
                << " resopt = " << resopt << std::endl;
  }


}

int main() {
  init();
  check();
  //std::cout << std::fixed << get_time_without_simd() << '\n';
  //std::cout << std::fixed << get_time_with_simd() << '\n';
  std::cout << "density_0 time = " << get_time_of(solve_0) << std::endl;
  std::cout << "density_vec_0 time = " << get_time_of(solve_vec_0) << std::endl;
  std::cout << "density_opt_0 time = " << get_time_of(solve_opt_0) << std::endl;
  std::cout << "density_vec_opt_0 time = " << get_time_of(solve_vec_opt_0) << std::endl;
  std::cout << "density_1 time = " << get_time_of(solve_1) << std::endl;
  std::cout << "density_vec_1 time = " << get_time_of(solve_vec_1) << std::endl;
  std::cout << "minus_density_opt_1 time = " << get_time_of(solve_what_opt_1) << std::endl;
  std::cout << "density_direct_opt_1 time = " << get_time_of(solve_direct_opt_1) << std::endl;
  std::cout << "density_opt_1 time = " << get_time_of(solve_opt_1) << std::endl;
  std::cout << "density_vec_opt_1 time = " << get_time_of(solve_vec_opt_1) << std::endl;
  std::cout << "density_vec_opt_unroll_1 time = " << get_time_of(solve_vec_opt_unroll_1) << std::endl;
  std::cout << "density_2 time = " << get_time_of(solve_2) << std::endl;
  std::cout << "density_vec_2 time = " << get_time_of(solve_vec_2) << std::endl;
  std::cout << "density_opt_2 time = " << get_time_of(solve_opt_2) << std::endl;
  std::cout << "density_vec_opt_2 time = " << get_time_of(solve_vec_opt_2) << std::endl;
  return 0;
}
