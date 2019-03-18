//
// Created by MacBook on 13.03.2019.
//

#ifndef PROJECT_FIND_DENSITY_H
#define PROJECT_FIND_DENSITY_H


//#define M_PI 3.14159265358979323846264338327950288

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
static inline T find_density_opt_0(T &x, T &mean, T &deviation) {
  T one_over_div = 1. / deviation;
  T one_over_sqrt_two_pi = 1. / sqrt(2. * M_PI);
  T arg_minus_mean = x - mean;

  return exp(-0.5 * arg_minus_mean * arg_minus_mean
             * one_over_div * one_over_div)
         * one_over_div * one_over_sqrt_two_pi;
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
  // добавил знак, который был потерян
  // Теперь функция правильная
  return -ret;
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
T find_density_2(T &x, T &mean, T &deviation) {

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

  // убрал -1
  T ret2 = my_exp;

  ret2 /= sqrt(2 * M_PI);


  ret2 /= sigma * sigma;
  ret2 /= sigma;
  // поменял знак
  // Теперь функция правильная
  return ret2 - ret1;
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
#endif //PROJECT_FIND_DENSITY_H
