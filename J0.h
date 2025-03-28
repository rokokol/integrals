#pragma once
#include "integral_methods.h"
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

namespace J0 {
constexpr static double FACTOR = 2 / std::sqrt(M_PI);

inline double q(double x, int N)
{
    double numerator = x * x * (2 * N + 1);
    double denominator = (N + 1) * (2 * N + 3);

    return -numerator / denominator;
}

inline double erf(double x, double e = 1E-6)
{
    double a = x;
    double result = a;
    int n = 0;

    while (std::abs(a) > e) {
        a *= q(x, n);
        result += a;
        n += 1;
    }

    return FACTOR * result;
}

inline double f(double t)
{
    return std::exp(-t * t);
}

inline std::vector<double> linspace(double a, double b, int N)
{
    std::vector<double> result;
    result.resize(N + 1);

    double step = (b - a) / static_cast<double>(N);

    result[0] = a;
    for (int i = 1; i <= N; i++) {
        result[i] = a + step * i;
    }
    return result;
}

inline std::pair<ulong, double> calculate_integral(
    const std::function<double(std::vector<double>, std::function<double(double)>)> &method,
    const std::function<double(double, double)> &j0,
    const std::function<double(double)> &function,

    double a = 0,
    double x = 2,
    double e = 1E-6,
    double factor = 1)
{
    int n = 1;
    ulong i = 0;
    std::vector<double> dots = linspace(a, x, n);
    double approx = j0(x, e);
    double result = factor * method(dots, function);
    while (std::abs(approx - result) > e) {
        i++;
        n *= 2;
        dots = linspace(a, x, n);
        result = factor * method(dots, function);
    }

    return {n, result};
}

inline void test_integral_method_with_erf(
    const std::function<double(std::vector<double>, std::function<double(double)>)> &integral)
{
    printf("| $x$    | $J_0(x)$    | $J_n(x)$    | $n$    |\n");
    printf("|:------:|:----------:|:----------:|:------|\n");
    for (auto &&t : linspace(0, 2, 10)) {
        auto value = erf(t);
        auto res = calculate_integral(integral, erf, f, 0, t, 1E-6, FACTOR);
        printf("|$%-2.2f$ | $%-6.6f$ | $%-6.6f$ | $%-8lu$|\n", t, value, res.second, res.first);
    }
}
} // namespace J0
