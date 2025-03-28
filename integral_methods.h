#pragma once
#include <complex>
#include <functional>

namespace integral_methods {
constexpr static double GAUSSIAN_FACTOR_L = 1 - 1/std::sqrt(3);
constexpr static double GAUSSIAN_FACTOR_G = 1 + 1/std::sqrt(3);

inline double calculate_left_integral(
    std::vector<double> dots, const std::function<double(double)> &f)
{
    double result = 0;
    for (int i = 0; i < dots.size() - 1; i++) {
        result += (dots[i + 1] - dots[i]) * f(dots[i]);
    }

    return result;
}

inline double calculate_right_integral(
    std::vector<double> dots, const std::function<double(double)> &f)
{
    double result = 0;
    for (int i = 1; i < dots.size(); i++) {
        result += (dots[i + 1] - dots[i]) * f(dots[i + 1]);
    }

    return result;
}

inline double calculate_central_integral(
    std::vector<double> dots, const std::function<double(double)> &f)
{
    double result = 0;
    for (int i = 0; i < dots.size() - 1; i++) {
        auto left = dots[i];
        auto right = dots[i + 1];
        auto central = (left + right) / 2;

        result += (right - left) * f(central);
    }

    return result;
}

inline double calculate_trapezoid_integral(
    std::vector<double> dots, const std::function<double(double)> &f)
{
    double result = 0;
    for (int i = 0; i < dots.size() - 1; i++) {
        auto left = dots[i];
        auto right = dots[i + 1];

        result += (right - left) * (f(left) + f(right)) / 2;
    }

    return result;
}

inline double calculate_simpson_integral(
    std::vector<double> dots, const std::function<double(double)> &f)
{
    double result = 0;
    for (int i = 0; i < dots.size() - 1; i++) {
        auto left = dots[i];
        auto right = dots[i + 1];
        auto f_factor = f(left) + 4 * f((left + right) / 2) + f(right);

        result += (right - left) / 6 * f_factor;
    }

    return result;
}

inline double calculate_gaussian_integral(
    std::vector<double> dots, const std::function<double(double)> &f)
{
    double result = 0;

    for (int i = 0; i < dots.size() - 1; i++) {
        auto left = dots[i];
        auto right = dots[i + 1];
        auto h_half = (right - left) / 2;
        auto f_factor = f(left + h_half * GAUSSIAN_FACTOR_L)
                        + f(left + h_half * GAUSSIAN_FACTOR_G);

        result += h_half * f_factor;
    }

    return result;
}
} // namespace integral_methods
