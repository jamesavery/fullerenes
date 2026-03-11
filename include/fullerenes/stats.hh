// Two-pass variance/standard deviation for numerical stability.
// The one-pass formula Var = E[X^2] - E[X]^2 suffers from catastrophic
// cancellation when the variance is small relative to the mean squared.
#pragma once

#include <cmath>
#include <numeric>
#include <vector>
#include <functional>

// Standard deviation of values produced by a generator function.
// Calls gen(i) for i in [0, n) twice: once for the mean, once for deviations.
template<typename F>
double stddev_twopass(int n, F gen) {
    if (n <= 0) return 0;
    double sum = 0;
    for (int i = 0; i < n; i++) sum += gen(i);
    double mean = sum / n;
    double var = 0;
    for (int i = 0; i < n; i++) { double d = gen(i) - mean; var += d*d; }
    return std::sqrt(var / n);
}

// CV = stddev / mean.
template<typename F>
double cv_twopass(int n, F gen) {
    if (n <= 0) return 0;
    double sum = 0;
    for (int i = 0; i < n; i++) sum += gen(i);
    double mean = sum / n;
    if (mean == 0) return 0;
    double var = 0;
    for (int i = 0; i < n; i++) { double d = gen(i) - mean; var += d*d; }
    return std::sqrt(var / n) / mean;
}

// Standard deviation of a vector.
inline double stddev_twopass(const std::vector<double>& v) {
    return stddev_twopass((int)v.size(), [&](int i){ return v[i]; });
}

// CV of a vector.
inline double cv_twopass(const std::vector<double>& v) {
    return cv_twopass((int)v.size(), [&](int i){ return v[i]; });
}
