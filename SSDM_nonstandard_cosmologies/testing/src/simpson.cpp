#include <cmath>
#include <functional>
#include <array>
#include <iostream>
#include "simpson.h"

// Adaptación de r_simpson_arg para std::function
void r_simpson_arg(
    const std::function<double(double)>& func,
    double a, double b, double eps,
    std::array<double, 9>& f,
    double& ans, double& absAns, double& dErr,
    int depth, int depth1, int& nErr)
{
    std::array<double, 9> f1;
    double s1, s2, s3, e_err;

    s1 = (f[0] + 4 * f[4] + f[8]) / 6;
    s2 = (f[0] + 4 * f[2] + 2 * f[4] + 4 * f[6] + f[8]) / 12;
    s3 = (f[0] + 4 * f[1] + 2 * f[2] + 4 * f[3] + 2 * f[4] + 4 * f[5] + 2 * f[6] + 4 * f[7] + f[8]) / 24;

    e_err = eps * std::fabs(s3);

    int ok = 0;
    if (std::fabs(s3 - s2) <= e_err && std::fabs(s3 - s1) <= 16 * e_err) {
        absAns += std::fabs(f[0])+4*std::fabs(f[1])+2*std::fabs(f[2])+4*std::fabs(f[3])+2*std::fabs(f[4])+4*std::fabs(f[5])+2*std::fabs(f[6])+4*std::fabs(f[7])+std::fabs(f[8]);
        absAns *= std::fabs(b - a) / 24;
        ok = 1;
    } else if (std::fabs((s3 - s2) * (b - a)) <= 0.1 * eps * absAns &&
               std::fabs((s3 - s1) * (b - a)) <= 1.6 * eps * absAns) {
        ok = 1;
        absAns -= std::fabs((s3 - s2) * (b - a)) / eps;
    }

    if (ok) {
        ans += s3 * (b - a);
        return;
    }

    if (depth >= 10) {
        int c = 0, inc = f[0] < f[1];
        for (int i = 1; i < 7; ++i) {
            if (inc) { if (f[i] > f[i + 1]) { inc = 0; c++; } }
            else if (f[i] < f[i + 1]) { inc = 1; c++; }
        }
        if (c >= 3) {
            double ff = 0.5 * (f[0] + f[8]);
            for (int i = 1; i < 8; ++i) ff += f[i];
            ff /= 8;
            double ff2 = 0.5 * (f[0] * f[0] + f[8] * f[8]);
            for (int i = 1; i < 8; ++i) ff2 += f[i] * f[i];
            ff2 /= 8;
            double df = std::sqrt(ff2 - ff * ff) / std::pow(2, 0.5 * depth);
            if (df < 0.01 * std::fabs(ff)) {
                nErr |= 4;
                ans += ff * (b - a);
                return;
            }
        }
    }

    if (depth > depth1) {
        ans += s3 * (b - a);
        dErr += std::fabs((s3 - s2) * (b - a));
        return;
    }

    for (int i = 0; i < 5; ++i) f1[2 * i] = f[4 + i];
    for (int i = 8; i > 0; i -= 2) f[i] = f[i / 2];

    for (int i = 1; i < 8; i += 2) {
        f[i] = func(a + i * (b - a) / 16);
        if (!std::isfinite(f[i])) { f[i] = 0; nErr |= 1; }
        f1[i] = func((a + b) / 2 + i * (b - a) / 16);
        if (!std::isfinite(f1[i])) { f1[i] = 0; nErr |= 1; }
    }

    double as = 0, as1 = 0;
    for (int i = 0; i < 9; ++i) { as += std::fabs(f[i]); as1 += std::fabs(f1[i]); }
    if (std::fabs(as) > std::fabs(as1)) {
        r_simpson_arg(func, a, (a + b) / 2, eps, f, ans, absAns, dErr, depth + 1, depth1, nErr);
        r_simpson_arg(func, (a + b) / 2, b, eps, f1, ans, absAns, dErr, depth + 1, depth1, nErr);
    } else {
        r_simpson_arg(func, (a + b) / 2, b, eps, f1, ans, absAns, dErr, depth + 1, depth1, nErr);
        r_simpson_arg(func, a, (a + b) / 2, eps, f, ans, absAns, dErr, depth + 1, depth1, nErr);
    }
}

// Interfaz principal tipo simpson
double simpson(const std::function<double(double)>& func, double a, double b, double eps, int* err = nullptr)
{
    if (a == b) return 0;

    int nErr = 0;
    int depth1 = 50 - std::log(eps);

    double ans = 0, absAns = 0, dErr = 0;
    std::array<double, 9> f;
    for (int i = 0; i < 9; ++i) {
        f[i] = func(a + i * (b - a) / 8);
        if (!std::isfinite(f[i])) { f[i] = 0; nErr |= 1; }
    }

    r_simpson_arg(func, a, b, eps, f, ans, absAns, dErr, 0, depth1, nErr);

    if (dErr > 0.1 * eps * absAns) { nErr |= 2; }

    if (err) *err = nErr;
    if (nErr) {
        std::cout << "simpson warnings: ";
        if (nErr & 1) std::cout << "NaN in integrand; ";
        if (nErr & 2) std::cout << "Too deep recursion; ";
        if (nErr & 4) std::cout << "Lost of precision.";
        std::cout << std::endl;
    }
    return ans;
}