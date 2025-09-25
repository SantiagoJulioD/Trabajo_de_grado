#include <cmath>
#include <functional>
#include <array>
#include <iostream>

void r_simpson_arg(
    const std::function<double(double)>& func,
    double a, double b, double eps,
    std::array<double, 9>& f,
    double& ans, double& absAns, double& dErr,
    int depth, int depth1, int& nErr);

double simpson(const std::function<double(double)>& func, double a, double b, double eps, int* err);

