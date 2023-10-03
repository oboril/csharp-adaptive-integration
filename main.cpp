#include <iostream>
#include <chrono>
#include <cmath>

#include "numerical_integration.hpp"

using namespace std;

double func1(const double x)
{
    double f = std::pow(std::fmod(x, 1234) / 100, 2);
    f += std::sqrt(x) * 0.2;
    f += std::sin(x / 58 * std::exp(-x / 2000)) * 30;
    return f;
}


double func2(double x)
{
    double f = 1e-5 + std::pow(x, 2);
    return f;
}

double func3(double x)
{
    double f = 1 / (1 + std::fmod(x, 12));
    f += std::floor(std::fmod(x, 23) / 5) * 0.3;
    return f;
}

double func4(double x)
{
    double f = std::sin(x) * x;
    return f;
}

int main()
{
    constexpr double lims1[] = {14.9, 4534.453};
    constexpr double lims2[] = {-1e-6, 3e-6};
    constexpr double lims3[] = {474564, 474599};
    constexpr double lims4[] = {-50, 50};

    cout << "Hello world!" << endl;

    auto start = chrono::high_resolution_clock().now();

    const IntegrationResult integral1 = integrate(func1, lims1[0], lims1[1], 0, 1e-13, 99e99);
    const IntegrationResult integral2 = integrate(func2, lims2[0], lims2[1], 0, 1e-13, 99e99);
    const IntegrationResult integral3 = integrate(func3, lims3[0], lims3[1], 0, 1e-13, 99e99);
    const IntegrationResult integral4 = integrate(func4, lims4[0], lims4[1], 0, 1e-13, 99e99);

    auto end = chrono::high_resolution_clock().now();

    cout.precision(15);
    cout << integral1.integral << ", error " << integral1.error << endl;
    cout << integral2.integral << ", error " << integral2.error << endl;
    cout << integral3.integral << ", error " << integral3.error << endl;
    cout << integral4.integral << ", error " << integral4.error << endl;

    cout.precision(5);
    cout << "Elapsed: " << (end - start).count() / 1000000.0 << "ms" << endl;
}