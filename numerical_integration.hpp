#pragma once

struct IntegrationResult{
    double integral;
    double error;
};

IntegrationResult integrate(double (*func)(double), const double a, const double b, const double epsabs, const double epsrel, const double max_step);
