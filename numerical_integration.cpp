#include "numerical_integration.hpp"

#include <queue>
#include <iostream>
#include <cmath>

struct IntegrationRegion
{
    double a;
    double b;
    double integral;
    double error;
};

constexpr double points[] = {
    0.991455371120813,
    0.949107912342759,
    0.864864423359769,
    0.741531185599394,
    0.586087235467691,
    0.405845151377397,
    0.207784955007898,
    0.000000000000000,
    -0.991455371120813,
    -0.949107912342759,
    -0.864864423359769,
    -0.741531185599394,
    -0.586087235467691,
    -0.405845151377397,
    -0.207784955007898
};

constexpr double gauss_weights[] = {
    0.0,
    0.129484966168870,
    0.0,
    0.279705391489277,
    0.0,
    0.381830050505119,
    0.0,
    0.417959183673469,
    0.0,
    0.129484966168870,
    0.0,
    0.279705391489277,
    0.0,
    0.381830050505119,
    0.0,
};

constexpr double kronrod_weights[] = {
    0.022935322010529,
    0.063092092629979,
    0.104790010322250,
    0.140653259715525,
    0.169004726639267,
    0.190350578064785,
    0.204432940075298,
    0.209482141084728,
    0.022935322010529,
    0.063092092629979,
    0.104790010322250,
    0.140653259715525,
    0.169004726639267,
    0.190350578064785,
    0.204432940075298
};

IntegrationRegion fixed_quadrature(double (*func)(double), const double a, const double b)
{
    double ba2 = (b - a) / 2;

    double gauss = 0;
    double kronrod = 0;
    for (int i = 0; i < 15; i++)
    {
        double x = (points[i] + 1) * ba2 + a;
        double val = func(x);
        gauss += val * gauss_weights[i];
        kronrod += val * kronrod_weights[i];
    }

    gauss *= ba2;
    kronrod *= ba2;

    double error = std::abs(kronrod - gauss);

    return IntegrationRegion{a, b, kronrod, error};
}

IntegrationResult integrate(
    double (*func)(double),
    const double a,
    const double b,
    const double epsabs,
    const double epsrel,
    const double max_step)
{
    // The quadrature is 15-point, thus adjust minimal and maximal step
    double max_region_size = max_step * 15;

    // Get priority queue for storing integration regions
    // Region with largest error pops first
    auto cmp = [](IntegrationRegion left, IntegrationRegion right) { return left.error < right.error; };
    std::priority_queue<IntegrationRegion, std::vector<IntegrationRegion>, decltype(cmp)> regions(cmp);

    double integral = 0;
    double error = 0;

    // Create initial integration regions according to maximum step size
    double ba = b - a;
    int init_regions = (int)std::ceil(ba / max_region_size);
    if (init_regions < 1)
        init_regions = 1;
    double init_region_size = ba / init_regions;
    double ai = a;
    double bi = ai + init_region_size;
    IntegrationRegion region;
    for (int i = 0; i < init_regions - 1; i++)
    {
        region = fixed_quadrature(func, ai, bi);
        regions.push(region);

        integral += region.integral;
        error += region.error;

        (ai, bi) = (bi, bi + init_region_size);
    }
    region = fixed_quadrature(func, ai, b);
    regions.push(region);
    integral += region.integral;
    error += region.error;

    // Check convergence and subdivide region with largest error
    while (error > epsabs && error / std::abs(integral) > epsrel)
    {
        const IntegrationRegion old_region = regions.top();
        regions.pop();
        

        // Check that the region size is large enough for numerical accuracy
        if (old_region.b - old_region.a < std::max(std::abs(old_region.a), std::abs(old_region.b)) * 1e-14)
        {
            std::cout << "WARNING: The integration step is smaller that spacing between numbers!" << std::endl;
            break;
        }
        double midpoint = (old_region.a + old_region.b) / 2;
        IntegrationRegion r1 = fixed_quadrature(func, old_region.a, midpoint);
        IntegrationRegion r2 = fixed_quadrature(func, midpoint, old_region.b);
        regions.push(r1);
        regions.push(r2);

        integral += r1.integral + r2.integral - old_region.integral;
        error += r1.error + r2.error - old_region.error;
    }

    std::cout << "Integration regions: " << regions.size() << std::endl;

    // Use the regions to recalculate integral and error to avoid rounding error
    integral = 0;
    error = 0;
    while (regions.size() > 0)
    {
        const IntegrationRegion final_region = regions.top();
        regions.pop();
        integral += final_region.integral;
        error += final_region.error;
    }

    return IntegrationResult{integral, error};
}