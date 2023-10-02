### Adaptive quadrature in C#

The file `numerical_integration.cs` contains a standalone general-purpose adaptive integrator based on Gauss-Kronrod (G7K15) quadrature.

The accuracy of the integrator was validated on various functions against the integrator in `scipy`.

For the test functions and example use see `Program.cs`.

### Calling the integrator

```
using NumericalIntegration;

var integral = AdaptiveQuadrature.integrate(func, a, b, epsabs, epsrel, max_step);
Console.WriteLine("Integral: {0}, estimated error: {1}", integral.integral, integral.error);
```

The parameters are:
 - func: function double -> double that will be integrated
 - a: lower integration limit
 - b: upper integration limit
 - epsabs: tolerable absolute error
 - epsrel: tolerable relative error
 - max_step: size of maximum allowed step

The adaptive quadrature terminates when epsabs or epsrel are satisfied, or the required integration step is too small.

The parameter max_step can be useful when the function is very smooth except for a small region, which the integrator might miss when using large step size.