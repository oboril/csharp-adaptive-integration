using NumericalIntegration;

var lims1 = (14.9, 4534.453);
double func1(double x)
{
    double f = Math.Pow((x % 1234) / 100, 2);
    f += Math.Sqrt(x) * 0.2;
    f += Math.Sin(x / 58 * Math.Exp(-x / 2000)) * 30;
    return f;
}


var lims2 = (-1e-6, 3e-6);
double func2(double x)
{
    double f = 1e-5 + Math.Pow(x, 2);
    return f;
}

var lims3 = (474564, 474599);
double func3(double x)
{
    double f = 1 / (1 + (x % 12));
    f += Math.Floor((x % 23) / 5) * 0.3;
    return f;
}

var integral1 = AdaptiveQuadrature.integrate(func1, lims1.Item1, lims1.Item2, 0, 1e-13, 99e99);
Console.WriteLine("{0} {1}", integral1.integral, integral1.error);

var integral2 = AdaptiveQuadrature.integrate(func2, lims2.Item1, lims2.Item2, 0, 1e-13, 99e99);
Console.WriteLine("{0} {1}", integral2.integral, integral2.error);

var integral3 = AdaptiveQuadrature.integrate(func3, lims3.Item1, lims3.Item2, 0, 1e-13, 99e99);
Console.WriteLine("{0} {1}", integral3.integral, integral3.error);


// OUTPUT FROM C#:
// Integration regions: 112
// 219799.25080829102 1.611878881396469E-08
// Integration regions: 1
// 4.000000933333322E-11 1.357093192469811E-25
// WARNING: The integration step is smaller that spacing between numbers!
// Integration regions: 266
// 26.514805362870966 1.3626203079323724E-09

// (OUTPUT FROM PYTHON:
// 219799.2508082828 1.017497627091681e-08
// 4.000000933333334e-11 4.440893134708784e-25
// /home/honza/Desktop/quadrature/test_integrals.py:29: IntegrationWarning: Extremely bad integrand behavior occurs at some points of the
//   integration interval.
//   integral, error = quad(func, *lim, epsrel=1e-13, limit=10000, epsabs=0)
// 26.514805364706493 3.936406756110955e-11
