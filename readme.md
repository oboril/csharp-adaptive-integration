### Adaptive quadrature in C#

The file `numerical_integration.cs` contains a standalone general-purpose adaptive integrator based on Gauss-Kronrod (G7K15) quadrature.

The accuracy of the integrator was validated on various functions against the integrator in `scipy`.

For the test functions and example use see `Program.cs`.