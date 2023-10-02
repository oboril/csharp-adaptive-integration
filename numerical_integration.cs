using System.Runtime.InteropServices;
using System.Security.Cryptography.X509Certificates;

namespace NumericalIntegration
{
    public struct IntegrationResult
    {
        public double integral;
        public double error;

        public IntegrationResult(double _integral, double _error) 
        {
            integral = _integral;
            error = _error;
        }
    }

    struct IntegrationRegion
    {
        public double a;
        public double b;
        public double integral;
        public double error;

        public IntegrationRegion(double _a, double _b, double _integral, double _error)
        {
            a = _a;
            b = _b;
            integral = _integral;
            error = _error;
        }
    }

    public static class AdaptiveQuadrature
    {
        private static readonly double[] points = {
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
        private static readonly double[] gauss_weights = {
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
        private static readonly double[] kronrod_weights = {
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

        private static IntegrationRegion fixed_quadrature(Func<double, double> func, double a, double b)
        {
            double ba2 = (b-a)/2;

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

            double error = Math.Abs(kronrod-gauss);

            return new IntegrationRegion(a, b, kronrod, error);
        }

        /// <summary>
        /// Adaptive integration using Gauss-Kronrod quadrature (G7K15). Integration terminates when epsrel or epsabs is satisfied.
        /// </summary>
        /// <param name="func">Function to be integrated</param>
        /// <param name="a">SLower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="epsabs">Maximum absolute error</param>
        /// <param name="epsrel">Maximum relative error</param>
        /// <param name="max_step">Size of largest allowed step</param>
        /// <returns>Value and estimated error of the integration</returns>
        public static IntegrationResult integrate(
            Func<double, double> func,
            double a,
            double b,
            double epsabs,
            double epsrel,
            double max_step
        )
        {
            // The quadrature is 15-point, thus adjust minimal and maximal step
            max_step *= 15;

            PriorityQueue<IntegrationRegion, double> regions = new PriorityQueue<IntegrationRegion, double>();

            double integral = 0;
            double error = 0;

            // Create initial integration regions according to maximum step size
            double ba = b-a;
            int init_regions = (int)Math.Ceiling(ba/max_step);
            if (init_regions < 1) init_regions = 1;
            double init_region_size = ba/init_regions;
            double ai = a;
            double bi = ai + init_region_size;
            IntegrationRegion region;
            for (int i = 0; i < init_regions - 1; i++)
            {
                region = fixed_quadrature(func, ai, bi);
                regions.Enqueue(region, -region.error);

                integral += region.integral;
                error += region.error;

                (ai, bi) = (bi, bi+init_region_size);
            }
            region = fixed_quadrature(func, ai, b);
            regions.Enqueue(region, -region.error);
            integral += region.integral;
            error += region.error;


            // Check convergence and subdivide region with largest error
            while (error > epsabs && error / Math.Abs(integral) > epsrel)
            {
                IntegrationRegion old_region = regions.Dequeue();

                // Check that the region size is large enough for numerical accuracy
                if (old_region.b - old_region.a < Math.Max(Math.Abs(old_region.a), Math.Abs(old_region.b)) * 1e-14)
                {
                    Console.WriteLine("WARNING: The integration step is smaller that spacing between numbers!");
                    break;
                }
                double midpoint = (old_region.a + old_region.b)/2;
                IntegrationRegion r1 = fixed_quadrature(func, old_region.a, midpoint);
                IntegrationRegion r2 = fixed_quadrature(func, midpoint, old_region.b);
                regions.Enqueue(r1, -r1.error);
                regions.Enqueue(r2, -r2.error);

                integral += r1.integral + r2.integral - old_region.integral;
                error += r1.error + r2.error - old_region.error;
            }
            
            Console.WriteLine("Integration regions: " + regions.Count.ToString());

            // Use the regions to recalculate integral and error to avoid rounding error
            integral = 0;
            error = 0;
            IntegrationRegion final_region;
            while (regions.TryDequeue(out final_region, out _))
            {
                integral += final_region.integral;
                error += final_region.error;
            }


            return new IntegrationResult(integral, error);
        }
    }
}