using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Express an r-matrix as an r-vector.
        /// </summary>
        /// <param name="r">double[3][3]    rotation matrix</param>
        /// <param name="w">double[3]       rotation vector</param>
        public static void iauRm2v(in double[][] r, out double[] w)
        {
            w = new double[3];

            double x, y, z, s2, c2, phi, f;


            x = r[1][2] - r[2][1];
            y = r[2][0] - r[0][2];
            z = r[0][1] - r[1][0];
            s2 = Math.Sqrt(x * x + y * y + z * z);
            if (s2 > 0)
            {
                c2 = r[0][0] + r[1][1] + r[2][2] - 1.0;
                phi = Math.Atan2(s2, c2);
                f = phi / s2;
                w[0] = x * f;
                w[1] = y * f;
                w[2] = z * f;
            }
            else
            {
                w[0] = 0.0;
                w[1] = 0.0;
                w[2] = 0.0;
            }
        }


        /// <summary>
        /// Form the r-matrix corresponding to a given r-vector.
        /// </summary>
        /// <param name="w">double[3]      rotation vector</param>
        /// <param name="r">double[3][3]    rotation matrix</param>
        public static void iauRv2m(in double[] w, out double[][] r)
        {
            r = new double[][] { new double[3], new double[3], new double[3] };

            double x, y, z, phi, s, c, f;


            /* Euler angle (magnitude of rotation vector) and functions. */
            x = w[0];
            y = w[1];
            z = w[2];
            phi = Math.Sqrt(x * x + y * y + z * z);
            s = Math.Sin(phi);
            c = Math.Cos(phi);
            f = 1.0 - c;

            /* Euler axis (direction of rotation vector), perhaps null. */
            if (phi > 0.0)
            {
                x /= phi;
                y /= phi;
                z /= phi;
            }

            /* Form the rotation matrix. */
            r[0][0] = x * x * f + c;
            r[0][1] = x * y * f + z * s;
            r[0][2] = x * z * f - y * s;
            r[1][0] = y * x * f - z * s;
            r[1][1] = y * y * f + c;
            r[1][2] = y * z * f + x * s;
            r[2][0] = z * x * f + y * s;
            r[2][1] = z * y * f - x * s;
            r[2][2] = z * z * f + c;
        }


    }
}
