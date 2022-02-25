using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Rotate an r-matrix about the x-axis.
        /// </summary>
        /// <param name="phi">angle (radians)</param>
        /// <param name="r">double[3][3]    r-matrix, rotated</param>
        public static void iauRx(double phi, double[][] r)
        {
            double s, c, a10, a11, a12, a20, a21, a22;


            s = Math.Sin(phi);
            c = Math.Cos(phi);

            a10 = c * r[1][0] + s * r[2][0];
            a11 = c * r[1][1] + s * r[2][1];
            a12 = c * r[1][2] + s * r[2][2];
            a20 = -s * r[1][0] + c * r[2][0];
            a21 = -s * r[1][1] + c * r[2][1];
            a22 = -s * r[1][2] + c * r[2][2];

            r[1][0] = a10;
            r[1][1] = a11;
            r[1][2] = a12;
            r[2][0] = a20;
            r[2][1] = a21;
            r[2][2] = a22;
        }


        /// <summary>
        /// Rotate an r-matrix about the y-axis.
        /// </summary>
        /// <param name="theta">angle (radians)</param>
        /// <param name="r">double[3][3]    r-matrix, rotated</param>
        public static void iauRy(double theta, double[][] r)
        {
            double s, c, a00, a01, a02, a20, a21, a22;


            s = Math.Sin(theta);
            c = Math.Cos(theta);

            a00 = c * r[0][0] - s * r[2][0];
            a01 = c * r[0][1] - s * r[2][1];
            a02 = c * r[0][2] - s * r[2][2];
            a20 = s * r[0][0] + c * r[2][0];
            a21 = s * r[0][1] + c * r[2][1];
            a22 = s * r[0][2] + c * r[2][2];

            r[0][0] = a00;
            r[0][1] = a01;
            r[0][2] = a02;
            r[2][0] = a20;
            r[2][1] = a21;
            r[2][2] = a22;
        }


        /// <summary>
        /// Rotate an r-matrix about the z-axis.
        /// </summary>
        /// <param name="psi">angle (radians)</param>
        /// <param name="r">double[3][3]    r-matrix, rotated</param>
        public static void iauRz(double psi, double[][] r)
        {
            double s, c, a00, a01, a02, a10, a11, a12;


            s = Math.Sin(psi);
            c = Math.Cos(psi);

            a00 = c * r[0][0] + s * r[1][0];
            a01 = c * r[0][1] + s * r[1][1];
            a02 = c * r[0][2] + s * r[1][2];
            a10 = -s * r[0][0] + c * r[1][0];
            a11 = -s * r[0][1] + c * r[1][1];
            a12 = -s * r[0][2] + c * r[1][2];

            r[0][0] = a00;
            r[0][1] = a01;
            r[0][2] = a02;
            r[1][0] = a10;
            r[1][1] = a11;
            r[1][2] = a12;
        }


    }
}
