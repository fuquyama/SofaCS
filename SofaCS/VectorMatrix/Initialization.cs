using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Initialize an r-matrix to the identity matrix.
        /// </summary>
        /// <param name="r">double[3][3]    r-matrix</param>
        public static void iauIr(out double[][] r)
        {
            r = new double[][] { new double[3], new double[3], new double[3] };

            r[0][0] = 1.0;
            r[0][1] = 0.0;
            r[0][2] = 0.0;
            r[1][0] = 0.0;
            r[1][1] = 1.0;
            r[1][2] = 0.0;
            r[2][0] = 0.0;
            r[2][1] = 0.0;
            r[2][2] = 1.0;
        }


        /// <summary>
        /// Zero a p-vector.
        /// </summary>
        /// <param name="p">double[3]      zero p-vector</param>
        public static void iauZp(out double[] p)
        {
            p = new double[3];
            p[0] = 0.0;
            p[1] = 0.0;
            p[2] = 0.0;
        }


        /// <summary>
        /// Zero a pv-vector.
        /// </summary>
        /// <param name="pv">double[2][3]      zero pv-vector</param>
        public static void iauZpv(out double[][] pv)
        {
            pv = new double[][] { new double[3], new double[3] };
            iauZp(out pv[0]);
            iauZp(out pv[1]);
        }


        /// <summary>
        /// Initialize an r-matrix to the null matrix.
        /// </summary>
        /// <param name="r">double[3][3]    r-matrix</param>
        public static void iauZr(out double[][] r)
        {
            r = new double[][] { new double[3], new double[3], new double[3] };

            r[0][0] = 0.0;
            r[0][1] = 0.0;
            r[0][2] = 0.0;
            r[1][0] = 0.0;
            r[1][1] = 0.0;
            r[1][2] = 0.0;
            r[2][0] = 0.0;
            r[2][1] = 0.0;
            r[2][2] = 0.0;
        }


    }
}
