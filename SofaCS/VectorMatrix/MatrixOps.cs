using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Multiply two r-matrices.
        /// </summary>
        /// <param name="a">double[3][3]    first r-matrix</param>
        /// <param name="b">double[3][3]    second r-matrix</param>
        /// <param name="atb">double[3][3]    a * b</param>
        public static void iauRxr(in double[][] a, in double[][] b, ref double[][] atb)
        {
            int i, j, k;
            double w;
            double[][] wm = { new double[3], new double[3], new double[3] };


            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    w = 0.0;
                    for (k = 0; k < 3; k++)
                    {
                        w += a[i][k] * b[k][j];
                    }
                    wm[i][j] = w;
                }
            }
            iauCr(wm, atb);
        }


        /// <summary>
        /// Transpose an r-matrix.
        /// </summary>
        /// <param name="r">double[3][3]    r-matrix</param>
        /// <param name="rt">double[3][3]    transpose</param>
        public static void iauTr(in double[][] r, ref double[][] rt)
        {
            double[][] wm = { new double[3], new double[3], new double[3] };
            int i, j;


            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    wm[i][j] = r[j][i];
                }
            }
            iauCr(wm, rt);
        }


    }
}
