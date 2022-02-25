using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Multiply a p-vector by an r-matrix.
        /// </summary>
        /// <param name="r">double[3][3]    r-matrix</param>
        /// <param name="p">double[3]       p-vector</param>
        /// <param name="rp">double[3]       r * p</param>
        public static void iauRxp(in double[][] r, in double[] p, ref double[] rp)
        {
            double w;
            double[] wrp = new double[3];
            int i, j;


            /* Matrix r * vector p. */
            for (j = 0; j < 3; j++)
            {
                w = 0.0;
                for (i = 0; i < 3; i++)
                {
                    w += r[j][i] * p[i];
                }
                wrp[j] = w;
            }

            /* Return the result. */
            iauCp(wrp, rp);
        }


        /// <summary>
        /// Multiply a pv-vector by an r-matrix.
        /// </summary>
        /// <param name="r">double[3][3]    r-matrix</param>
        /// <param name="pv">double[2][3]    pv-vector</param>
        /// <param name="rpv">double[2][3]    r * pv</param>
        public static void iauRxpv(in double[][] r, in double[][] pv, ref double[][] rpv)
        {
            iauRxp(r, pv[0], ref rpv[0]);
            iauRxp(r, pv[1], ref rpv[1]);
        }


        /// <summary>
        /// Multiply a p-vector by the transpose of an r-matrix.
        /// </summary>
        /// <param name="r">double[3][3]   r-matrix</param>
        /// <param name="p">double[3]      p-vector</param>
        /// <param name="trp">double[3]      r^T * p</param>
        public static void iauTrxp(in double[][] r, in double[] p, ref double[] trp)
        {
            double[][] tr = { new double[3], new double[3], new double[3] };


            /* Transpose of matrix r. */
            iauTr(r, ref tr);

            /* Matrix tr * vector p -> vector trp. */
            iauRxp(tr, p, ref trp);
        }


        /// <summary>
        /// Multiply a pv-vector by the transpose of an r-matrix.
        /// </summary>
        /// <param name="r">double[3][3]    r-matrix</param>
        /// <param name="pv">double[2][3]    pv-vector</param>
        /// <param name="trpv">double[2][3]    r^T * pv</param>
        public static void iauTrxpv(in double[][] r, in double[][] pv, ref double[][] trpv)
        {
            double[][] tr = { new double[3], new double[3], new double[3] };


            /* Transpose of matrix r. */
            iauTr(r, ref tr);

            /* Matrix tr * vector pv -> vector trpv. */
            iauRxpv(tr, pv, ref trpv);
        }


    }
}
