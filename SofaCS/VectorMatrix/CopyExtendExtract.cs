using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Copy a p-vector.
        /// </summary>
        /// <param name="p">double[3]     p-vector to be copied</param>
        /// <param name="c">double[3]     copy</param>
        /// TODO: out
        public static void iauCp(in double[] p, double[] c)
        {
            c[0] = p[0];
            c[1] = p[1];
            c[2] = p[2];
        }


        /// <summary>
        /// Copy a position/velocity vector.
        /// </summary>
        /// <param name="pv">double[2][3]    position/velocity vector to be copied</param>
        /// <param name="c">double[2][3]    copy</param>
        /// TODO: out
        public static void iauCpv(in double[][] pv, double[][] c)
        {
            iauCp(pv[0], c[0]);
            iauCp(pv[1], c[1]);
        }


        /// <summary>
        /// Copy an r-matrix.
        /// </summary>
        /// <param name="r">double[3][3]    r-matrix to be copied</param>
        /// <param name="c">double[3][3]    copy</param>
        /// TODO: out
        public static void iauCr(in double[][] r, double[][] c)
        {
            iauCp(r[0], c[0]);
            iauCp(r[1], c[1]);
            iauCp(r[2], c[2]);
        }


        /// <summary>
        /// Extend a p-vector to a pv-vector by appending a zero velocity.
        /// </summary>
        /// <param name="p">double[3]       p-vector</param>
        /// <param name="pv">double[2][3]    pv-vector</param>
        public static void iauP2pv(in double[] p, out double[][] pv)
        {
            pv = new double[][] { new double[3], new double[3] };
            iauCp(p, pv[0]);
            iauZp(out pv[1]);
        }


        /// <summary>
        /// Discard velocity component of a pv-vector.
        /// </summary>
        /// <param name="pv">double[2][3]     pv-vector</param>
        /// <param name="p">double[3]        p-vector</param>
        public static void iauPv2p(in double[][] pv, out double[] p)
        {
            p = new double[3];

            iauCp(pv[0], p);
        }


    }
}
