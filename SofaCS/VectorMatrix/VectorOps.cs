using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// p-vector inner (=scalar=dot) product.
        /// </summary>
        /// <param name="a">double[3]     first p-vector</param>
        /// <param name="b">double[3]     second p-vector</param>
        /// <returns>a . b</returns>
        public static double iauPdp(in double[] a, in double[] b)
        {
            double w;

            w = a[0] * b[0]
               + a[1] * b[1]
               + a[2] * b[2];

            return w;
        }


        /// <summary>
        /// Modulus of p-vector.
        /// </summary>
        /// <param name="p">double[3]     p-vector</param>
        /// <returns>modulus</returns>
        public static double iauPm(in double[] p)
        {
            return Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        }


        /// <summary>
        /// P-vector subtraction.
        /// </summary>
        /// <param name="a">first p-vector</param>
        /// <param name="b">second p-vector</param>
        /// <param name="amb">double[3]      a - b</param>
        public static void iauPmp(in double[] a, in double[] b, ref double[] amb)
        {
            amb[0] = a[0] - b[0];
            amb[1] = a[1] - b[1];
            amb[2] = a[2] - b[2];
        }


        /// <summary>
        /// Convert a p-vector into modulus and unit vector.
        /// </summary>
        /// <param name="p">double[3]      p-vector</param>
        /// <param name="r">double         modulus</param>
        /// <param name="u">double[3]      unit vector</param>
        public static void iauPn(in double[] p, out double r, ref double[] u)
        {
            double w;

            /* Obtain the modulus and test for zero. */
            w = iauPm(p);
            if (w == 0.0)
            {
                /* Null vector. */
                iauZp(out u);
            }
            else
            {
                /* Unit vector. */
                iauSxp(1.0 / w, p, ref u);
            }

            /* Return the modulus. */
            r = w;
        }


        /// <summary>
        /// P-vector addition.
        /// </summary>
        /// <param name="a">double[3]      first p-vector</param>
        /// <param name="b">double[3]      second p-vector</param>
        /// <param name="apb">double[3]      a + b</param>
        public static void iauPpp(in double[] a, in double[] b, ref double[] apb)
        {
            apb[0] = a[0] + b[0];
            apb[1] = a[1] + b[1];
            apb[2] = a[2] + b[2];
        }


        /// <summary>
        /// P-vector plus scaled p-vector.
        /// </summary>
        /// <param name="a">double[3]     first p-vector</param>
        /// <param name="s">scalar (multiplier for b)</param>
        /// <param name="b">double[3]     second p-vector</param>
        /// <param name="apsb"></param>
        public static void iauPpsp(in double[] a, double s, double[] b, ref double[] apsb)
        {
            double[] sb = new double[3];

            /* s*b. */
            iauSxp(s, b, ref sb);

            /* a + s*b. */
            iauPpp(a, sb, ref apsb);
        }


        /// <summary>
        /// Inner (=scalar=dot) product of two pv-vectors.
        /// </summary>
        /// <param name="a">double[2][3]      first pv-vector</param>
        /// <param name="b">double[2][3]      second pv-vector</param>
        /// <param name="adb">double[2]         a . b (see note)</param>
        public static void iauPvdpv(in double[][] a, in double[][] b, out double[] adb)
        {
            double adbd, addb;

            adb = new double[2];

            /* a . b = constant part of result. */
            adb[0] = iauPdp(a[0], b[0]);

            /* a . bdot */
            adbd = iauPdp(a[0], b[1]);

            /* adot . b */
            addb = iauPdp(a[1], b[0]);

            /* Velocity part of result. */
            adb[1] = adbd + addb;
        }


        /// <summary>
        /// Modulus of pv-vector.
        /// </summary>
        /// <param name="pv">double[2][3]   pv-vector</param>
        /// <param name="r">modulus of position component</param>
        /// <param name="s">modulus of velocity component</param>
        public static void iauPvm(in double[][] pv, out double r, out double s)
        {
            /* Distance. */
            r = iauPm(pv[0]);

            /* Speed. */
            s = iauPm(pv[1]);
        }


        /// <summary>
        /// Subtract one pv-vector from another.
        /// </summary>
        /// <param name="a">double[2][3]      first pv-vector</param>
        /// <param name="b">double[2][3]      second pv-vector</param>
        /// <param name="amb">double[2][3]      a - b</param>
        public static void iauPvmpv(in double[][] a, in double[][] b, ref double[][] amb)
        {
            iauPmp(a[0], b[0], ref amb[0]);
            iauPmp(a[1], b[1], ref amb[1]);
        }


        /// <summary>
        /// Add one pv-vector to another.
        /// </summary>
        /// <param name="a">double[2][3]      first pv-vector</param>
        /// <param name="b">double[2][3]      second pv-vector</param>
        /// <param name="apb">double[2][3]      a + b</param>
        public static void iauPvppv(in double[][] a, in double[][] b, ref double[][] apb)
        {
            iauPpp(a[0], b[0], ref apb[0]);
            iauPpp(a[1], b[1], ref apb[1]);
        }


        /// <summary>
        /// Update a pv-vector.
        /// </summary>
        /// <param name="dt">time interval</param>
        /// <param name="pv">double[2][3]     pv-vector</param>
        /// <param name="upv">double[2][3]     p updated, v unchanged</param>
        public static void iauPvu(double dt, in double[][] pv, ref double[][] upv)
        {
            iauPpsp(pv[0], dt, pv[1], ref upv[0]);
            iauCp(pv[1], upv[1]);
        }


        /// <summary>
        /// Update a pv-vector, discarding the velocity component.
        /// </summary>
        /// <param name="dt">time interval</param>
        /// <param name="pv">double[2][3]      pv-vector</param>
        /// <param name="p">double[3]         p-vector</param>
        public static void iauPvup(double dt, in double[][] pv, out double[] p)
        {
            p = new double[3];

            p[0] = pv[0][0] + dt * pv[1][0];
            p[1] = pv[0][1] + dt * pv[1][1];
            p[2] = pv[0][2] + dt * pv[1][2];
        }


        /// <summary>
        /// Outer (=vector=cross) product of two pv-vectors.
        /// </summary>
        /// <param name="a">double[2][3]      first pv-vector</param>
        /// <param name="b">double[2][3]      second pv-vector</param>
        /// <param name="axb">double[2][3]      a x b</param>
        public static void iauPvxpv(in double[][] a, in double[][] b, ref double[][] axb)
        {
            double[] axbd = new double[3];
            double[] adxb = new double[3];
            double[][] wa = { new double[3], new double[3] };
            double[][] wb = { new double[3], new double[3] };

            /* Make copies of the inputs. */
            iauCpv(a, wa);
            iauCpv(b, wb);

            /* a x b = position part of result. */
            iauPxp(wa[0], wb[0], ref axb[0]);

            /* a x bdot + adot x b = velocity part of result. */
            iauPxp(wa[0], wb[1], ref axbd);
            iauPxp(wa[1], wb[0], ref adxb);
            iauPpp(axbd, adxb, ref axb[1]);
        }


        /// <summary>
        /// p-vector outer (=vector=cross) product.
        /// </summary>
        /// <param name="a">double[3]      first p-vector</param>
        /// <param name="b">double[3]      second p-vector</param>
        /// <param name="axb">double[3]      a x b</param>
        public static void iauPxp(in double[] a, in double[] b, ref double[] axb)
        {
            double xa, ya, za, xb, yb, zb;


            xa = a[0];
            ya = a[1];
            za = a[2];
            xb = b[0];
            yb = b[1];
            zb = b[2];
            axb[0] = ya * zb - za * yb;
            axb[1] = za * xb - xa * zb;
            axb[2] = xa * yb - ya * xb;
        }


        /// <summary>
        /// Multiply a pv-vector by two scalars.
        /// </summary>
        /// <param name="s1">scalar to multiply position component by</param>
        /// <param name="s2">scalar to multiply velocity component by</param>
        /// <param name="pv">double[2][3]   pv-vector</param>
        /// <param name="spv">double[2][3]   pv-vector: p scaled by s1, v scaled by s2</param>
        public static void iauS2xpv(double s1, double s2, in double[][] pv, ref double[][] spv)
        {
            iauSxp(s1, pv[0], ref spv[0]);
            iauSxp(s2, pv[1], ref spv[1]);
        }


        /// <summary>
        /// Multiply a p-vector by a scalar.
        /// </summary>
        /// <param name="s">scalar</param>
        /// <param name="p">double[3]     p-vector</param>
        /// <param name="sp">double[3]     s * p</param>
        public static void iauSxp(double s, in double[] p, ref double[] sp)
        {
            sp[0] = s * p[0];
            sp[1] = s * p[1];
            sp[2] = s * p[2];
        }


        /// <summary>
        /// Multiply a pv-vector by a scalar.
        /// </summary>
        /// <param name="s">scalar</param>
        /// <param name="pv">double[2][3]    pv-vector</param>
        /// <param name="spv">double[2][3]    s * pv</param>
        public static void iauSxpv(double s, in double[][] pv, ref double[][] spv)
        {
            iauS2xpv(s, s, pv, ref spv);
        }


    }
}
