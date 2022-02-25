using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Position-angle from two p-vectors.
        /// </summary>
        /// <param name="a">double[3]  direction of reference point</param>
        /// <param name="b">double[3]  direction of point whose PA is required</param>
        /// <returns>position angle of b with respect to a (radians)</returns>
        public static double iauPap(in double[] a, in double[] b)
        {
            double am, bm, st, ct, xa, ya, za, pa;
            double[] au = new double[3];
            double[] eta = new double[3];
            double[] xi = new double[3];
            double[] a2b = new double[3];


            /* Modulus and direction of the a vector. */
            iauPn(a, out am, ref au);

            /* Modulus of the b vector. */
            bm = iauPm(b);

            /* Deal with the case of a null vector. */
            if ((am == 0.0) || (bm == 0.0))
            {
                st = 0.0;
                ct = 1.0;
            }
            else
            {

                /* The "north" axis tangential from a (arbitrary length). */
                xa = a[0];
                ya = a[1];
                za = a[2];
                eta[0] = -xa * za;
                eta[1] = -ya * za;
                eta[2] = xa * xa + ya * ya;

                /* The "east" axis tangential from a (same length). */
                iauPxp(eta, au, ref xi);

                /* The vector from a to b. */
                iauPmp(b, a, ref a2b);

                /* Resolve into components along the north and east axes. */
                st = iauPdp(a2b, xi);
                ct = iauPdp(a2b, eta);

                /* Deal with degenerate cases. */
                if ((st == 0.0) && (ct == 0.0)) ct = 1.0;
            }

            /* Position angle. */
            pa = Math.Atan2(st, ct);

            return pa;
        }


        /// <summary>
        /// Position-angle from spherical coordinates.
        /// </summary>
        /// <param name="al">longitude of point A (e.g. RA) in radians</param>
        /// <param name="ap">latitude of point A (e.g. Dec) in radians</param>
        /// <param name="bl">longitude of point B</param>
        /// <param name="bp">latitude of point B</param>
        /// <returns>position angle of B with respect to A</returns>
        public static double iauPas(double al, double ap, double bl, double bp)
        {
            double dl, x, y, pa;


            dl = bl - al;
            y = Math.Sin(dl) * Math.Cos(bp);
            x = Math.Sin(bp) * Math.Cos(ap) - Math.Cos(bp) * Math.Sin(ap) * Math.Cos(dl);
            pa = ((x != 0.0) || (y != 0.0)) ? Math.Atan2(y, x) : 0.0;

            return pa;
        }


        /// <summary>
        /// Angular separation between two p-vectors.
        /// </summary>
        /// <param name="a">double[3]    first p-vector (not necessarily unit length)</param>
        /// <param name="b">double[3]    second p-vector (not necessarily unit length)</param>
        /// <returns>angular separation (radians, always positive)</returns>
        public static double iauSepp(in double[] a, in double[] b)
        {
            double[] axb = new double[3];
            double ss, cs, s;


            /* Sine of angle between the vectors, multiplied by the two moduli. */
            iauPxp(a, b, ref axb);
            ss = iauPm(axb);

            /* Cosine of the angle, multiplied by the two moduli. */
            cs = iauPdp(a, b);

            /* The angle. */
            s = ((ss != 0.0) || (cs != 0.0)) ? Math.Atan2(ss, cs) : 0.0;

            return s;
        }


        /// <summary>
        /// Angular separation between two sets of spherical coordinates.
        /// </summary>
        /// <param name="al">first longitude (radians)</param>
        /// <param name="ap">first latitude (radians)</param>
        /// <param name="bl">second longitude (radians)</param>
        /// <param name="bp">second latitude (radians)</param>
        /// <returns>angular separation (radians)</returns>
        public static double iauSeps(double al, double ap, double bl, double bp)
        {
            double[] ac, bc;
            double s;


            /* Spherical to Cartesian. */
            iauS2c(al, ap, out ac);
            iauS2c(bl, bp, out bc);

            /* Angle between the vectors. */
            s = iauSepp(ac, bc);

            return s;
        }


    }
}
