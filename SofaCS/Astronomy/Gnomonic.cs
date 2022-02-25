using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// In the tangent plane projection, given the rectangular coordinates
        /// of a star and its spherical coordinates, determine the spherical
        /// coordinates of the tangent point.
        /// </summary>
        /// <param name="xi">rectangular coordinates of star image</param>
        /// <param name="eta">rectangular coordinates of star image</param>
        /// <param name="a">star's spherical coordinates (radians)</param>
        /// <param name="b">star's spherical coordinates (radians)</param>
        /// <param name="a01">tangent point's spherical coordinates, Soln. 1 (radians)</param>
        /// <param name="b01">tangent point's spherical coordinates, Soln. 1 (radians)</param>
        /// <param name="a02">tangent point's spherical coordinates, Soln. 2 (radians)</param>
        /// <param name="b02">tangent point's spherical coordinates, Soln. 2 (radians)</param>
        /// <returns>number of solutions:
        /// <para>0 = no solutions returned</para>
        /// <para>1 = only the first solution is useful</para>
        /// <para>2 = both solutions are useful</para>
        /// </returns>
        public static int iauTpors(double xi, double eta, double a, double b,
             out double a01, out double b01, out double a02, out double b02)
        {
            double xi2, r, sb, cb, rsb, rcb, w2, w, s, c;

            a01 = 0; b01 = 0; a02 = 0; b02 = 0;

            xi2 = xi * xi;
            r = Math.Sqrt(1.0 + xi2 + eta * eta);
            sb = Math.Sin(b);
            cb = Math.Cos(b);
            rsb = r * sb;
            rcb = r * cb;
            w2 = rcb * rcb - xi2;
            if (w2 >= 0.0)
            {
                w = Math.Sqrt(w2);
                s = rsb - eta * w;
                c = rsb * eta + w;
                if (xi == 0.0 && w == 0.0) w = 1.0;
                a01 = iauAnp(a - Math.Atan2(xi, w));
                b01 = Math.Atan2(s, c);
                w = -w;
                s = rsb - eta * w;
                c = rsb * eta + w;
                a02 = iauAnp(a - Math.Atan2(xi, w));
                b02 = Math.Atan2(s, c);
                return (Math.Abs(rsb) < 1.0) ? 1 : 2;
            }
            else
            {
                return 0;
            }
        }


        /// <summary>
        /// In the tangent plane projection, given the rectangular coordinates
        /// of a star and its direction cosines, determine the direction
        /// cosines of the tangent point.
        /// </summary>
        /// <param name="xi">rectangular coordinates of star image</param>
        /// <param name="eta">rectangular coordinates of star image</param>
        /// <param name="v">double[3] star's direction cosines</param>
        /// <param name="v01">double[3] tangent point's direction cosines, Solution 1</param>
        /// <param name="v02">double[3] tangent point's direction cosines, Solution 2</param>
        /// <returns>number of solutions:
        /// <para>0 = no solutions returned</para>
        /// <para>1 = only the first solution is useful</para>
        /// <para>2 = both solutions are useful</para>
        /// </returns>
        public static int iauTporv(double xi, double eta, in double[] v,
             out double[] v01, out double[] v02)
        {
            double x, y, z, rxy2, xi2, eta2p1, r, rsb, rcb, w2, w, c;
            v01 = new double[3];
            v02 = new double[3];

            x = v[0];
            y = v[1];
            z = v[2];
            rxy2 = x * x + y * y;
            xi2 = xi * xi;
            eta2p1 = eta * eta + 1.0;
            r = Math.Sqrt(xi2 + eta2p1);
            rsb = r * z;
            rcb = r * Math.Sqrt(x * x + y * y);
            w2 = rcb * rcb - xi2;
            if (w2 > 0.0)
            {
                w = Math.Sqrt(w2);
                c = (rsb * eta + w) / (eta2p1 * Math.Sqrt(rxy2 * (w2 + xi2)));
                v01[0] = c * (x * w + y * xi);
                v01[1] = c * (y * w - x * xi);
                v01[2] = (rsb - eta * w) / eta2p1;
                w = -w;
                c = (rsb * eta + w) / (eta2p1 * Math.Sqrt(rxy2 * (w2 + xi2)));
                v02[0] = c * (x * w + y * xi);
                v02[1] = c * (y * w - x * xi);
                v02[2] = (rsb - eta * w) / eta2p1;
                return (Math.Abs(rsb) < 1.0) ? 1 : 2;
            }
            else
            {
                return 0;
            }
        }


        /// <summary>
        /// In the tangent plane projection, given the star's rectangular
        /// coordinates and the spherical coordinates of the tangent point,
        /// solve for the spherical coordinates of the star.
        /// </summary>
        /// <param name="xi">rectangular coordinates of star image</param>
        /// <param name="eta">rectangular coordinates of star image</param>
        /// <param name="a0">tangent point's spherical coordinates</param>
        /// <param name="b0">tangent point's spherical coordinates</param>
        /// <param name="a">star's spherical coordinates</param>
        /// <param name="b">star's spherical coordinates</param>
        public static void iauTpsts(double xi, double eta, double a0, double b0,
              out double a, out double b)
        {
            double sb0, cb0, d;

            sb0 = Math.Sin(b0);
            cb0 = Math.Cos(b0);
            d = cb0 - eta * sb0;
            a = iauAnp(Math.Atan2(xi, d) + a0);
            b = Math.Atan2(sb0 + eta * cb0, Math.Sqrt(xi * xi + d * d));
        }


        /// <summary>
        /// In the tangent plane projection, given the star's rectangular
        /// coordinates and the direction cosines of the tangent point, solve
        /// for the direction cosines of the star.
        /// </summary>
        /// <param name="xi">rectangular coordinates of star image</param>
        /// <param name="eta">rectangular coordinates of star image</param>
        /// <param name="v0">double[3]  tangent point's direction cosines</param>
        /// <param name="v">double[3]  star's direction cosines</param>
        public static void iauTpstv(double xi, double eta, in double[] v0, out double[] v)
        {
            double x, y, z, f, r;

            v = new double[3];

            /* Tangent point. */
            x = v0[0];
            y = v0[1];
            z = v0[2];

            /* Deal with polar case. */
            r = Math.Sqrt(x * x + y * y);
            if (r == 0.0)
            {
                r = 1e-20;
                x = r;
            }

            /* Star vector length to tangent plane. */
            f = Math.Sqrt(1.0 + xi * xi + eta * eta);

            /* Apply the transformation and normalize. */
            v[0] = (x - (xi * y + eta * x * z) / r) / f;
            v[1] = (y + (xi * x - eta * y * z) / r) / f;
            v[2] = (z + eta * r) / f;
        }


        /// <summary>
        /// In the tangent plane projection, given celestial spherical
        /// coordinates for a star and the tangent point, solve for the star's
        /// rectangular coordinates in the tangent plane.
        /// </summary>
        /// <param name="a">star's spherical coordinates</param>
        /// <param name="b">star's spherical coordinates</param>
        /// <param name="a0">tangent point's spherical coordinates</param>
        /// <param name="b0">tangent point's spherical coordinates</param>
        /// <param name="xi">rectangular coordinates of star image</param>
        /// <param name="eta">rectangular coordinates of star image</param>
        /// <returns>status:
        /// <para>0 = OK</para>
        /// <para>1 = star too far from axis</para>
        /// <para>2 = antistar on tangent plane</para>
        /// <para>3 = antistar too far from axis</para>
        /// </returns>
        public static int iauTpxes(double a, double b, double a0, double b0,
             out double xi, out double eta)
        {
            double TINY = 1e-6;
            int j;
            double sb0, sb, cb0, cb, da, sda, cda, d;


            /* Functions of the spherical coordinates. */
            sb0 = Math.Sin(b0);
            sb = Math.Sin(b);
            cb0 = Math.Cos(b0);
            cb = Math.Cos(b);
            da = a - a0;
            sda = Math.Sin(da);
            cda = Math.Cos(da);

            /* Reciprocal of star vector length to tangent plane. */
            d = sb * sb0 + cb * cb0 * cda;

            /* Check for error cases. */
            if (d > TINY)
            {
                j = 0;
            }
            else if (d >= 0.0)
            {
                j = 1;
                d = TINY;
            }
            else if (d > -TINY)
            {
                j = 2;
                d = -TINY;
            }
            else
            {
                j = 3;
            }

            /* Return the tangent plane coordinates (even in dubious cases). */
            xi = cb * sda / d;
            eta = (sb * cb0 - cb * sb0 * cda) / d;

            /* Return the status. */
            return j;
        }


        /// <summary>
        /// In the tangent plane projection, given celestial direction cosines
        /// for a star and the tangent point, solve for the star's rectangular
        /// coordinates in the tangent plane.
        /// </summary>
        /// <param name="v">double[3]  direction cosines of star</param>
        /// <param name="v0">double[3]  direction cosines of tangent point</param>
        /// <param name="xi">tangent plane coordinates of star</param>
        /// <param name="eta">tangent plane coordinates of star</param>
        /// <returns>status:
        /// <para>0 = OK</para>
        /// <para>1 = star too far from axis</para>
        /// <para>2 = antistar on tangent plane</para>
        /// <para>3 = antistar too far from axis</para>
        /// </returns>
        public static int iauTpxev(in double[] v, in double[] v0, out double xi, out double eta)
        {
            double TINY = 1e-6;
            int j;
            double x, y, z, x0, y0, z0, r2, r, w, d;


            /* Star and tangent point. */
            x = v[0];
            y = v[1];
            z = v[2];
            x0 = v0[0];
            y0 = v0[1];
            z0 = v0[2];

            /* Deal with polar case. */
            r2 = x0 * x0 + y0 * y0;
            r = Math.Sqrt(r2);
            if (r == 0.0)
            {
                r = 1e-20;
                x0 = r;
            }

            /* Reciprocal of star vector length to tangent plane. */
            w = x * x0 + y * y0;
            d = w + z * z0;

            /* Check for error cases. */
            if (d > TINY)
            {
                j = 0;
            }
            else if (d >= 0.0)
            {
                j = 1;
                d = TINY;
            }
            else if (d > -TINY)
            {
                j = 2;
                d = -TINY;
            }
            else
            {
                j = 3;
            }

            /* Return the tangent plane coordinates (even in dubious cases). */
            d *= r;
            xi = (y * x0 - x * y0) / d;
            eta = (z * r2 - z0 * w) / d;

            /* Return the status. */
            return j;
        }


    }
}
