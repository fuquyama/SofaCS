using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Earth reference ellipsoids.
        /// </summary>
        /// <param name="n">ellipsoid identifier</param>
        /// <param name="a">equatorial radius (meters)</param>
        /// <param name="f">flattening</param>
        /// <returns>status:  0 = OK
        /// <para>-1 = illegal identifier</para>
        /// </returns>
        public static int iauEform(ReferenceEllipsoids n, out double a, out double f)
        {
            /* Look up a and f for the specified reference ellipsoid. */
            switch (n)
            {

                case ReferenceEllipsoids.WGS84:
                    a = 6378137.0;
                    f = 1.0 / 298.257223563;
                    break;

                case ReferenceEllipsoids.GRS80:
                    a = 6378137.0;
                    f = 1.0 / 298.257222101;
                    break;

                case ReferenceEllipsoids.WGS72:
                    a = 6378135.0;
                    f = 1.0 / 298.26;
                    break;

                default:

                    /* Invalid identifier. */
                    a = 0.0;
                    f = 0.0;
                    return -1;

            }

            /* OK status. */
            return 0;
        }


        /// <summary>
        /// Transform geocentric coordinates to geodetic using the specified
        /// reference ellipsoid.
        /// </summary>
        /// <param name="n">ellipsoid identifier</param>
        /// <param name="xyz">douvle[3] geocentric vector</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="height">height above ellipsoid (geodetic)</param>
        /// <returns>status:  0 = OK
        /// <para>-1 = illegal identifier</para>
        /// <para>-2 = internal error</para>
        /// </returns>
        public static int iauGc2gd(ReferenceEllipsoids n, double[] xyz,
               out double elong, out double phi, out double height)
        {
            int j;
            double a, f;

            elong = 0;
            phi = 0;
            height = 0;

            /* Obtain reference ellipsoid parameters. */
            j = iauEform(n, out a, out f);

            /* If OK, transform x,y,z to longitude, geodetic latitude, height. */
            if (j == 0)
            {
                j = iauGc2gde(a, f, xyz, out elong, out phi, out height);
                if (j < 0) j = -2;
            }

            /* Deal with any errors. */
            if (j < 0)
            {
                elong = -1e9;
                phi = -1e9;
                height = -1e9;
            }

            /* Return the status. */
            return j;
        }


        /// <summary>
        /// Transform geocentric coordinates to geodetic for a reference
        /// ellipsoid of specified form.
        /// </summary>
        /// <param name="a">equatorial radius</param>
        /// <param name="f">flattening</param>
        /// <param name="xyz">double[3]  geocentric vector</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="height">height above ellipsoid (geodetic)</param>
        /// <returns>status:  0 = OK
        /// <para>-1 = illegal f</para>
        /// <para>-2 = illegal a</para>
        /// </returns>
        public static int iauGc2gde(double a, double f, double[] xyz,
                out double elong, out double phi, out double height)
        {
            double aeps2, e2, e4t, ec2, ec, b, x, y, z, p2, absz, p, s0, pn, zc,
                          c0, c02, c03, s02, s03, a02, a0, a03, d0, f0, b0, s1,
                          cc, s12, cc2;

            elong = 0;
            phi = 0;
            height = 0;


            /* ------------- */
            /* Preliminaries */
            /* ------------- */

            /* Validate ellipsoid parameters. */
            if (f < 0.0 || f >= 1.0) return -1;
            if (a <= 0.0) return -2;

            /* Functions of ellipsoid parameters (with further validation of f). */
            aeps2 = a * a * 1e-32;
            e2 = (2.0 - f) * f;
            e4t = e2 * e2 * 1.5;
            ec2 = 1.0 - e2;
            if (ec2 <= 0.0) return -1;
            ec = Math.Sqrt(ec2);
            b = a * ec;

            /* Cartesian components. */
            x = xyz[0];
            y = xyz[1];
            z = xyz[2];

            /* Distance from polar axis squared. */
            p2 = x * x + y * y;

            /* Longitude. */
            elong = p2 > 0.0 ? Math.Atan2(y, x) : 0.0;

            /* Unsigned z-coordinate. */
            absz = Math.Abs(z);

            /* Proceed unless polar case. */
            if (p2 > aeps2)
            {

                /* Distance from polar axis. */
                p = Math.Sqrt(p2);

                /* Normalization. */
                s0 = absz / a;
                pn = p / a;
                zc = ec * s0;

                /* Prepare Newton correction factors. */
                c0 = ec * pn;
                c02 = c0 * c0;
                c03 = c02 * c0;
                s02 = s0 * s0;
                s03 = s02 * s0;
                a02 = c02 + s02;
                a0 = Math.Sqrt(a02);
                a03 = a02 * a0;
                d0 = zc * a03 + e2 * s03;
                f0 = pn * a03 - e2 * c03;

                /* Prepare Halley correction factor. */
                b0 = e4t * s02 * c02 * pn * (a0 - ec);
                s1 = d0 * f0 - b0 * s0;
                cc = ec * (f0 * f0 - b0 * c0);

                /* Evaluate latitude and height. */
                phi = Math.Atan(s1 / cc);
                s12 = s1 * s1;
                cc2 = cc * cc;
                height = (p * cc + absz * s1 - a * Math.Sqrt(ec2 * s12 + cc2)) /
                                                                  Math.Sqrt(s12 + cc2);
            }
            else
            {

                /* Exception: pole. */
                phi = DPI / 2.0;
                height = absz - b;
            }

            /* Restore sign of latitude. */
            if (z < 0) phi = -phi;

            /* OK status. */
            return 0;
        }


        /// <summary>
        /// Transform geodetic coordinates to geocentric using the specified
        /// reference ellipsoid.
        /// </summary>
        /// <param name="n">ellipsoid identifier</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="height">height above ellipsoid (geodetic)</param>
        /// <param name="xyz">double[3]  geocentric vector</param>
        /// <returns>status:  0 = OK
        /// <para>-1 = illegal identifier</para>
        /// <para>-2 = internal case</para>
        /// </returns>
        public static int iauGd2gc(ReferenceEllipsoids n, double elong, double phi, double height,
               out double[] xyz)
        {
            xyz = new double[3];

            int j;
            double a, f;


            /* Obtain reference ellipsoid parameters. */
            j = iauEform(n, out a, out f);

            /* If OK, transform longitude, geodetic latitude, height to x,y,z. */
            if (j == 0)
            {
                j = iauGd2gce(a, f, elong, phi, height, out xyz);
                if (j != 0) j = -2;
            }

            /* Deal with any errors. */
            if (j != 0) iauZp(out xyz);

            /* Return the status. */
            return j;
        }


        /// <summary>
        /// Transform geodetic coordinates to geocentric for a reference
        /// ellipsoid of specified form.
        /// </summary>
        /// <param name="a">equatorial radius</param>
        /// <param name="f">flattening</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="height">height above ellipsoid (geodetic)</param>
        /// <param name="xyz">double[3]  geocentric vector</param>
        /// <returns>status:  0 = OK
        /// <para>-1 = illegal case</para>
        /// </returns>
        public static int iauGd2gce(double a, double f, double elong, double phi,
                double height, out double[] xyz)
        {
            xyz = new double[3];

            double sp, cp, w, d, ac, As, r;


            /* Functions of geodetic latitude. */
            sp = Math.Sin(phi);
            cp = Math.Cos(phi);
            w = 1.0 - f;
            w = w * w;
            d = cp * cp + w * sp * sp;
            if (d <= 0.0) return -1;
            ac = a / Math.Sqrt(d);
            As = w * ac;

            /* Geocentric vector. */
            r = (ac + height) * cp;
            xyz[0] = r * Math.Cos(elong);
            xyz[1] = r * Math.Sin(elong);
            xyz[2] = (As + height) * sp;

            /* Success. */
            return 0;
        }


    }
}
