using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Transformation from ecliptic coordinates (mean equinox and ecliptic
        /// of date) to ICRS RA,Dec, using the IAU 2006 precession model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian date</param>
        /// <param name="date2">TT as a 2-part Julian date</param>
        /// <param name="dl">ecliptic longitude</param>
        /// <param name="db">ecliptic latitude</param>
        /// <param name="dr">ICRS right ascension (radians)</param>
        /// <param name="dd">ICRS declination (radians)</param>
        public static void iauEceq06(double date1, double date2, double dl, double db,
               out double dr, out double dd)
        {
            double a, b;
            double[][] rm;
            double[] v1;
            double[] v2 = new double[3];

            /* Spherical to Cartesian. */
            iauS2c(dl, db, out v1);

            /* Rotation matrix, ICRS equatorial to ecliptic. */
            iauEcm06(date1, date2, out rm);

            /* The transformation from ecliptic to ICRS. */
            iauTrxp(rm, v1, ref v2);

            /* Cartesian to spherical. */
            iauC2s(v2, out a, out b);

            /* Express in conventional ranges. */
            dr = iauAnp(a);
            dd = iauAnpm(b);
        }


        /// <summary>
        /// ICRS equatorial to ecliptic rotation matrix, IAU 2006.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian date</param>
        /// <param name="date2">TT as a 2-part Julian date</param>
        /// <param name="rm">double[3][3]   ICRS to ecliptic rotation matrix</param>
        public static void iauEcm06(double date1, double date2, out double[][] rm)
        {
            double ob;
            double[][] bp;
            double[][] e;

            rm = new double[][] { new double[3], new double[3], new double[3] };

            /* Obliquity, IAU 2006. */
            ob = iauObl06(date1, date2);

            /* Precession-bias matrix, IAU 2006. */
            iauPmat06(date1, date2, out bp);

            /* Equatorial of date to ecliptic matrix. */
            iauIr(out e);
            iauRx(ob, e);

            /* ICRS to ecliptic coordinates rotation matrix, IAU 2006. */
            iauRxr(e, bp, ref rm);
        }


        /// <summary>
        /// Transformation from ICRS equatorial coordinates to ecliptic
        /// coordinates (mean equinox and ecliptic of date) using IAU 2006
        /// precession model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian date</param>
        /// <param name="date2">TT as a 2-part Julian date</param>
        /// <param name="dr">ICRS right ascension (radians)</param>
        /// <param name="dd">ICRS right declination (radians)</param>
        /// <param name="dl">ecliptic longitude (radians)</param>
        /// <param name="db">ecliptic latitude (radians)</param>
        public static void iauEqec06(double date1, double date2, double dr, double dd,
               out double dl, out double db)
        {
            double a, b;
            double[][] rm = { new double[3], new double[3], new double[3] };
            double[] v1;
            double[] v2 = new double[3];

            /* Spherical to Cartesian. */
            iauS2c(dr, dd, out v1);

            /* Rotation matrix, ICRS equatorial to ecliptic. */
            iauEcm06(date1, date2, out rm);

            /* The transformation from ICRS to ecliptic. */
            iauRxp(rm, v1, ref v2);

            /* Cartesian to spherical. */
            iauC2s(v2, out a, out b);

            /* Express in conventional ranges. */
            dl = iauAnp(a);
            db = iauAnpm(b);
        }


        /// <summary>
        /// Transformation from ecliptic coordinates (mean equinox and ecliptic
        /// of date) to ICRS RA,Dec, using a long-term precession model.
        /// </summary>
        /// <param name="epj">Julian epoch (TT)</param>
        /// <param name="dl">ecliptic longitude (radians)</param>
        /// <param name="db">ecliptic latitude (radians)</param>
        /// <param name="dr">ICRS right ascension (radians)</param>
        /// <param name="dd">ICRS right declination (radians)</param>
        public static void iauLteceq(double epj, double dl, double db, out double dr, out double dd)
        {
            double a, b;
            double[][] rm = { new double[3], new double[3], new double[3] };
            double[] v1;
            double[] v2 = new double[3];


            /* Spherical to Cartesian. */
            iauS2c(dl, db, out v1);

            /* Rotation matrix, ICRS equatorial to ecliptic. */
            iauLtecm(epj, out rm);

            /* The transformation from ecliptic to ICRS. */
            iauTrxp(rm, v1, ref v2);

            /* Cartesian to spherical. */
            iauC2s(v2, out a, out b);

            /* Express in conventional ranges. */
            dr = iauAnp(a);
            dd = iauAnpm(b);
        }


        /// <summary>
        /// ICRS equatorial to ecliptic rotation matrix, long-term.
        /// </summary>
        /// <param name="epj">Julian epoch (TT)</param>
        /// <param name="rm">double[3][3]   ICRS to ecliptic rotation matrix</param>
        public static void iauLtecm(double epj, out double[][] rm)
        {
            /* Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33) */
            double dx = -0.016617 * DAS2R,
                         de = -0.0068192 * DAS2R,
                         dr = -0.0146 * DAS2R;

            double s;
            double[] p;
            double[] z;
            double[] w = new double[3];
            double[] x = new double[3];
            double[] y = new double[3];

            rm = new double[][] { new double[3], new double[3], new double[3] };

            /* Equator pole. */
            iauLtpequ(epj, out p);

            /* Ecliptic pole (bottom row of equatorial to ecliptic matrix). */
            iauLtpecl(epj, out z);

            /* Equinox (top row of matrix). */
            iauPxp(p, z, ref w);
            iauPn(w, out s, ref x);

            /* Middle row of matrix. */
            iauPxp(z, x, ref y);

            /* Combine with frame bias. */
            rm[0][0] = x[0] - x[1] * dr + x[2] * dx;
            rm[0][1] = x[0] * dr + x[1] + x[2] * de;
            rm[0][2] = -x[0] * dx - x[1] * de + x[2];
            rm[1][0] = y[0] - y[1] * dr + y[2] * dx;
            rm[1][1] = y[0] * dr + y[1] + y[2] * de;
            rm[1][2] = -y[0] * dx - y[1] * de + y[2];
            rm[2][0] = z[0] - z[1] * dr + z[2] * dx;
            rm[2][1] = z[0] * dr + z[1] + z[2] * de;
            rm[2][2] = -z[0] * dx - z[1] * de + z[2];
        }


        /// <summary>
        /// Transformation from ICRS equatorial coordinates to ecliptic
        /// coordinates (mean equinox and ecliptic of date) using a long-term
        /// precession model.
        /// </summary>
        /// <param name="epj">Julian epoch (TT)</param>
        /// <param name="dr">ICRS right ascension (radians)</param>
        /// <param name="dd">ICRS right declination (radians)</param>
        /// <param name="dl">ecliptic longitude (radians)</param>
        /// <param name="db">ecliptic latitude (radians)</param>
        public static void iauLteqec(double epj, double dr, double dd,
            out double dl, out double db)
        {
            double a, b;
            double[] v1;
            double[] v2 = new double[3];
            double[][] rm = { new double[3], new double[3], new double[3] };


            /* Spherical to Cartesian. */
            iauS2c(dr, dd, out v1);

            /* Rotation matrix, ICRS equatorial to ecliptic. */
            iauLtecm(epj, out rm);

            /* The transformation from ICRS to ecliptic. */
            iauRxp(rm, v1, ref v2);

            /* Cartesian to spherical. */
            iauC2s(v2, out a, out b);

            /* Express in conventional ranges. */
            dl = iauAnp(a);
            db = iauAnpm(b);
        }


    }
}
