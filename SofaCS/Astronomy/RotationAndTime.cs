using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// The equation of the equinoxes, compatible with IAU 2000 resolutions,
        /// given the nutation in longitude and the mean obliquity.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="epsa">mean obliquity</param>
        /// <param name="dpsi">nutation in longitude</param>
        /// <returns>equation of the equinoxes</returns>
        public static double iauEe00(double date1, double date2, double epsa, double dpsi)
        {
            double ee;


            /* Equation of the equinoxes. */
            ee = dpsi * Math.Cos(epsa) + iauEect00(date1, date2);

            return ee;
        }


        /// <summary>
        /// Equation of the equinoxes, compatible with IAU 2000 resolutions.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>equation of the equinoxes</returns>
        public static double iauEe00a(double date1, double date2)
        {
            double dpsipr, depspr, epsa, dpsi, deps, ee;


            /* IAU 2000 precession-rate adjustments. */
            iauPr00(date1, date2, out dpsipr, out depspr);

            /* Mean obliquity, consistent with IAU 2000 precession-nutation. */
            epsa = iauObl80(date1, date2) + depspr;

            /* Nutation in longitude. */
            iauNut00a(date1, date2, out dpsi, out deps);

            /* Equation of the equinoxes. */
            ee = iauEe00(date1, date2, epsa, dpsi);

            return ee;
        }


        /// <summary>
        /// Equation of the equinoxes, compatible with IAU 2000 resolutions but
        /// using the truncated nutation model IAU 2000B.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>equation of the equinoxes</returns>
        public static double iauEe00b(double date1, double date2)
        {
            double dpsipr, depspr, epsa, dpsi, deps, ee;


            /* IAU 2000 precession-rate adjustments. */
            iauPr00(date1, date2, out dpsipr, out depspr);

            /* Mean obliquity, consistent with IAU 2000 precession-nutation. */
            epsa = iauObl80(date1, date2) + depspr;

            /* Nutation in longitude. */
            iauNut00b(date1, date2, out dpsi, out deps);

            /* Equation of the equinoxes. */
            ee = iauEe00(date1, date2, epsa, dpsi);

            return ee;
        }


        /// <summary>
        /// Equation of the equinoxes, compatible with IAU 2000 resolutions and
        /// IAU 2006/2000A precession-nutation.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>equation of the equinoxes</returns>
        public static double iauEe06a(double date1, double date2)
        {
            double gst06a, gmst06, ee;


            /* Apparent and mean sidereal times. */
            gst06a = iauGst06a(0.0, 0.0, date1, date2);
            gmst06 = iauGmst06(0.0, 0.0, date1, date2);

            /* Equation of the equinoxes. */
            ee = iauAnpm(gst06a - gmst06);

            return ee;
        }


        /* ----------------------------------------- */
        /* The series for the EE complementary terms */
        /* ----------------------------------------- */

        readonly struct TERM
        {
            public readonly int[] nfa;      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
            public readonly double s, c;     /* sine and cosine coefficients */

            public TERM(int[] nfa, double s, double c)
            {
                this.nfa = nfa;
                this.s = s;
                this.c = c;
            }
        };

        /* Terms of order t^0 */
        static readonly TERM[] e0 = {
            /* 1-10 */
            new TERM(new int[] { 0,  0,  0,  0,  1,  0,  0,  0}, 2640.96e-6, -0.39e-6 ),
            new TERM(new int[] { 0,  0,  0,  0,  2,  0,  0,  0},   63.52e-6, -0.02e-6 ),
            new TERM(new int[] { 0,  0,  2, -2,  3,  0,  0,  0},   11.75e-6,  0.01e-6 ),
            new TERM(new int[] { 0,  0,  2, -2,  1,  0,  0,  0},   11.21e-6,  0.01e-6 ),
            new TERM(new int[] { 0,  0,  2, -2,  2,  0,  0,  0},   -4.55e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  2,  0,  3,  0,  0,  0},    2.02e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  2,  0,  1,  0,  0,  0},    1.98e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  0,  0,  3,  0,  0,  0},   -1.72e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  1,  0,  0,  1,  0,  0,  0},   -1.41e-6, -0.01e-6 ),
            new TERM(new int[] { 0,  1,  0,  0, -1,  0,  0,  0},   -1.26e-6, -0.01e-6 ),

            /* 11-20 */
            new TERM(new int[] { 1,  0,  0,  0, -1,  0,  0,  0},   -0.63e-6,  0.00e-6 ),
            new TERM(new int[] { 1,  0,  0,  0,  1,  0,  0,  0},   -0.63e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  1,  2, -2,  3,  0,  0,  0},    0.46e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  1,  2, -2,  1,  0,  0,  0},    0.45e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  4, -4,  4,  0,  0,  0},    0.36e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  1, -1,  1, -8, 12,  0},   -0.24e-6, -0.12e-6 ),
            new TERM(new int[] { 0,  0,  2,  0,  0,  0,  0,  0},    0.32e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  2,  0,  2,  0,  0,  0},    0.28e-6,  0.00e-6 ),
            new TERM(new int[] { 1,  0,  2,  0,  3,  0,  0,  0},    0.27e-6,  0.00e-6 ),
            new TERM(new int[] { 1,  0,  2,  0,  1,  0,  0,  0},    0.26e-6,  0.00e-6 ),

            /* 21-30 */
            new TERM(new int[] { 0,  0,  2, -2,  0,  0,  0,  0},   -0.21e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  1, -2,  2, -3,  0,  0,  0},    0.19e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  1, -2,  2, -1,  0,  0,  0},    0.18e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  0,  0,  0,  8,-13, -1},   -0.10e-6,  0.05e-6 ),
            new TERM(new int[] { 0,  0,  0,  2,  0,  0,  0,  0},    0.15e-6,  0.00e-6 ),
            new TERM(new int[] { 2,  0, -2,  0, -1,  0,  0,  0},   -0.14e-6,  0.00e-6 ),
            new TERM(new int[] { 1,  0,  0, -2,  1,  0,  0,  0},    0.14e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  1,  2, -2,  2,  0,  0,  0},   -0.14e-6,  0.00e-6 ),
            new TERM(new int[] { 1,  0,  0, -2, -1,  0,  0,  0},    0.14e-6,  0.00e-6 ),
            new TERM(new int[] { 0,  0,  4, -2,  4,  0,  0,  0},    0.13e-6,  0.00e-6 ),

            /* 31-33 */
            new TERM(new int[] { 0,  0,  2, -2,  4,  0,  0,  0},   -0.11e-6,  0.00e-6 ),
            new TERM(new int[] { 1,  0, -2,  0, -3,  0,  0,  0},    0.11e-6,  0.00e-6 ),
            new TERM(new int[] { 1,  0, -2,  0, -1,  0,  0,  0},    0.11e-6,  0.00e-6 )
        };

        /* Terms of order t^1 */
        static TERM[] e1 = {
            new TERM(new int[] { 0,  0,  0,  0,  1,  0,  0,  0},    -0.87e-6,  0.00e-6 )
        };


        /// <summary>
        /// Equation of the equinoxes complementary terms, consistent with
        /// IAU 2000 resolutions.
        /// </summary>
        /// <param name="date1">TT</param>
        /// <param name="date2">TT</param>
        /// <returns>complementary terms</returns>
        public static double iauEect00(double date1, double date2)
        {
            /* Time since J2000.0, in Julian centuries */
            double t;

            /* Miscellaneous */
            int i, j;
            double a, s0, s1;

            /* Fundamental arguments */
            double[] fa = new double[14];

            /* Returned value. */
            double eect;


            /* Number of terms in the series */
            int NE0 = (int)(e0.Length);
            int NE1 = (int)(e1.Length);

            /* ------------------------------------------------------------------ */

            /* Interval between fundamental epoch J2000.0 and current date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Fundamental Arguments (from IERS Conventions 2003) */

            /* Mean anomaly of the Moon. */
            fa[0] = iauFal03(t);

            /* Mean anomaly of the Sun. */
            fa[1] = iauFalp03(t);

            /* Mean longitude of the Moon minus that of the ascending node. */
            fa[2] = iauFaf03(t);

            /* Mean elongation of the Moon from the Sun. */
            fa[3] = iauFad03(t);

            /* Mean longitude of the ascending node of the Moon. */
            fa[4] = iauFaom03(t);

            /* Mean longitude of Venus. */
            fa[5] = iauFave03(t);

            /* Mean longitude of Earth. */
            fa[6] = iauFae03(t);

            /* General precession in longitude. */
            fa[7] = iauFapa03(t);

            /* Evaluate the EE complementary terms. */
            s0 = 0.0;
            s1 = 0.0;

            for (i = NE0 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)(e0[i].nfa[j]) * fa[j];
                }
                s0 += e0[i].s * Math.Sin(a) + e0[i].c * Math.Cos(a);
            }

            for (i = NE1 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)(e1[i].nfa[j]) * fa[j];
                }
                s1 += e1[i].s * Math.Sin(a) + e1[i].c * Math.Cos(a);
            }

            eect = (s0 + s1 * t) * DAS2R;

            return eect;
        }


        /// <summary>
        /// Equation of the equinoxes, IAU 1994 model.
        /// </summary>
        /// <param name="date1">TDB</param>
        /// <param name="date2">TDB</param>
        /// <returns></returns>
        public static double iauEqeq94(double date1, double date2)
        {
            double t, om, dpsi, deps, eps0, ee;


            /* Interval between fundamental epoch J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Longitude of the mean ascending node of the lunar orbit on the */
            /* ecliptic, measured from the mean equinox of date. */
            om = iauAnpm((450160.280 + (-482890.539
                    + (7.455 + 0.008 * t) * t) * t) * DAS2R
                    + ((-5.0 * t) % 1.0) * D2PI);

            /* Nutation components and mean obliquity. */
            iauNut80(date1, date2, out dpsi, out deps);
            eps0 = iauObl80(date1, date2);

            /* Equation of the equinoxes. */
            ee = dpsi * Math.Cos(eps0) + DAS2R * (0.00264 * Math.Sin(om) + 0.000063 * Math.Sin(om + om));

            return ee;
        }


        /// <summary>
        /// Earth rotation angle (IAU 2000 model).
        /// </summary>
        /// <param name="dj1">UT1</param>
        /// <param name="dj2">UT1</param>
        /// <returns>Earth rotation angle (radians), range 0-2pi</returns>
        public static double iauEra00(double dj1, double dj2)
        {
            double d1, d2, t, f, theta;


            /* Days since fundamental epoch. */
            if (dj1 < dj2)
            {
                d1 = dj1;
                d2 = dj2;
            }
            else
            {
                d1 = dj2;
                d2 = dj1;
            }
            t = d1 + (d2 - DJ00);

            /* Fractional part of T (days). */
            f = (d1 % 1.0) + (d2 % 1.0);

            /* Earth rotation angle at this UT1. */
            theta = iauAnp(D2PI * (f + 0.7790572732640
                                     + 0.00273781191135448 * t));

            return theta;
        }


        /// <summary>
        /// Greenwich mean sidereal time (model consistent with IAU 2000
        /// resolutions).
        /// </summary>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <returns>Greenwich mean sidereal time (radians)</returns>
        public static double iauGmst00(double uta, double utb, double tta, double ttb)
        {
            double t, gmst;


            /* TT Julian centuries since J2000.0. */
            t = ((tta - DJ00) + ttb) / DJC;

            /* Greenwich Mean Sidereal Time, IAU 2000. */
            gmst = iauAnp(iauEra00(uta, utb) +
                            (0.014506 +
                            (4612.15739966 +
                            (1.39667721 +
                            (-0.00009344 +
                            (0.00001882)
                   * t) * t) * t) * t) * DAS2R);

            return gmst;
        }


        /// <summary>
        /// Greenwich mean sidereal time (consistent with IAU 2006 precession).
        /// </summary>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <returns>Greenwich mean sidereal time (radians)</returns>
        public static double iauGmst06(double uta, double utb, double tta, double ttb)
        {
            double t, gmst;


            /* TT Julian centuries since J2000.0. */
            t = ((tta - DJ00) + ttb) / DJC;

            /* Greenwich mean sidereal time, IAU 2006. */
            gmst = iauAnp(iauEra00(uta, utb) +
                           (0.014506 +
                           (4612.156534 +
                           (1.3915817 +
                           (-0.00000044 +
                           (-0.000029956 +
                           (-0.0000000368)
                   * t) * t) * t) * t) * t) * DAS2R);

            return gmst;
        }


        /// <summary>
        /// Universal Time to Greenwich mean sidereal time (IAU 1982 model).
        /// </summary>
        /// <param name="dj1">UT1 Julian Date</param>
        /// <param name="dj2">UT1 Julian Date</param>
        /// <returns>Greenwich mean sidereal time (radians)</returns>
        public static double iauGmst82(double dj1, double dj2)
        {
            /* Coefficients of IAU 1982 GMST-UT1 model */
            double A = 24110.54841 - DAYSEC / 2.0;
            double B = 8640184.812866;
            double C = 0.093104;
            double D = -6.2e-6;

            /* The first constant, A, has to be adjusted by 12 hours because the */
            /* UT1 is supplied as a Julian date, which begins at noon.           */

            double d1, d2, t, f, gmst;


            /* Julian centuries since fundamental epoch. */
            if (dj1 < dj2)
            {
                d1 = dj1;
                d2 = dj2;
            }
            else
            {
                d1 = dj2;
                d2 = dj1;
            }
            t = (d1 + (d2 - DJ00)) / DJC;

            /* Fractional part of JD(UT1), in seconds. */
            f = DAYSEC * ((d1 % 1.0) + (d2 % 1.0));

            /* GMST at this UT1. */
            gmst = iauAnp(DS2R * ((A + (B + (C + D * t) * t) * t) + f));

            return gmst;
        }


        /// <summary>
        /// Greenwich apparent sidereal time (consistent with IAU 2000
        /// resolutions).
        /// </summary>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <returns>Greenwich apparent sidereal time (radians)</returns>
        public static double iauGst00a(double uta, double utb, double tta, double ttb)
        {
            double gmst00, ee00a, gst;


            gmst00 = iauGmst00(uta, utb, tta, ttb);
            ee00a = iauEe00a(tta, ttb);
            gst = iauAnp(gmst00 + ee00a);

            return gst;
        }


        /// <summary>
        /// Greenwich apparent sidereal time (consistent with IAU 2000
        /// resolutions but using the truncated nutation model IAU 2000B).
        /// </summary>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <returns>Greenwich apparent sidereal time (radians)</returns>
        public static double iauGst00b(double uta, double utb)
        {
            double gmst00, ee00b, gst;


            gmst00 = iauGmst00(uta, utb, uta, utb);
            ee00b = iauEe00b(uta, utb);
            gst = iauAnp(gmst00 + ee00b);

            return gst;
        }


        /// <summary>
        /// Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.
        /// </summary>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <param name="rnpb">double[3][3]  nutation x precession x bias matrix</param>
        /// <returns>Greenwich apparent sidereal time (radians)</returns>
        public static double iauGst06(double uta, double utb, double tta, double ttb,
                in double[][] rnpb)
        {
            double x, y, s, era, eors, gst;


            /* Extract CIP coordinates. */
            iauBpn2xy(rnpb, out x, out y);

            /* The CIO locator, s. */
            s = iauS06(tta, ttb, x, y);

            /* Greenwich apparent sidereal time. */
            era = iauEra00(uta, utb);
            eors = iauEors(rnpb, s);
            gst = iauAnp(era - eors);

            return gst;
        }


        /// <summary>
        /// Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
        /// resolutions).
        /// </summary>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <returns>Greenwich apparent sidereal time (radians)</returns>
        public static double iauGst06a(double uta, double utb, double tta, double ttb)
        {
            double[][] rnpb;
            double gst;


            /* Classical nutation x precession x bias matrix, IAU 2000A. */
            iauPnm06a(tta, ttb, out rnpb);

            /* Greenwich apparent sidereal time. */
            gst = iauGst06(uta, utb, tta, ttb, rnpb);

            return gst;
        }


        /// <summary>
        /// Greenwich apparent sidereal time (consistent with IAU 1982/94
        /// resolutions).
        /// </summary>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <returns>Greenwich apparent sidereal time (radians)</returns>
        public static double iauGst94(double uta, double utb)
        {
            double gmst82, eqeq94, gst;


            gmst82 = iauGmst82(uta, utb);
            eqeq94 = iauEqeq94(uta, utb);
            gst = iauAnp(gmst82 + eqeq94);

            return gst;
        }


    }
}
