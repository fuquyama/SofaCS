using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Frame bias components of IAU 2000 precession-nutation models;  part
        /// of the Mathews-Herring-Buffett (MHB2000) nutation series, with
        /// additions.
        /// </summary>
        /// <param name="dpsibi">longitude corrections</param>
        /// <param name="depsbi">obliquity corrections</param>
        /// <param name="dra">the ICRS RA of the J2000.0 mean equinox</param>
        public static void iauBi00(out double dpsibi, out double depsbi, out double dra)
        {
            /* The frame bias corrections in longitude and obliquity */
            double DPBIAS = -0.041775 * DAS2R,
                         DEBIAS = -0.0068192 * DAS2R;

            /* The ICRS RA of the J2000.0 equinox (Chapront et al., 2002) */
            double DRA0 = -0.0146 * DAS2R;


            /* Return the results (which are fixed). */
            dpsibi = DPBIAS;
            depsbi = DEBIAS;
            dra = DRA0;
        }


        /// <summary>
        /// Frame bias and precession, IAU 2000.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rb">double[3][3]   frame bias matrix</param>
        /// <param name="rp">double[3][3]   precession matrix</param>
        /// <param name="rbp">double[3][3]   bias-precession matrix</param>
        public static void iauBp00(double date1, double date2,
             out double[][] rb, out double[][] rp, out double[][] rbp)
        {
            /* J2000.0 obliquity (Lieske et al. 1977) */
            double EPS0 = 84381.448 * DAS2R;

            double t, dpsibi, depsbi, dra0, psia77, oma77, chia,
                   dpsipr, depspr, psia, oma;
            double[][] rbw;

            rb = new double[][] { new double[3], new double[3], new double[3] };
            rbp = new double[][] { new double[3], new double[3], new double[3] };


            /* Interval between fundamental epoch J2000.0 and current date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Frame bias. */
            iauBi00(out dpsibi, out depsbi, out dra0);

            /* Precession angles (Lieske et al. 1977) */
            psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * DAS2R;
            oma77 = EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * DAS2R;
            chia = (10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * DAS2R;

            /* Apply IAU 2000 precession corrections. */
            iauPr00(date1, date2, out dpsipr, out depspr);
            psia = psia77 + dpsipr;
            oma = oma77 + depspr;

            /* Frame bias matrix: GCRS to J2000.0. */
            iauIr(out rbw);
            iauRz(dra0, rbw);
            iauRy(dpsibi * Math.Sin(EPS0), rbw);
            iauRx(-depsbi, rbw);
            iauCr(rbw, rb);

            /* Precession matrix: J2000.0 to mean of date. */
            iauIr(out rp);
            iauRx(EPS0, rp);
            iauRz(-psia, rp);
            iauRx(-oma, rp);
            iauRz(chia, rp);

            /* Bias-precession matrix: GCRS to mean of date. */
            iauRxr(rp, rbw, ref rbp);
        }


        /// <summary>
        /// Frame bias and precession, IAU 2006.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rb">double[3][3]   frame bias matrix</param>
        /// <param name="rp">double[3][3]   precession matrix</param>
        /// <param name="rbp">double[3][3]   bias-precession matrix</param>
        public static void iauBp06(double date1, double date2,
             out double[][] rb, out double[][] rp, out double[][] rbp)
        {
            double gamb, phib, psib, epsa;
            double[][] rbpw;
            double[][] rbt = { new double[3], new double[3], new double[3] };

            rp = new double[][] { new double[3], new double[3], new double[3] };
            rbp = new double[][] { new double[3], new double[3], new double[3] };


            /* B matrix. */
            iauPfw06(DJM0, DJM00, out gamb, out phib, out psib, out epsa);
            iauFw2m(gamb, phib, psib, epsa, out rb);

            /* PxB matrix (temporary). */
            iauPmat06(date1, date2, out rbpw);

            /* P matrix. */
            iauTr(rb, ref rbt);
            iauRxr(rbpw, rbt, ref rp);

            /* PxB matrix. */
            iauCr(rbpw, rbp);
        }


        /// <summary>
        /// Extract from the bias-precession-nutation matrix the X,Y coordinates
        /// of the Celestial Intermediate Pole.
        /// </summary>
        /// <param name="rbpn">double[3][3]  celestial-to-true matrix</param>
        /// <param name="x">Celestial Intermediate Pole</param>
        /// <param name="y">Celestial Intermediate Pole</param>
        public static void iauBpn2xy(in double[][] rbpn, out double x, out double y)
        {
            /* Extract the X,Y coordinates. */
            x = rbpn[2][0];
            y = rbpn[2][1];
        }


        /// <summary>
        /// Form the celestial-to-intermediate matrix for a given date using the
        /// IAU 2000A precession-nutation model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rc2i">double[3][3] celestial-to-intermediate matrix</param>
        public static void iauC2i00a(double date1, double date2, out double[][] rc2i)
        {
            double[][] rbpn;

            /* Obtain the celestial-to-true matrix (IAU 2000A). */
            iauPnm00a(date1, date2, out rbpn);

            /* Form the celestial-to-intermediate matrix. */
            iauC2ibpn(date1, date2, rbpn, out rc2i);
        }


        /// <summary>
        /// Form the celestial-to-intermediate matrix for a given date using the
        /// IAU 2000B precession-nutation model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rc2i">double[3][3] celestial-to-intermediate matrix</param>
        public static void iauC2i00b(double date1, double date2, out double[][] rc2i)
        {
            double[][] rbpn;

            /* Obtain the celestial-to-true matrix (IAU 2000B). */
            iauPnm00b(date1, date2, out rbpn);

            /* Form the celestial-to-intermediate matrix. */
            iauC2ibpn(date1, date2, rbpn, out rc2i);
        }


        /// <summary>
        /// Form the celestial-to-intermediate matrix for a given date using the
        /// IAU 2006 precession and IAU 2000A nutation models.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rc2i">double[3][3] celestial-to-intermediate matrix</param>
        public static void iauC2i06a(double date1, double date2, out double[][] rc2i)
        {
            double x, y, s;
            double[][] rbpn;

            /* Obtain the celestial-to-true matrix (IAU 2006/2000A). */
            iauPnm06a(date1, date2, out rbpn);

            /* Extract the X,Y coordinates. */
            iauBpn2xy(rbpn, out x, out y);

            /* Obtain the CIO locator. */
            s = iauS06(date1, date2, x, y);

            /* Form the celestial-to-intermediate matrix. */
            iauC2ixys(x, y, s, out rc2i);
        }


        /// <summary>
        /// Form the celestial-to-intermediate matrix for a given date given
        /// the bias-precession-nutation matrix.  IAU 2000.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rbpn">double[3][3] celestial-to-true matrix</param>
        /// <param name="rc2i">double[3][3] celestial-to-intermediate matrix</param>
        public static void iauC2ibpn(double date1, double date2, in double[][] rbpn,
            out double[][] rc2i)
        {
            double x, y;


            /* Extract the X,Y coordinates. */
            iauBpn2xy(rbpn, out x, out y);

            /* Form the celestial-to-intermediate matrix (n.b. IAU 2000 specific). */
            iauC2ixy(date1, date2, x, y, out rc2i);
        }


        /// <summary>
        /// Form the celestial to intermediate-frame-of-date matrix for a given
        /// date when the CIP X,Y coordinates are known.  IAU 2000.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="x">Celestial Intermediate Pole</param>
        /// <param name="y">Celestial Intermediate Pole</param>
        /// <param name="rc2i">double[3][3] celestial-to-intermediate matrix</param>
        public static void iauC2ixy(double date1, double date2, double x, double y,
            out double[][] rc2i)
        {
            /* Compute s and then the matrix. */
            iauC2ixys(x, y, iauS00(date1, date2, x, y), out rc2i);
        }


        /// <summary>
        /// Form the celestial to intermediate-frame-of-date matrix given the CIP
        /// X,Y and the CIO locator s.
        /// </summary>
        /// <param name="x">Celestial Intermediate Pole</param>
        /// <param name="y">Celestial Intermediate Pole</param>
        /// <param name="s">the CIO locator s</param>
        /// <param name="rc2i">double[3][3]   celestial-to-intermediate matrix</param>
        public static void iauC2ixys(double x, double y, double s, out double[][] rc2i)
        {
            double r2, e, d;


            /* Obtain the spherical angles E and d. */
            r2 = x * x + y * y;
            e = (r2 > 0.0) ? Math.Atan2(y, x) : 0.0;
            d = Math.Atan(Math.Sqrt(r2 / (1.0 - r2)));

            /* Form the matrix. */
            iauIr(out rc2i);
            iauRz(e, rc2i);
            iauRy(d, rc2i);
            iauRz(-(e + s), rc2i);
        }


        /// <summary>
        /// Form the celestial to terrestrial matrix given the date, the UT1 and
        /// the polar motion, using the IAU 2000A precession-nutation model.
        /// </summary>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="xp">CIP coordinates (radians)</param>
        /// <param name="yp">CIP coordinates (radians)</param>
        /// <param name="rc2t">double[3][3]   celestial-to-terrestrial matrix</param>
        public static void iauC2t00a(double tta, double ttb, double uta, double utb,
               double xp, double yp, out double[][] rc2t)
        {
            double era, sp;
            double[][] rc2i;
            double[][] rpom;


            /* Form the celestial-to-intermediate matrix for this TT (IAU 2000A). */
            iauC2i00a(tta, ttb, out rc2i);

            /* Predict the Earth rotation angle for this UT1. */
            era = iauEra00(uta, utb);

            /* Estimate s'. */
            sp = iauSp00(tta, ttb);

            /* Form the polar motion matrix. */
            iauPom00(xp, yp, sp, out rpom);

            /* Combine to form the celestial-to-terrestrial matrix. */
            iauC2tcio(rc2i, era, rpom, out rc2t);
        }


        /// <summary>
        /// Form the celestial to terrestrial matrix given the date, the UT1 and
        /// the polar motion, using the IAU 2000B precession-nutation model.
        /// </summary>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="xp">coordinates of the pole (radians)</param>
        /// <param name="yp">coordinates of the pole (radians)</param>
        /// <param name="rc2t">double[3][3]   celestial-to-terrestrial matrix</param>
        public static void iauC2t00b(double tta, double ttb, double uta, double utb,
               double xp, double yp, out double[][] rc2t)
        {
            double era;
            double[][] rc2i;
            double[][] rpom;


            /* Form the celestial-to-intermediate matrix for this TT (IAU 2000B). */
            iauC2i00b(tta, ttb, out rc2i);

            /* Predict the Earth rotation angle for this UT1. */
            era = iauEra00(uta, utb);

            /* Form the polar motion matrix (neglecting s'). */
            iauPom00(xp, yp, 0.0, out rpom);

            /* Combine to form the celestial-to-terrestrial matrix. */
            iauC2tcio(rc2i, era, rpom, out rc2t);
        }


        /// <summary>
        /// Form the celestial to terrestrial matrix given the date, the UT1 and
        /// the polar motion, using the IAU 2006/2000A precession-nutation model.
        /// </summary>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="xp">coordinates of the pole (radians)</param>
        /// <param name="yp">coordinates of the pole (radians)</param>
        /// <param name="rc2t">double[3][3]   celestial-to-terrestrial matrix</param>
        public static void iauC2t06a(double tta, double ttb, double uta, double utb,
               double xp, double yp, out double[][] rc2t)
        {
            double era, sp;
            double[][] rc2i;
            double[][] rpom;


            /* Form the celestial-to-intermediate matrix for this TT. */
            iauC2i06a(tta, ttb, out rc2i);

            /* Predict the Earth rotation angle for this UT1. */
            era = iauEra00(uta, utb);

            /* Estimate s'. */
            sp = iauSp00(tta, ttb);

            /* Form the polar motion matrix. */
            iauPom00(xp, yp, sp, out rpom);

            /* Combine to form the celestial-to-terrestrial matrix. */
            iauC2tcio(rc2i, era, rpom, out rc2t);
        }


        /// <summary>
        /// Assemble the celestial to terrestrial matrix from CIO-based
        /// components (the celestial-to-intermediate matrix, the Earth Rotation
        /// Angle and the polar motion matrix).
        /// </summary>
        /// <param name="rc2i">double[3][3]    celestial-to-intermediate matrix</param>
        /// <param name="era">Earth rotation angle (radians)</param>
        /// <param name="rpom">double[3][3]    polar-motion matrix</param>
        /// <param name="rc2t">double[3][3]    celestial-to-terrestrial matrix</param>
        public static void iauC2tcio(in double[][] rc2i, double era, in double[][] rpom,
            out double[][] rc2t)
        {
            double[][] r = { new double[3], new double[3], new double[3] };

            rc2t = new double[][] { new double[3], new double[3], new double[3] };

            /* Construct the matrix. */
            iauCr(rc2i, r);
            iauRz(era, r);
            iauRxr(rpom, r, ref rc2t);
        }


        /// <summary>
        /// Assemble the celestial to terrestrial matrix from equinox-based
        /// components (the celestial-to-true matrix, the Greenwich Apparent
        /// Sidereal Time and the polar motion matrix).
        /// </summary>
        /// <param name="rbpn">double[3][3]  celestial-to-true matrix</param>
        /// <param name="gst">Greenwich (apparent) Sidereal Time (radians)</param>
        /// <param name="rpom">double[3][3]  polar-motion matrix</param>
        /// <param name="rc2t">double[3][3]  celestial-to-terrestrial matrix</param>
        public static void iauC2teqx(in double[][] rbpn, double gst, in double[][] rpom,
            out double[][] rc2t)
        {
            double[][] r = { new double[3], new double[3], new double[3] };

            rc2t = new double[][]{ new double[3], new double[3], new double[3] };

            /* Construct the matrix. */
            iauCr(rbpn, r);
            iauRz(gst, r);
            iauRxr(rpom, r, ref rc2t);
        }


        /// <summary>
        /// Form the celestial to terrestrial matrix given the date, the UT1,
        /// the nutation and the polar motion.  IAU 2000.
        /// </summary>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation</param>
        /// <param name="deps">nutation</param>
        /// <param name="xp">coordinates of the pole (radians)</param>
        /// <param name="yp">coordinates of the pole (radians)</param>
        /// <param name="rc2t">double[3][3]  celestial-to-terrestrial matrix</param>
        public static void iauC2tpe(double tta, double ttb, double uta, double utb,
              double dpsi, double deps, double xp, double yp,
              out double[][] rc2t)
        {
            double epsa, gmst, ee, sp;
            double[][] rb, rp, rbp, rn, rbpn, rpom;

            /* Form the celestial-to-true matrix for this TT. */
            iauPn00(tta, ttb, dpsi, deps, out epsa, out rb, out rp, out rbp, out rn, out rbpn);

            /* Predict the Greenwich Mean Sidereal Time for this UT1 and TT. */
            gmst = iauGmst00(uta, utb, tta, ttb);

            /* Predict the equation of the equinoxes given TT and nutation. */
            ee = iauEe00(tta, ttb, epsa, dpsi);

            /* Estimate s'. */
            sp = iauSp00(tta, ttb);

            /* Form the polar motion matrix. */
            iauPom00(xp, yp, sp, out rpom);

            /* Combine to form the celestial-to-terrestrial matrix. */
            iauC2teqx(rbpn, gmst + ee, rpom, out rc2t);
        }


        /// <summary>
        /// Form the celestial to terrestrial matrix given the date, the UT1,
        /// the CIP coordinates and the polar motion.  IAU 2000.
        /// </summary>
        /// <param name="tta">TT as a 2-part Julian Date</param>
        /// <param name="ttb">TT as a 2-part Julian Date</param>
        /// <param name="uta">UT1 as a 2-part Julian Date</param>
        /// <param name="utb">UT1 as a 2-part Julian Date</param>
        /// <param name="x">Celestial Intermediate Pole</param>
        /// <param name="y">Celestial Intermediate Pole</param>
        /// <param name="xp">coordinates of the pole (radians)</param>
        /// <param name="yp">coordinates of the pole (radians)</param>
        /// <param name="rc2t">double[3][3]   celestial-to-terrestrial matrix</param>
        public static void iauC2txy(double tta, double ttb, double uta, double utb,
              double x, double y, double xp, double yp,
              out double[][] rc2t)
        {
            double era, sp;
            double[][] rc2i, rpom;


            /* Form the celestial-to-intermediate matrix for this TT. */
            iauC2ixy(tta, ttb, x, y, out rc2i);

            /* Predict the Earth rotation angle for this UT1. */
            era = iauEra00(uta, utb);

            /* Estimate s'. */
            sp = iauSp00(tta, ttb);

            /* Form the polar motion matrix. */
            iauPom00(xp, yp, sp, out rpom);

            /* Combine to form the celestial-to-terrestrial matrix. */
            iauC2tcio(rc2i, era, rpom, out rc2t);
        }


        /// <summary>
        /// Equation of the origins, IAU 2006 precession and IAU 2000A nutation.
        /// </summary>
        /// <param name="date1">TT</param>
        /// <param name="date2">TT</param>
        /// <returns>the equation of the origins in radians</returns>
        public static double iauEo06a(double date1, double date2)
        {
            double x, y, s, eo;
            double[][] r;

            /* Classical nutation x precession x bias matrix. */
            iauPnm06a(date1, date2, out r);

            /* Extract CIP coordinates. */
            iauBpn2xy(r, out x, out y);

            /* The CIO locator, s. */
            s = iauS06(date1, date2, x, y);

            /* Solve for the EO. */
            eo = iauEors(r, s);

            return eo;
        }


        /// <summary>
        /// Equation of the origins, given the classical NPB matrix and the quantity s.
        /// </summary>
        /// <param name="rnpb">double[3][3]  classical nutation x precession x bias matrix</param>
        /// <param name="s">the quantity s (the CIO locator) in radians</param>
        /// <returns>the equation of the origins in radians</returns>
        public static double iauEors(in double[][] rnpb, double s)
        {
            double x, ax, xs, ys, zs, p, q, eo;


            /* Evaluate Wallace & Capitaine (2006) expression (16). */
            x = rnpb[2][0];
            ax = x / (1.0 + rnpb[2][2]);
            xs = 1.0 - ax * x;
            ys = -ax * rnpb[2][1];
            zs = -x;
            p = rnpb[0][0] * xs + rnpb[0][1] * ys + rnpb[0][2] * zs;
            q = rnpb[1][0] * xs + rnpb[1][1] * ys + rnpb[1][2] * zs;
            eo = ((p != 0) || (q != 0)) ? s - Math.Atan2(q, p) : s;

            return eo;
        }


        /// <summary>
        /// Form rotation matrix given the Fukushima-Williams angles.
        /// </summary>
        /// <param name="gamb">F-W angle gamma_bar (radians)</param>
        /// <param name="phib">F-W angle phi_bar (radians)</param>
        /// <param name="psi">F-W angle psi (radians)</param>
        /// <param name="eps">F-W angle epsilon (radians)</param>
        /// <param name="r">double[3][3]   rotation matrix</param>
        public static void iauFw2m(double gamb, double phib, double psi, double eps,
             out double[][] r)
        {
            /* Construct the matrix. */
            iauIr(out r);
            iauRz(gamb, r);
            iauRx(phib, r);
            iauRz(-psi, r);
            iauRx(-eps, r);
        }


        /// <summary>
        /// CIP X,Y given Fukushima-Williams bias-precession-nutation angles.
        /// </summary>
        /// <param name="gamb">F-W angle gamma_bar (radians)</param>
        /// <param name="phib">F-W angle phi_bar (radians)</param>
        /// <param name="psi">F-W angle psi (radians)</param>
        /// <param name="eps">F-W angle epsilon (radians)</param>
        /// <param name="x">CIP unit vector X</param>
        /// <param name="y">CIP unit vector Y</param>
        public static void iauFw2xy(double gamb, double phib, double psi, double eps,
              out double x, out double y)
        {
            double[][] r;


            /* Form NxPxB matrix. */
            iauFw2m(gamb, phib, psi, eps, out r);

            /* Extract CIP X,Y. */
            iauBpn2xy(r, out x, out y);
        }


        /// <summary>
        /// Long-term precession matrix.
        /// </summary>
        /// <param name="epj">Julian epoch (TT)</param>
        /// <param name="rp">double[3][3]   precession matrix, J2000.0 to date</param>
        public static void iauLtp(double epj, out double[][] rp)
        {
            int i;
            double w;
            double[] peqr;
            double[] pecl;
            double[] v = new double[3];
            double[] eqx = new double[3];

            rp = new double[][] { new double[3], new double[3], new double[3] };

            /* Equator pole (bottom row of matrix). */
            iauLtpequ(epj, out peqr);

            /* Ecliptic pole. */
            iauLtpecl(epj, out pecl);

            /* Equinox (top row of matrix). */
            iauPxp(peqr, pecl, ref v);
            iauPn(v, out w, ref eqx);

            /* Middle row of matrix. */
            iauPxp(peqr, eqx, ref v);

            /* Assemble the matrix. */
            for (i = 0; i < 3; i++)
            {
                rp[0][i] = eqx[i];
                rp[1][i] = v[i];
                rp[2][i] = peqr[i];
            }
        }


        /// <summary>
        /// Long-term precession matrix, including ICRS frame bias.
        /// </summary>
        /// <param name="epj">Julian epoch (TT)</param>
        /// <param name="rpb">double[3][3]   precession-bias matrix, J2000.0 to date</param>
        public static void iauLtpb(double epj, out double[][] rpb)
        {
            /* Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33) */
            double dx = -0.016617 * DAS2R,
                         de = -0.0068192 * DAS2R,
                         dr = -0.0146 * DAS2R;

            int i;
            double[][] rp;

            rpb = new double[][] { new double[3], new double[3], new double[3] };

            /* Precession matrix. */
            iauLtp(epj, out rp);

            /* Apply the bias. */
            for (i = 0; i < 3; i++)
            {
                rpb[i][0] = rp[i][0] - rp[i][1] * dr + rp[i][2] * dx;
                rpb[i][1] = rp[i][0] * dr + rp[i][1] + rp[i][2] * de;
                rpb[i][2] = -rp[i][0] * dx - rp[i][1] * de + rp[i][2];
            }
        }


        /// <summary>
        /// Long-term precession of the ecliptic.
        /// </summary>
        /// <param name="epj">Julian epoch (TT)</param>
        /// <param name="vec">double[3]      ecliptic pole unit vector</param>
        public static void iauLtpecl(double epj, out double[] vec)
        {
            vec = new double[3];

            /* Obliquity at J2000.0 (radians). */
            double eps0 = 84381.406 * DAS2R;

            /* Polynomial coefficients */
            //enum { NPOL = 4 };
            double[][] pqpol = {
                new double[]
                {
                    5851.607687,
                    -0.1189000,
                    -0.00028913,
                    0.000000101
                },
                new double[]
                {
                    -1600.886300,
                    1.1689818,
                    -0.00000020,
                    -0.000000437
                }
            };
            int NPOL = pqpol[0].Length;

            /* Periodic coefficients */
            double[][] pqper = {
                new double[]{ 708.15,-5486.751211,-684.661560,  667.666730,-5523.863691},
                new double[]{2309.00,  -17.127623,2446.283880,-2354.886252, -549.747450},
                new double[]{1620.00, -617.517403, 399.671049, -428.152441, -310.998056},
                new double[]{ 492.20,  413.442940,-356.652376,  376.202861,  421.535876},
                new double[]{1183.00,   78.614193,-186.387003,  184.778874,  -36.776172},
                new double[]{ 622.00, -180.732815,-316.800070,  335.321713, -145.278396},
                new double[]{ 882.00,  -87.676083, 198.296701, -185.138669,  -34.744450},
                new double[]{ 547.00,   46.140315, 101.135679, -120.972830,   22.885731}
            };

            //static int NPER = (int) ( sizeof pqper / 5 / sizeof (double) );
            int NPER = (int)(pqper.Length);

            /* Miscellaneous */
            int i;
            double t, p, q, w, a, s, c;


            /* Centuries since J2000. */
            t = (epj - 2000.0) / 100.0;

            /* Initialize P_A and Q_A accumulators. */
            p = 0.0;
            q = 0.0;

            /* Periodic terms. */
            w = D2PI * t;
            for (i = 0; i < NPER; i++)
            {
                a = w / pqper[i][0];
                s = Math.Sin(a);
                c = Math.Cos(a);
                p += c * pqper[i][1] + s * pqper[i][3];
                q += c * pqper[i][2] + s * pqper[i][4];
            }

            /* Polynomial terms. */
            w = 1.0;
            for (i = 0; i < NPOL; i++)
            {
                p += pqpol[0][i] * w;
                q += pqpol[1][i] * w;
                w *= t;
            }

            /* P_A and Q_A (radians). */
            p *= DAS2R;
            q *= DAS2R;

            /* Form the ecliptic pole vector. */
            w = 1.0 - p * p - q * q;
            w = w < 0.0 ? 0.0 : Math.Sqrt(w);
            s = Math.Sin(eps0);
            c = Math.Cos(eps0);
            vec[0] = p;
            vec[1] = -q * c - w * s;
            vec[2] = -q * s + w * c;
        }


        /// <summary>
        /// Long-term precession of the equator.
        /// </summary>
        /// <param name="epj">Julian epoch (TT)</param>
        /// <param name="veq">double[3]      equator pole unit vector</param>
        public static void iauLtpequ(double epj, out double[] veq)
        {
            veq = new double[3];

            /* Polynomial coefficients */
            //enum { NPOL = 4 };
            double[][] xypol =
            {
                new double[] {
                    5453.282155,
                    0.4252841,
                    -0.00037173,
                    -0.000000152
                },
                new double[] {
                    -73750.930350,
                    -0.7675452,
                    -0.00018725,
                    0.000000231
                }
            };
            int NPOL = xypol[0].Length;

            /* Periodic coefficients */
            double[][] xyper = {
                new double[] { 256.75, -819.940624,75004.344875,81491.287984, 1558.515853},
                new double[] { 708.15,-8444.676815,  624.033993,  787.163481, 7774.939698},
                new double[] { 274.20, 2600.009459, 1251.136893, 1251.296102,-2219.534038},
                new double[] { 241.45, 2755.175630,-1102.212834,-1257.950837,-2523.969396},
                new double[] {2309.00, -167.659835,-2660.664980,-2966.799730,  247.850422},
                new double[] { 492.20,  871.855056,  699.291817,  639.744522, -846.485643},
                new double[] { 396.10,   44.769698,  153.167220,  131.600209,-1393.124055},
                new double[] { 288.90, -512.313065, -950.865637, -445.040117,  368.526116},
                new double[] { 231.10, -819.415595,  499.754645,  584.522874,  749.045012},
                new double[] {1610.00, -538.071099, -145.188210,  -89.756563,  444.704518},
                new double[] { 620.00, -189.793622,  558.116553,  524.429630,  235.934465},
                new double[] { 157.87, -402.922932,  -23.923029,  -13.549067,  374.049623},
                new double[] { 220.30,  179.516345, -165.405086, -210.157124, -171.330180},
                new double[] {1200.00,   -9.814756,    9.344131,  -44.919798,  -22.899655}
            };
            //static const int NPER = (int)(sizeof xyper / 5 / sizeof(double));
            int NPER = xyper.Length;

            /* Miscellaneous */
            int i;
            double t, x, y, w, a, s, c;


            /* Centuries since J2000. */
            t = (epj - 2000.0) / 100.0;

            /* Initialize X and Y accumulators. */
            x = 0.0;
            y = 0.0;

            /* Periodic terms. */
            w = D2PI * t;
            for (i = 0; i < NPER; i++)
            {
                a = w / xyper[i][0];
                s = Math.Sin(a);
                c = Math.Cos(a);
                x += c * xyper[i][1] + s * xyper[i][3];
                y += c * xyper[i][2] + s * xyper[i][4];
            }

            /* Polynomial terms. */
            w = 1.0;
            for (i = 0; i < NPOL; i++)
            {
                x += xypol[0][i] * w;
                y += xypol[1][i] * w;
                w *= t;
            }

            /* X and Y (direction cosines). */
            x *= DAS2R;
            y *= DAS2R;

            /* Form the equator pole vector. */
            veq[0] = x;
            veq[1] = y;
            w = 1.0 - x * x - y * y;
            veq[2] = w < 0.0 ? 0.0 : Math.Sqrt(w);
        }


        /// <summary>
        /// Form the matrix of nutation for a given date, IAU 2000A model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rmatn">double[3][3]    nutation matrix</param>
        public static void iauNum00a(double date1, double date2, out double[][] rmatn)
        {
            double dpsi, deps, epsa;
            double[][] rb;
            double[][] rp;
            double[][] rbp;
            double[][] rbpn;


            /* Obtain the required matrix (discarding other results). */
            iauPn00a(date1, date2,
                     out dpsi, out deps, out epsa, out rb, out rp, out rbp, out rmatn, out rbpn);
        }


        /// <summary>
        /// Form the matrix of nutation for a given date, IAU 2000B model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rmatn">double[3][3]   nutation matrix</param>
        public static void iauNum00b(double date1, double date2, out double[][] rmatn)
        {
            double dpsi, deps, epsa;
            double[][] rb;
            double[][] rp;
            double[][] rbp;
            double[][] rbpn;


            /* Obtain the required matrix (discarding other results). */
            iauPn00b(date1, date2,
                     out dpsi, out deps, out epsa, out rb, out rp, out rbp, out rmatn, out rbpn);
        }


        /// <summary>
        /// Form the matrix of nutation for a given date, IAU 2006/2000A model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rmatn">double[3][3]    nutation matrix</param>
        public static void iauNum06a(double date1, double date2, out double[][] rmatn)
        {
            double eps, dp, de;


            /* Mean obliquity. */
            eps = iauObl06(date1, date2);

            /* Nutation components. */
            iauNut06a(date1, date2, out dp, out de);

            /* Nutation matrix. */
            iauNumat(eps, dp, de, out rmatn);
        }


        /// <summary>
        /// Form the matrix of nutation.
        /// </summary>
        /// <param name="epsa">mean obliquity of date</param>
        /// <param name="dpsi">longitude (radians)</param>
        /// <param name="deps">obliquity (radians)</param>
        /// <param name="rmatn">double[3][3]   nutation matrix</param>
        public static void iauNumat(double epsa, double dpsi, double deps, out double[][] rmatn)
        {
            /* Build the rotation matrix. */
            iauIr(out rmatn);
            iauRx(epsa, rmatn);
            iauRz(-dpsi, rmatn);
            iauRx(-(epsa + deps), rmatn);
        }


        readonly struct LuniSolar
        {
            public readonly int nl, nlp, nf, nd, nom; /* coefficients of l,l',F,D,Om */
            public readonly double sp, spt, cp;     /* longitude sin, t*sin, cos coefficients */
            public readonly double ce, cet, se;     /* obliquity cos, t*cos, sin coefficients */

            public LuniSolar(int nl, int nlp, int nf, int nd, int nom, double sp, double spt, double cp, double ce, double cet, double se)
            {
                this.nl = nl;
                this.nlp = nlp;
                this.nf = nf;
                this.nd = nd;
                this.nom = nom;
                this.sp = sp;
                this.spt = spt;
                this.cp = cp;
                this.ce = ce;
                this.cet = cet;
                this.se = se;
            }
        }

        readonly struct PlanetaryNutation
        {
            public readonly int nl,               /* coefficients of l, F, D and Omega */
                nf,
                nd,
                nom,
                nme,              /* coefficients of planetary longitudes */
                nve,
                nea,
                nma,
                nju,
                nsa,
                nur,
                nne,
                npa;              /* coefficient of general precession */
            public readonly int sp, cp;            /* longitude sin, cos coefficients */
            public readonly int se, ce;            /* obliquity sin, cos coefficients */

            public PlanetaryNutation(int nl, int nf, int nd, int nom, int nme, int nve, int nea, int nma, int nju, int nsa, int nur, int nne, int npa, int sp, int cp, int se, int ce)
            {
                this.nl = nl;
                this.nf = nf;
                this.nd = nd;
                this.nom = nom;
                this.nme = nme;
                this.nve = nve;
                this.nea = nea;
                this.nma = nma;
                this.nju = nju;
                this.nsa = nsa;
                this.nur = nur;
                this.nne = nne;
                this.npa = npa;
                this.sp = sp;
                this.cp = cp;
                this.se = se;
                this.ce = ce;
            }
        }

        /* ------------------------- */
        /* Luni-Solar nutation model */
        /* ------------------------- */

        /* The units for the sine and cosine coefficients are */
        /* 0.1 microarcsecond and the same per Julian century */

        static readonly LuniSolar[] Nut00axls =
        {
            /* 1- 10 */
            new LuniSolar( 0, 0, 0, 0, 1,
            -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0),
            new LuniSolar( 0, 0, 2,-2, 2,
            -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0),
            new LuniSolar( 0, 0, 2, 0, 2,-2276413.0,-234.0,2796.0,978459.0,-485.0, 1374.0),
            new LuniSolar( 0, 0, 0, 0, 2,2074554.0, 207.0, -698.0,-897492.0,470.0, -291.0),
            new LuniSolar( 0, 1, 0, 0, 0,1475877.0,-3633.0,11817.0,73871.0,-184.0,-1924.0),
            new LuniSolar( 0, 1, 2,-2, 2,-516821.0,1226.0, -524.0,224386.0,-677.0, -174.0),
            new LuniSolar( 1, 0, 0, 0, 0, 711159.0,  73.0, -872.0,  -6750.0,  0.0,  358.0),
            new LuniSolar( 0, 0, 2, 0, 1,-387298.0,-367.0,  380.0, 200728.0, 18.0,  318.0),
            new LuniSolar( 1, 0, 2, 0, 2,-301461.0, -36.0,  816.0, 129025.0,-63.0,  367.0),
            new LuniSolar( 0,-1, 2,-2, 2, 215829.0,-494.0,  111.0, -95929.0,299.0,  132.0),

            /* 11-20 */
            new LuniSolar( 0, 0, 2,-2, 1, 128227.0, 137.0,  181.0, -68982.0, -9.0,   39.0),
            new LuniSolar(-1, 0, 2, 0, 2, 123457.0,  11.0,   19.0, -53311.0, 32.0,   -4.0),
            new LuniSolar(-1, 0, 0, 2, 0, 156994.0,  10.0, -168.0,  -1235.0,  0.0,   82.0),
            new LuniSolar( 1, 0, 0, 0, 1,  63110.0,  63.0,   27.0, -33228.0,  0.0,   -9.0),
            new LuniSolar(-1, 0, 0, 0, 1, -57976.0, -63.0, -189.0,  31429.0,  0.0,  -75.0),
            new LuniSolar(-1, 0, 2, 2, 2, -59641.0, -11.0,  149.0,  25543.0,-11.0,   66.0),
            new LuniSolar( 1, 0, 2, 0, 1, -51613.0, -42.0,  129.0,  26366.0,  0.0,   78.0),
            new LuniSolar(-2, 0, 2, 0, 1,  45893.0,  50.0,   31.0, -24236.0,-10.0,   20.0),
            new LuniSolar( 0, 0, 0, 2, 0,  63384.0,  11.0, -150.0,  -1220.0,  0.0,   29.0),
            new LuniSolar( 0, 0, 2, 2, 2, -38571.0,  -1.0,  158.0,  16452.0,-11.0,   68.0),

            /* 21-30 */
            new LuniSolar( 0,-2, 2,-2, 2,  32481.0,   0.0,    0.0, -13870.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 2, 0, -47722.0,   0.0,  -18.0,    477.0,  0.0,  -25.0),
            new LuniSolar( 2, 0, 2, 0, 2, -31046.0,  -1.0,  131.0,  13238.0,-11.0,   59.0),
            new LuniSolar( 1, 0, 2,-2, 2,  28593.0,   0.0,   -1.0, -12338.0, 10.0,   -3.0),
            new LuniSolar(-1, 0, 2, 0, 1,  20441.0,  21.0,   10.0, -10758.0,  0.0,   -3.0),
            new LuniSolar( 2, 0, 0, 0, 0,  29243.0,   0.0,  -74.0,   -609.0,  0.0,   13.0),
            new LuniSolar( 0, 0, 2, 0, 0,  25887.0,   0.0,  -66.0,   -550.0,  0.0,   11.0),
            new LuniSolar( 0, 1, 0, 0, 1, -14053.0, -25.0,   79.0,   8551.0, -2.0,  -45.0),
            new LuniSolar(-1, 0, 0, 2, 1,  15164.0,  10.0,   11.0,  -8001.0,  0.0,   -1.0),
            new LuniSolar( 0, 2, 2,-2, 2, -15794.0,  72.0,  -16.0,   6850.0,-42.0,   -5.0),

            /* 31-40 */
            new LuniSolar( 0, 0,-2, 2, 0,  21783.0,   0.0,   13.0,   -167.0,  0.0,   13.0),
            new LuniSolar( 1, 0, 0,-2, 1, -12873.0, -10.0,  -37.0,   6953.0,  0.0,  -14.0),
            new LuniSolar( 0,-1, 0, 0, 1, -12654.0,  11.0,   63.0,   6415.0,  0.0,   26.0),
            new LuniSolar(-1, 0, 2, 2, 1, -10204.0,   0.0,   25.0,   5222.0,  0.0,   15.0),
            new LuniSolar( 0, 2, 0, 0, 0,  16707.0, -85.0,  -10.0,    168.0, -1.0,   10.0),
            new LuniSolar( 1, 0, 2, 2, 2,  -7691.0,   0.0,   44.0,   3268.0,  0.0,   19.0),
            new LuniSolar(-2, 0, 2, 0, 0, -11024.0,   0.0,  -14.0,    104.0,  0.0,    2.0),
            new LuniSolar( 0, 1, 2, 0, 2,   7566.0, -21.0,  -11.0,  -3250.0,  0.0,   -5.0),
            new LuniSolar( 0, 0, 2, 2, 1,  -6637.0, -11.0,   25.0,   3353.0,  0.0,   14.0),
            new LuniSolar( 0,-1, 2, 0, 2,  -7141.0,  21.0,    8.0,   3070.0,  0.0,    4.0),

            /* 41-50 */
            new LuniSolar( 0, 0, 0, 2, 1,  -6302.0, -11.0,    2.0,   3272.0,  0.0,    4.0),
            new LuniSolar( 1, 0, 2,-2, 1,   5800.0,  10.0,    2.0,  -3045.0,  0.0,   -1.0),
            new LuniSolar( 2, 0, 2,-2, 2,   6443.0,   0.0,   -7.0,  -2768.0,  0.0,   -4.0),
            new LuniSolar(-2, 0, 0, 2, 1,  -5774.0, -11.0,  -15.0,   3041.0,  0.0,   -5.0),
            new LuniSolar( 2, 0, 2, 0, 1,  -5350.0,   0.0,   21.0,   2695.0,  0.0,   12.0),
            new LuniSolar( 0,-1, 2,-2, 1,  -4752.0, -11.0,   -3.0,   2719.0,  0.0,   -3.0),
            new LuniSolar( 0, 0, 0,-2, 1,  -4940.0, -11.0,  -21.0,   2720.0,  0.0,   -9.0),
            new LuniSolar(-1,-1, 0, 2, 0,   7350.0,   0.0,   -8.0,    -51.0,  0.0,    4.0),
            new LuniSolar( 2, 0, 0,-2, 1,   4065.0,   0.0,    6.0,  -2206.0,  0.0,    1.0),
            new LuniSolar( 1, 0, 0, 2, 0,   6579.0,   0.0,  -24.0,   -199.0,  0.0,    2.0),

            /* 51-60 */
            new LuniSolar( 0, 1, 2,-2, 1,   3579.0,   0.0,    5.0,  -1900.0,  0.0,    1.0),
            new LuniSolar( 1,-1, 0, 0, 0,   4725.0,   0.0,   -6.0,    -41.0,  0.0,    3.0),
            new LuniSolar(-2, 0, 2, 0, 2,  -3075.0,   0.0,   -2.0,   1313.0,  0.0,   -1.0),
            new LuniSolar( 3, 0, 2, 0, 2,  -2904.0,   0.0,   15.0,   1233.0,  0.0,    7.0),
            new LuniSolar( 0,-1, 0, 2, 0,   4348.0,   0.0,  -10.0,    -81.0,  0.0,    2.0),
            new LuniSolar( 1,-1, 2, 0, 2,  -2878.0,   0.0,    8.0,   1232.0,  0.0,    4.0),
            new LuniSolar( 0, 0, 0, 1, 0,  -4230.0,   0.0,    5.0,    -20.0,  0.0,   -2.0),
            new LuniSolar(-1,-1, 2, 2, 2,  -2819.0,   0.0,    7.0,   1207.0,  0.0,    3.0),
            new LuniSolar(-1, 0, 2, 0, 0,  -4056.0,   0.0,    5.0,     40.0,  0.0,   -2.0),
            new LuniSolar( 0,-1, 2, 2, 2,  -2647.0,   0.0,   11.0,   1129.0,  0.0,    5.0),

            /* 61-70 */
            new LuniSolar(-2, 0, 0, 0, 1,  -2294.0,   0.0,  -10.0,   1266.0,  0.0,   -4.0),
            new LuniSolar( 1, 1, 2, 0, 2,   2481.0,   0.0,   -7.0,  -1062.0,  0.0,   -3.0),
            new LuniSolar( 2, 0, 0, 0, 1,   2179.0,   0.0,   -2.0,  -1129.0,  0.0,   -2.0),
            new LuniSolar(-1, 1, 0, 1, 0,   3276.0,   0.0,    1.0,     -9.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 0, 0, 0,  -3389.0,   0.0,    5.0,     35.0,  0.0,   -2.0),
            new LuniSolar( 1, 0, 2, 0, 0,   3339.0,   0.0,  -13.0,   -107.0,  0.0,    1.0),
            new LuniSolar(-1, 0, 2,-2, 1,  -1987.0,   0.0,   -6.0,   1073.0,  0.0,   -2.0),
            new LuniSolar( 1, 0, 0, 0, 2,  -1981.0,   0.0,    0.0,    854.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 1, 0,   4026.0,   0.0, -353.0,   -553.0,  0.0, -139.0),
            new LuniSolar( 0, 0, 2, 1, 2,   1660.0,   0.0,   -5.0,   -710.0,  0.0,   -2.0),

            /* 71-80 */
            new LuniSolar(-1, 0, 2, 4, 2,  -1521.0,   0.0,    9.0,    647.0,  0.0,    4.0),
            new LuniSolar(-1, 1, 0, 1, 1,   1314.0,   0.0,    0.0,   -700.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 2,-2, 1,  -1283.0,   0.0,    0.0,    672.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 2, 1,  -1331.0,   0.0,    8.0,    663.0,  0.0,    4.0),
            new LuniSolar(-2, 0, 2, 2, 2,   1383.0,   0.0,   -2.0,   -594.0,  0.0,   -2.0),
            new LuniSolar(-1, 0, 0, 0, 2,   1405.0,   0.0,    4.0,   -610.0,  0.0,    2.0),
            new LuniSolar( 1, 1, 2,-2, 2,   1290.0,   0.0,    0.0,   -556.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2, 4, 2,  -1214.0,   0.0,    5.0,    518.0,  0.0,    2.0),
            new LuniSolar(-1, 0, 4, 0, 2,   1146.0,   0.0,   -3.0,   -490.0,  0.0,   -1.0),
            new LuniSolar( 2, 0, 2,-2, 1,   1019.0,   0.0,   -1.0,   -527.0,  0.0,   -1.0),

            /* 81-90 */
            new LuniSolar( 2, 0, 2, 2, 2,  -1100.0,   0.0,    9.0,    465.0,  0.0,    4.0),
            new LuniSolar( 1, 0, 0, 2, 1,   -970.0,   0.0,    2.0,    496.0,  0.0,    1.0),
            new LuniSolar( 3, 0, 0, 0, 0,   1575.0,   0.0,   -6.0,    -50.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 2,-2, 2,    934.0,   0.0,   -3.0,   -399.0,  0.0,   -1.0),
            new LuniSolar( 0, 0, 4,-2, 2,    922.0,   0.0,   -1.0,   -395.0,  0.0,   -1.0),
            new LuniSolar( 0, 1, 2, 0, 1,    815.0,   0.0,   -1.0,   -422.0,  0.0,   -1.0),
            new LuniSolar( 0, 0,-2, 2, 1,    834.0,   0.0,    2.0,   -440.0,  0.0,    1.0),
            new LuniSolar( 0, 0, 2,-2, 3,   1248.0,   0.0,    0.0,   -170.0,  0.0,    1.0),
            new LuniSolar(-1, 0, 0, 4, 0,   1338.0,   0.0,   -5.0,    -39.0,  0.0,    0.0),
            new LuniSolar( 2, 0,-2, 0, 1,    716.0,   0.0,   -2.0,   -389.0,  0.0,   -1.0),

            /* 91-100 */
            new LuniSolar(-2, 0, 0, 4, 0,   1282.0,   0.0,   -3.0,    -23.0,  0.0,    1.0),
            new LuniSolar(-1,-1, 0, 2, 1,    742.0,   0.0,    1.0,   -391.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 1, 1,   1020.0,   0.0,  -25.0,   -495.0,  0.0,  -10.0),
            new LuniSolar( 0, 1, 0, 0, 2,    715.0,   0.0,   -4.0,   -326.0,  0.0,    2.0),
            new LuniSolar( 0, 0,-2, 0, 1,   -666.0,   0.0,   -3.0,    369.0,  0.0,   -1.0),
            new LuniSolar( 0,-1, 2, 0, 1,   -667.0,   0.0,    1.0,    346.0,  0.0,    1.0),
            new LuniSolar( 0, 0, 2,-1, 2,   -704.0,   0.0,    0.0,    304.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 4, 2,   -694.0,   0.0,    5.0,    294.0,  0.0,    2.0),
            new LuniSolar(-2,-1, 0, 2, 0,  -1014.0,   0.0,   -1.0,      4.0,  0.0,   -1.0),
            new LuniSolar( 1, 1, 0,-2, 1,   -585.0,   0.0,   -2.0,    316.0,  0.0,   -1.0),

            /* 101-110 */
            new LuniSolar(-1, 1, 0, 2, 0,   -949.0,   0.0,    1.0,      8.0,  0.0,   -1.0),
            new LuniSolar(-1, 1, 0, 1, 2,   -595.0,   0.0,    0.0,    258.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0, 0, 1,    528.0,   0.0,    0.0,   -279.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 2, 2,   -590.0,   0.0,    4.0,    252.0,  0.0,    2.0),
            new LuniSolar(-1, 1, 2, 2, 2,    570.0,   0.0,   -2.0,   -244.0,  0.0,   -1.0),
            new LuniSolar( 3, 0, 2, 0, 1,   -502.0,   0.0,    3.0,    250.0,  0.0,    2.0),
            new LuniSolar( 0, 1,-2, 2, 0,   -875.0,   0.0,    1.0,     29.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0,-2, 1,   -492.0,   0.0,   -3.0,    275.0,  0.0,   -1.0),
            new LuniSolar( 0, 1, 2, 2, 2,    535.0,   0.0,   -2.0,   -228.0,  0.0,   -1.0),
            new LuniSolar(-1,-1, 2, 2, 1,   -467.0,   0.0,    1.0,    240.0,  0.0,    1.0),

            /* 111-120 */
            new LuniSolar( 0,-1, 0, 0, 2,    591.0,   0.0,    0.0,   -253.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-4, 1,   -453.0,   0.0,   -1.0,    244.0,  0.0,   -1.0),
            new LuniSolar(-1, 0,-2, 2, 0,    766.0,   0.0,    1.0,      9.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2, 2, 1,   -446.0,   0.0,    2.0,    225.0,  0.0,    1.0),
            new LuniSolar( 2,-1, 2, 0, 2,   -488.0,   0.0,    2.0,    207.0,  0.0,    1.0),
            new LuniSolar( 0, 0, 0, 2, 2,   -468.0,   0.0,    0.0,    201.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 0, 1,   -421.0,   0.0,    1.0,    216.0,  0.0,    1.0),
            new LuniSolar(-1, 1, 2, 0, 2,    463.0,   0.0,    0.0,   -200.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 2, 0,   -673.0,   0.0,    2.0,     14.0,  0.0,    0.0),
            new LuniSolar( 0,-1,-2, 2, 0,    658.0,   0.0,    0.0,     -2.0,  0.0,    0.0),

            /* 121-130 */
            new LuniSolar( 0, 3, 2,-2, 2,   -438.0,   0.0,    0.0,    188.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 1, 1,   -390.0,   0.0,    0.0,    205.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 2, 0,    639.0, -11.0,   -2.0,    -19.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2, 0, 2,    412.0,   0.0,   -2.0,   -176.0,  0.0,   -1.0),
            new LuniSolar( 1, 1, 0, 0, 1,   -361.0,   0.0,    0.0,    189.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2, 0, 1,    360.0,   0.0,   -1.0,   -185.0,  0.0,   -1.0),
            new LuniSolar( 2, 0, 0, 2, 0,    588.0,   0.0,   -3.0,    -24.0,  0.0,    0.0),
            new LuniSolar( 1, 0,-2, 2, 0,   -578.0,   0.0,    1.0,      5.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 2, 2,   -396.0,   0.0,    0.0,    171.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 1, 0,    565.0,   0.0,   -1.0,     -6.0,  0.0,    0.0),

            /* 131-140 */
            new LuniSolar( 0, 1, 0,-2, 1,   -335.0,   0.0,   -1.0,    184.0,  0.0,   -1.0),
            new LuniSolar(-1, 0, 2,-2, 2,    357.0,   0.0,    1.0,   -154.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0,-1, 1,    321.0,   0.0,    1.0,   -174.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0, 0, 1,   -301.0,   0.0,   -1.0,    162.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-1, 2,   -334.0,   0.0,    0.0,    144.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0, 2, 0,    493.0,   0.0,   -2.0,    -15.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 4, 0,    494.0,   0.0,   -2.0,    -19.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 1, 2,    337.0,   0.0,   -1.0,   -143.0,  0.0,   -1.0),
            new LuniSolar( 0, 0, 2, 1, 1,    280.0,   0.0,   -1.0,   -144.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0,-2, 2,    309.0,   0.0,    1.0,   -134.0,  0.0,    0.0),

            /* 141-150 */
            new LuniSolar(-1, 0, 2, 4, 1,   -263.0,   0.0,    2.0,    131.0,  0.0,    1.0),
            new LuniSolar( 1, 0,-2, 0, 1,    253.0,   0.0,    1.0,   -138.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2,-2, 1,    245.0,   0.0,    0.0,   -128.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 2, 0,    416.0,   0.0,   -2.0,    -17.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2,-1, 1,   -229.0,   0.0,    0.0,    128.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2, 2, 1,    231.0,   0.0,    0.0,   -120.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 2, 0, 2,   -259.0,   0.0,    2.0,    109.0,  0.0,    1.0),
            new LuniSolar( 2,-1, 0, 0, 0,    375.0,   0.0,   -1.0,     -8.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2,-2, 2,    252.0,   0.0,    0.0,   -108.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2, 1, 2,   -245.0,   0.0,    1.0,    104.0,  0.0,    0.0),

            /* 151-160 */
            new LuniSolar( 1, 0, 4,-2, 2,    243.0,   0.0,   -1.0,   -104.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 0, 1,    208.0,   0.0,    1.0,   -112.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 2, 1,    199.0,   0.0,    0.0,   -102.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2, 4, 1,   -208.0,   0.0,    1.0,    105.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 0, 0,    335.0,   0.0,   -2.0,    -14.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0, 1, 0,   -325.0,   0.0,    1.0,      7.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 4, 1,   -187.0,   0.0,    0.0,     96.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 4, 0, 1,    197.0,   0.0,   -1.0,   -100.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 2, 1,   -192.0,   0.0,    2.0,     94.0,  0.0,    1.0),
            new LuniSolar( 0, 0, 2,-3, 2,   -188.0,   0.0,    0.0,     83.0,  0.0,    0.0),

            /* 161-170 */
            new LuniSolar(-1,-2, 0, 2, 0,    276.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 0, 0, 0,   -286.0,   0.0,    1.0,      6.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4, 0, 2,    186.0,   0.0,   -1.0,    -79.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 0, 3,   -219.0,   0.0,    0.0,     43.0,  0.0,    0.0),
            new LuniSolar( 0, 3, 0, 0, 0,    276.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2,-4, 1,   -153.0,   0.0,   -1.0,     84.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 0, 2, 1,   -156.0,   0.0,    0.0,     81.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 4, 1,   -154.0,   0.0,    1.0,     78.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 4, 2,   -174.0,   0.0,    1.0,     75.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 4, 2,   -163.0,   0.0,    2.0,     69.0,  0.0,    1.0),

            /* 171-180 */
            new LuniSolar(-2, 2, 0, 2, 0,   -228.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 0, 1,     91.0,   0.0,   -4.0,    -54.0,  0.0,   -2.0),
            new LuniSolar(-2, 0, 0, 2, 2,    175.0,   0.0,    0.0,    -75.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 0, 2,   -159.0,   0.0,    0.0,     69.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4,-2, 1,    141.0,   0.0,    0.0,    -72.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 2,-2, 1,    147.0,   0.0,    0.0,    -75.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 0, 2, 1,   -132.0,   0.0,    0.0,     69.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0,-1, 1,    159.0,   0.0,  -28.0,    -54.0,  0.0,   11.0),
            new LuniSolar( 0,-2, 0, 2, 0,    213.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 4, 1,    123.0,   0.0,    0.0,    -64.0,  0.0,    0.0),

            /* 181-190 */
            new LuniSolar(-3, 0, 0, 0, 1,   -118.0,   0.0,   -1.0,     66.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2, 2, 2,    144.0,   0.0,   -1.0,    -61.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 4, 1,   -121.0,   0.0,    1.0,     60.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 2, 2, 2,   -134.0,   0.0,    1.0,     56.0,  0.0,    1.0),
            new LuniSolar(-1, 1, 2,-2, 1,   -105.0,   0.0,    0.0,     57.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0,-4, 1,   -102.0,   0.0,    0.0,     56.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0,-2, 2,    120.0,   0.0,    0.0,    -52.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2,-4, 1,    101.0,   0.0,    0.0,    -54.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0, 2, 1,   -113.0,   0.0,    0.0,     59.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2,-1, 1,   -106.0,   0.0,    0.0,     61.0,  0.0,    0.0),

            /* 191-200 */
            new LuniSolar( 0,-2, 2, 2, 2,   -129.0,   0.0,    1.0,     55.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0, 2, 1,   -114.0,   0.0,    0.0,     57.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 2,-2, 2,    113.0,   0.0,   -1.0,    -49.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0,-2, 2,   -102.0,   0.0,    0.0,     44.0,  0.0,    0.0),
            new LuniSolar( 0, 2, 0, 0, 1,    -94.0,   0.0,    0.0,     51.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0,-4, 1,   -100.0,   0.0,   -1.0,     56.0,  0.0,    0.0),
            new LuniSolar( 0, 2, 2,-2, 1,     87.0,   0.0,    0.0,    -47.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0, 4, 0,    161.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2, 0, 1,     96.0,   0.0,    0.0,    -50.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 4, 0,    151.0,   0.0,   -1.0,     -5.0,  0.0,    0.0),

            /* 201-210 */
            new LuniSolar(-1,-2, 2, 2, 2,   -104.0,   0.0,    0.0,     44.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 4, 2,   -110.0,   0.0,    0.0,     48.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 2, 1,   -100.0,   0.0,    1.0,     50.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 0, 2, 0,     92.0,   0.0,   -5.0,     12.0,  0.0,   -2.0),
            new LuniSolar(-2, 1, 2, 0, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 0,-2, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 0, 1,    -78.0,   0.0,    0.0,     41.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2,-2, 1,    -77.0,   0.0,    0.0,     43.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0, 2, 2,      2.0,   0.0,    0.0,     54.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2,-1, 2,     94.0,   0.0,    0.0,    -40.0,  0.0,    0.0),

            /* 211-220 */
            new LuniSolar(-1, 0, 4,-2, 2,    -93.0,   0.0,    0.0,     40.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 2, 0, 2,    -83.0,   0.0,   10.0,     40.0,  0.0,   -2.0),
            new LuniSolar(-1, 0, 2, 1, 2,     83.0,   0.0,    0.0,    -36.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0, 0, 2,    -91.0,   0.0,    0.0,     39.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 0, 3,    128.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 4, 0, 2,    -79.0,   0.0,    0.0,     34.0,  0.0,    0.0),
            new LuniSolar(-1, 0,-2, 0, 1,    -83.0,   0.0,    0.0,     47.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -44.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 0, 0, 1,     83.0,   0.0,    0.0,    -43.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 3, 2,     91.0,   0.0,    0.0,    -39.0,  0.0,    0.0),

            /* 221-230 */
            new LuniSolar( 2,-1, 2, 0, 1,    -77.0,   0.0,    0.0,     39.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -43.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2, 4, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 2, 2, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0),
            new LuniSolar( 0, 2,-2, 2, 0,    -94.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2,-1, 1,     68.0,   0.0,    0.0,    -36.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 0, 0, 1,    -61.0,   0.0,    0.0,     32.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-4, 2,     71.0,   0.0,    0.0,    -31.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0,-2, 1,     62.0,   0.0,    0.0,    -34.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 0, 1,    -63.0,   0.0,    0.0,     33.0,  0.0,    0.0),

            /* 231-240 */
            new LuniSolar( 1,-1, 2,-2, 2,    -73.0,   0.0,    0.0,     32.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 0, 4, 0,    115.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 3, 0,   -103.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 2, 2,     63.0,   0.0,    0.0,    -28.0,  0.0,    0.0),
            new LuniSolar( 0, 2, 2, 0, 2,     74.0,   0.0,    0.0,    -32.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 0, 2, 0,   -103.0,   0.0,   -3.0,      3.0,  0.0,   -1.0),
            new LuniSolar( 2, 0, 2,-1, 2,    -69.0,   0.0,    0.0,     30.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 1, 1,     57.0,   0.0,    0.0,    -29.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 0, 0, 0,     94.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2, 0, 1,     64.0,   0.0,    0.0,    -33.0,  0.0,    0.0),

            /* 241-250 */
            new LuniSolar( 3,-1, 2, 0, 2,    -63.0,   0.0,    0.0,     26.0,  0.0,    0.0),
            new LuniSolar(-2, 2, 0, 2, 1,    -38.0,   0.0,    0.0,     20.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-3, 1,    -43.0,   0.0,    0.0,     24.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2,-4, 1,    -45.0,   0.0,    0.0,     23.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2,-2, 1,     47.0,   0.0,    0.0,    -24.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 0,-1, 1,    -48.0,   0.0,    0.0,     25.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 0,-2, 1,     45.0,   0.0,    0.0,    -26.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 0, 2,     56.0,   0.0,    0.0,    -25.0,  0.0,    0.0),
            new LuniSolar(-2, 0,-2, 2, 0,     88.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1, 0,-2, 4, 0,    -75.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 251-260 */
            new LuniSolar( 1,-2, 0, 0, 0,     85.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 1, 1,     49.0,   0.0,    0.0,    -26.0,  0.0,    0.0),
            new LuniSolar(-1, 2, 0, 2, 0,    -74.0,   0.0,   -3.0,     -1.0,  0.0,   -1.0),
            new LuniSolar( 1,-1, 2,-2, 1,    -39.0,   0.0,    0.0,     21.0,  0.0,    0.0),
            new LuniSolar( 1, 2, 2,-2, 2,     45.0,   0.0,    0.0,    -20.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 2,-2, 2,     51.0,   0.0,    0.0,    -22.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-1, 1,    -40.0,   0.0,    0.0,     21.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2,-2, 1,     41.0,   0.0,    0.0,    -21.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0,-2, 1,    -42.0,   0.0,    0.0,     24.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 2, 0, 2,    -51.0,   0.0,    0.0,     22.0,  0.0,    0.0),

            /* 261-270 */
            new LuniSolar( 0, 1, 2, 1, 1,    -42.0,   0.0,    0.0,     22.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 4,-2, 1,     39.0,   0.0,    0.0,    -21.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 4, 2, 2,     46.0,   0.0,    0.0,    -18.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2, 1, 2,    -53.0,   0.0,    0.0,     22.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0, 4, 0,     82.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 2, 0,     81.0,   0.0,   -1.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 1, 2,     47.0,   0.0,    0.0,    -19.0,  0.0,    0.0),
            new LuniSolar( 3, 1, 2, 0, 2,     53.0,   0.0,    0.0,    -23.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 2, 0, 1,    -45.0,   0.0,    0.0,     22.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 0, 0,    -44.0,   0.0,    0.0,     -2.0,  0.0,    0.0),

            /* 271-280 */
            new LuniSolar( 0, 1,-2, 2, 1,    -33.0,   0.0,    0.0,     16.0,  0.0,    0.0),
            new LuniSolar( 1, 0,-2, 1, 0,    -61.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 0,-1,-2, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 0,-2, 1,    -38.0,   0.0,    0.0,     19.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2,-1, 2,    -33.0,   0.0,    0.0,     21.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-3, 2,    -60.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2,-2, 3,     48.0,   0.0,    0.0,    -10.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2,-3, 1,     27.0,   0.0,    0.0,    -14.0,  0.0,    0.0),
            new LuniSolar(-1, 0,-2, 2, 1,     38.0,   0.0,    0.0,    -20.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2,-4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0),

            /* 281-290 */
            new LuniSolar(-2, 1, 0, 0, 1,    -29.0,   0.0,    0.0,     15.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0,-1, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2,-4, 2,    -32.0,   0.0,    0.0,     15.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4,-4, 4,     45.0,   0.0,    0.0,     -8.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4,-4, 2,    -44.0,   0.0,    0.0,     19.0,  0.0,    0.0),
            new LuniSolar(-1,-2, 0, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 3, 0,    -51.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0,-2, 2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 2, 2,     44.0,   0.0,    0.0,    -19.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 2, 1,     26.0,   0.0,    0.0,    -14.0,  0.0,    0.0),

            /* 291-300 */
            new LuniSolar(-2, 0, 2, 2, 0,    -60.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 0, 0, 1,     35.0,   0.0,    0.0,    -18.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 2, 2, 2,    -27.0,   0.0,    0.0,     11.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 0, 1, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 4,-2, 2,     36.0,   0.0,    0.0,    -15.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0,-2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0,-4, 1,    -35.0,   0.0,    0.0,     19.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0, 2, 1,    -37.0,   0.0,    0.0,     19.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 0, 2, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0),
            new LuniSolar(-1, 2, 2, 2, 2,     35.0,   0.0,    0.0,    -14.0,  0.0,    0.0),

            /* 301-310 */
            new LuniSolar( 3, 1, 2,-2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 0, 4, 0,     65.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 0, 2, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4, 0, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 4,-2, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 4, 1,    -30.0,   0.0,    0.0,     15.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0, 4, 1,    -32.0,   0.0,    0.0,     16.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 2, 2, 2,    -31.0,   0.0,    0.0,     13.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 3, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2, 4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0),

            /* 311-320 */
            new LuniSolar( 3, 0, 0, 2, 0,     49.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 4, 2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2, 2, 1,     23.0,   0.0,    0.0,    -12.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2, 6, 2,    -43.0,   0.0,    0.0,     18.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2, 2, 2,     26.0,   0.0,    0.0,    -11.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 6, 2,    -32.0,   0.0,    0.0,     14.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 4, 1,    -29.0,   0.0,    0.0,     14.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 4, 2,    -27.0,   0.0,    0.0,     12.0,  0.0,    0.0),
            new LuniSolar( 1, 1,-2, 1, 0,     30.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-3, 1, 2, 1, 2,    -11.0,   0.0,    0.0,      5.0,  0.0,    0.0),

            /* 321-330 */
            new LuniSolar( 2, 0,-2, 0, 2,    -21.0,   0.0,    0.0,     10.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 1, 2,    -34.0,   0.0,    0.0,     15.0,  0.0,    0.0),
            new LuniSolar(-4, 0, 2, 2, 1,    -10.0,   0.0,    0.0,      6.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 1, 0,    -36.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 0,-2, 2, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0,-1, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2,-2, 3,    -21.0,   0.0,    0.0,      5.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 2, 0, 0,    -29.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2,-2, 4,    -15.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-2,-2, 0, 2, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 331-340 */
            new LuniSolar(-2, 0,-2, 4, 0,     28.0,   0.0,    0.0,      0.0,  0.0,   -2.0),
            new LuniSolar( 0,-2,-2, 2, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 2, 0,-2, 1,    -22.0,   0.0,    0.0,     12.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 0,-4, 1,    -14.0,   0.0,    0.0,      7.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2,-2, 2,     24.0,   0.0,    0.0,    -11.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2,-4, 1,     11.0,   0.0,    0.0,     -6.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 0,-2, 2,     14.0,   0.0,    0.0,     -6.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 0, 0,     24.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 0, 2,     18.0,   0.0,    0.0,     -8.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 1, 0,    -38.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 341-350 */
            new LuniSolar( 0, 0,-2, 1, 0,    -31.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0, 2, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0),
            new LuniSolar(-1,-1,-2, 2, 0,     29.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2,-4, 1,    -18.0,   0.0,    0.0,     10.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 0,-4, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0),
            new LuniSolar( 0, 2, 0,-2, 1,    -17.0,   0.0,    0.0,     10.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0,-3, 1,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2,-2, 2,     16.0,   0.0,    0.0,     -6.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 0, 0, 1,     22.0,   0.0,    0.0,    -12.0,  0.0,    0.0),
            new LuniSolar(-4, 0, 0, 2, 0,     20.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 351-360 */
            new LuniSolar( 1, 1, 0,-4, 1,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2,-4, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4,-4, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0),
            new LuniSolar( 0, 3, 2,-2, 2,      0.0,   0.0,    0.0,     -7.0,  0.0,    0.0),
            new LuniSolar(-3,-1, 0, 4, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0, 4, 1,     19.0,   0.0,    0.0,    -10.0,  0.0,    0.0),
            new LuniSolar( 1,-1,-2, 2, 0,    -34.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 2, 2,    -20.0,   0.0,    0.0,      8.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 0, 0, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0, 0, 2,    -18.0,   0.0,    0.0,      7.0,  0.0,    0.0),

            /* 361-370 */
            new LuniSolar( 0, 0, 0, 1, 2,     13.0,   0.0,    0.0,     -6.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 0, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 2,-2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2,-1, 1,     15.0,   0.0,    0.0,     -8.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 0, 3,    -11.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 0, 0, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2, 0, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 2, 0, 0, 0,    -35.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 2, 2, 0, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 4,-2, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0),

            /* 371-380 */
            new LuniSolar( 3, 0, 2,-4, 2,    -26.0,   0.0,    0.0,     11.0,  0.0,    0.0),
            new LuniSolar( 1, 2, 2,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 4,-4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 0, 4, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 0, 2, 2,    -21.0,   0.0,    0.0,      9.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 0, 4, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 2, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0),
            new LuniSolar( 2, 0,-2, 2, 0,    -29.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0, 1, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0),

            /* 381-390 */
            new LuniSolar( 1,-1, 2,-1, 2,     22.0,   0.0,    0.0,     -9.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 4, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 0, 0, 1,    -20.0,   0.0,    0.0,     11.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2, 0, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 4,-2, 2,    -17.0,   0.0,    0.0,      7.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4,-2, 4,     15.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 0, 2, 2, 0, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0, 6, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 4, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 0, 2, 0,     25.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 391-400 */
            new LuniSolar(-1, 0, 0, 4, 2,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0),
            new LuniSolar(-1,-2, 2, 2, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0,-2, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0),
            new LuniSolar( 1, 0,-2,-2, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0),
            new LuniSolar( 0, 0,-2,-2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0),
            new LuniSolar(-2, 0,-2, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 3, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 3, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0, 4, 0,    -22.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 2, 0,     28.0,   0.0,    0.0,     -1.0,  0.0,    0.0),

            /* 401-410 */
            new LuniSolar(-2, 0, 2, 3, 2,     15.0,   0.0,    0.0,     -7.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0, 2, 2,     23.0,   0.0,    0.0,    -10.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2, 1, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0),
            new LuniSolar( 3,-1, 0, 0, 0,     29.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0, 1, 0,    -25.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 0, 0,     22.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 1, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 0, 3,     15.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 3, 1, 0, 0, 0,    -23.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 3,-1, 2,-2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0),

            /* 411-420 */
            new LuniSolar( 2, 0, 2,-1, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2, 0, 0,    -19.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4,-1, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar( 1, 2, 2, 0, 2,     21.0,   0.0,    0.0,     -9.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 6, 0,     23.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 0, 4, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 4, 1,    -19.0,   0.0,    0.0,      9.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 2, 2, 1,    -22.0,   0.0,    0.0,     10.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2, 2, 0,     27.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 3, 1,     16.0,   0.0,    0.0,     -8.0,  0.0,    0.0),

            /* 421-430 */
            new LuniSolar(-2, 1, 2, 4, 2,     19.0,   0.0,    0.0,     -8.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0, 2, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 2,-2, 2, 0, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2, 3, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 2,-1, 2,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 2,-2, 1,     18.0,   0.0,    0.0,     -9.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 6, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-1,-2, 2, 4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 6, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 4, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0),

            /* 431-440 */
            new LuniSolar( 3, 0, 0, 2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0),
            new LuniSolar( 3,-1, 2, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 2, 0, 0,     30.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 4, 0, 2,     24.0,   0.0,    0.0,    -10.0,  0.0,    0.0),
            new LuniSolar( 5, 0, 2,-2, 2,     10.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2, 4, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 2, 2, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2, 4, 2,     17.0,   0.0,    0.0,     -7.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 4, 2,    -24.0,   0.0,    0.0,     10.0,  0.0,    0.0),
            new LuniSolar( 3,-1, 2, 2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0),

            /* 441-450 */
            new LuniSolar( 3, 0, 2, 2, 1,    -24.0,   0.0,    0.0,     11.0,  0.0,    0.0),
            new LuniSolar( 5, 0, 2, 0, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 6, 2,    -13.0,   0.0,    0.0,      5.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 2, 2, 2,    -15.0,   0.0,    0.0,      7.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 1,-1, 1,      0.0,   0.0,-1988.0,      0.0,  0.0,-1679.0),
            new LuniSolar(-1, 0, 1, 0, 3,      0.0,   0.0,  -63.0,      0.0,  0.0,  -27.0),
            new LuniSolar( 0,-2, 2,-2, 3,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0,-1, 0, 1,      0.0,   0.0,    5.0,      0.0,  0.0,    4.0),
            new LuniSolar( 2,-2, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 1, 0, 2,      0.0,   0.0,  364.0,      0.0,  0.0,  176.0),

            /* 451-460 */
            new LuniSolar(-1, 0, 1, 0, 1,      0.0,   0.0,-1044.0,      0.0,  0.0, -891.0),
            new LuniSolar(-1,-1, 2,-1, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar(-2, 2, 0, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 1, 0, 0,      0.0,   0.0,  330.0,      0.0,  0.0,    0.0),
            new LuniSolar(-4, 1, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 1, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 1, 0,-2, 1, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 2,-1,-2, 0, 1,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-4, 0, 2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 461-470 */
            new LuniSolar(-3, 1, 0, 3, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 0,-1, 2, 0,      0.0,   0.0,    5.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 0, 0, 2,      0.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 0, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-1, 0,-2, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-4, 0, 0, 4, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2, 1,-2, 0, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 0,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),

            /* 471-480 */
            new LuniSolar( 0, 0, 1,-1, 0,     -5.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 2, 0, 1, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 2, 0, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 0,-1, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 1,-2, 1,      0.0,   0.0,  -12.0,      0.0,  0.0,  -10.0),
            new LuniSolar( 0, 2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2,-3, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2,-1, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 4,-2, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 4,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),

            /* 481-490 */
            new LuniSolar(-2,-2, 0, 2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar(-2, 0,-2, 4, 0,      0.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 2, 2,-4, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2,-4, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar(-1, 2, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0,-3, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-1, 2, 0, 0, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0,-2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),

            /* 491-500 */
            new LuniSolar( 0, 0, 0,-1, 2,     -8.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 0, 1, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 0,-2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 1, 0,-2, 0, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-3, 1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 1,-2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0, 2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-3,-1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),

            /* 501-510 */
            new LuniSolar( 0, 1, 2,-4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 2,-2, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2,-4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0,-2, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0,-2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2, 0,-2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-4, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0,-1, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0, 0,-2, 0, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0),

            /* 511-520 */
            new LuniSolar(-3, 0, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 0,-2, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2, 0,-2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 0, 0,-4, 2, 0,      8.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2,-1,-2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2,-4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0),

            /* 521-530 */
            new LuniSolar( 0, 1, 4,-4, 4,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 4,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar(-1,-1,-2, 4, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 0,-2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 0,-2, 3, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 3, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 0, 1, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 531-540 */
            new LuniSolar( 1, 1,-2, 2, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 2,-2, 1,     10.0,   0.0,   13.0,      6.0,  0.0,   -5.0),
            new LuniSolar( 0, 0, 1, 0, 2,      0.0,   0.0,   30.0,      0.0,  0.0,   14.0),
            new LuniSolar( 0, 0, 1, 0, 1,      0.0,   0.0, -162.0,      0.0,  0.0, -138.0),
            new LuniSolar( 0, 0, 1, 0, 0,      0.0,   0.0,   75.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 2, 0, 2, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0,-1, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),

            /* 541-550 */
            new LuniSolar( 3, 0, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2,-2, 3,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 2, 0, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2,-3, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 4,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-2,-2, 0, 4, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 0,-2, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),

            /* 551-560 */
            new LuniSolar(-1, 0, 0, 3, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 2,-2, 0, 0, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 2, 0, 1,     -6.0,   0.0,   -3.0,      3.0,  0.0,    1.0),
            new LuniSolar(-1, 0, 1, 2, 1,      0.0,   0.0,   -3.0,      0.0,  0.0,   -2.0),
            new LuniSolar(-1, 1, 0, 3, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 1, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 2, 0, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2, 1, 2, 2, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),

            /* 561-570 */
            new LuniSolar( 2,-2, 2,-2, 2,     -1.0,   0.0,    3.0,      3.0,  0.0,   -1.0),
            new LuniSolar( 1, 1, 0, 1, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 1, 0, 1,      0.0,   0.0,  -13.0,      0.0,  0.0,  -11.0),
            new LuniSolar( 1, 0, 1, 0, 0,      3.0,   0.0,    6.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 2, 0, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 2,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 0,-1, 4,-2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4,-2, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 2,-4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0),

            /* 571-580 */
            new LuniSolar( 2, 2, 2,-2, 2,      8.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 4,-4, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1,-2, 0, 4, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1,-3, 2, 2, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-3, 0,-2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0,-4, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0),

            /* 581-590 */
            new LuniSolar(-2, 1, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-4, 0, 0, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0,-4, 1,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 3, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 0, 4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 2, 0, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 2, 3,      6.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),

            /* 591-600 */
            new LuniSolar(-2, 0, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 0, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 0, 1, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 2,-1, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 0, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 3, 0, 3,      0.0,   0.0,  -26.0,      0.0,  0.0,  -11.0),
            new LuniSolar( 0, 0, 3, 0, 2,      0.0,   0.0,  -10.0,      0.0,  0.0,   -5.0),
            new LuniSolar(-1, 2, 2, 2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 4, 0, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 601-610 */
            new LuniSolar( 1, 2, 2, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 3, 1, 2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 4,-2, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 0, 6, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 0, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 0, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-2,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0,-3, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 2, 3, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),

            /* 611-620 */
            new LuniSolar(-2, 0, 2, 4, 0,     13.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 0, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 4, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 0, 4, 0,    -11.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 1, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 2, 3,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 4, 2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0),

            /* 621-630 */
            new LuniSolar( 2, 1, 0, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 0, 2, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 2, 0, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 1, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 0, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 2, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 3, 0, 3,      0.0,   0.0,   -5.0,      0.0,  0.0,   -2.0),
            new LuniSolar( 1, 1, 2, 1, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0),

            /* 631-640 */
            new LuniSolar( 0, 2, 2, 2, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2, 0, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 4, 1, 2,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar(-1,-1, 0, 6, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar(-3,-1, 2, 6, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 0, 6, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-3, 0, 2, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 0, 4, 0,     12.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 641-650 */
            new LuniSolar(-2, 0, 2, 5, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 1,-2, 2, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 3,-1, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 2, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 3, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar(-1, 1, 2, 4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 0, 1, 2, 3, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 4, 2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 1, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 5, 0, 0, 0, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 651-660 */
            new LuniSolar( 2, 1, 2, 1, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 4, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 3, 1, 2, 0, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 4,-2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar(-2,-1, 2, 6, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 0, 6, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar(-2, 0, 2, 6, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0, 4, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 0, 4, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0),

            /* 661-670 */
            new LuniSolar( 2,-2, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 2, 4, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 1, 0, 2, 3, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0),
            new LuniSolar( 4, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0),
            new LuniSolar( 0, 0, 4, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 4,-1, 2, 0, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 3, 0, 2, 1, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 2, 1, 2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 4, 1, 2, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0),

            /* 671-678 */
            new LuniSolar(-1,-1, 2, 6, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar(-1, 0, 2, 6, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 1,-1, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0),
            new LuniSolar( 1, 1, 2, 4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0),
            new LuniSolar( 3, 1, 2, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0),
            new LuniSolar( 5, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 2,-1, 2, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0),
            new LuniSolar( 2, 0, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0)
        };

        /* ------------------------ */
        /* Planetary nutation model */
        /* ------------------------ */

        /* The units for the sine and cosine coefficients are */
        /* 0.1 microarcsecond                                 */

        static readonly PlanetaryNutation[] xpl =
        {
            /* 1-10 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 0, 1440,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -8, 16,-4,-5, 0, 0, 2,   56,-117,  -42, -40),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 2,  125, -43,    0, -54),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 0,-1, 2, 2,    0,   5,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -4,  8,-1,-5, 0, 0, 2,    3,  -7,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 1,    3,   0,    0,  -2),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0,  3, -8, 3, 0, 0, 0, 0, -114,   0,    0,  61),
            new PlanetaryNutation(-1, 0, 0, 0, 0, 10, -3,  0, 0, 0, 0, 0, 0, -219,  89,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0,-2, 6,-3, 0, 2,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0, -462,1604,    0,   0),

            /* 11-20 */
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -5,  8,-3, 0, 0, 0, 0,   99,   0,    0, -53),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -4,  8,-3, 0, 0, 0, 1,   -3,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -8, 1, 5, 0, 0, 2,    0,   6,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  6,  4, 0, 0, 0, 0, 2,    3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 2,  -12,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 1,   14,-218,  117,   8),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 2,-5, 0, 0, 0,   31,-481, -257, -17),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 0, -491, 128,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0,-2, 5, 0, 0, 0,-3084,5123, 2735,1647),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 1,-1444,2409,-1286,-771),

            /* 21-30 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 2,   11, -24,  -11,  -9),
            new PlanetaryNutation( 2,-1,-1, 0, 0,  0,  3, -7, 0, 0, 0, 0, 0,   26,  -9,    0,   0),
            new PlanetaryNutation( 1, 0,-2, 0, 0, 19,-21,  3, 0, 0, 0, 0, 0,  103, -60,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  2, -4,  0,-3, 0, 0, 0, 0,    0, -13,   -7,   0),
            new PlanetaryNutation( 1, 0,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -26, -29,  -16,  14),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0,-4,10, 0, 0, 0,    9, -27,  -14,  -5),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  0,  2,  0, 0,-5, 0, 0, 0,   12,   0,    0,  -6),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -7,  4, 0, 0, 0, 0, 0,   -7,   0,    0,   0),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0, 1,-1, 0, 0, 0,    0,  24,    0,   0),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,  284,   0,    0,-151),

            /* 31-40 */
            new PlanetaryNutation(-1, 0, 0, 0, 0, 18,-16,  0, 0, 0, 0, 0, 0,  226, 101,    0,   0),
            new PlanetaryNutation(-2, 1, 1, 2, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -8,   -2,   0),
            new PlanetaryNutation(-1, 1,-1, 1, 0, 18,-17,  0, 0, 0, 0, 0, 0,    0,  -6,   -3,   0),
            new PlanetaryNutation(-1, 0, 1, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 2,  -41, 175,   76,  17),
            new PlanetaryNutation( 0, 2,-2, 2, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,  15,    6,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 1,  425, 212, -133, 269),
            new PlanetaryNutation( 0, 1,-1, 1, 0, -8, 12,  0, 0, 0, 0, 0, 0, 1200, 598,  319,-641),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 0,  235, 334,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  8,-14,  0, 0, 0, 0, 0, 0,   11, -12,   -7,  -6),

            /* 41-50 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 1,    5,  -6,    3,   3),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  0,  2,  0,-4, 5, 0, 0, 0,   -5,   0,    0,   3),
            new PlanetaryNutation(-2, 0, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0,-3, 1, 0, 0, 0,   15,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  3, -5,  0, 2, 0, 0, 0, 0,   13,   0,    0,  -7),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0,-4, 3, 0, 0, 0,   -6,  -9,    0,   0),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,  266, -78,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -1,  2, 0, 0, 0, 0, 0, -460,-435, -232, 246),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -2,  2, 0, 0, 0, 0, 0,    0,  15,    7,   0),
            new PlanetaryNutation(-1, 1, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2),

            /* 51-60 */
            new PlanetaryNutation(-1, 0, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0, 131,    0,   0),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0,-2,-2, 0, 0, 0,    4,   0,    0,   0),
            new PlanetaryNutation(-2, 2, 0, 2, 0,  0, -5,  9, 0, 0, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0, 0,-1, 0, 0,    0,   4,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 1, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 2, 0,  -17, -19,  -10,   9),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 1,   -9, -11,    6,  -5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 2,   -6,   0,    0,   3),
            new PlanetaryNutation(-1, 0, 1, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -16,   8,    0,   0),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0, 0, 2, 0, 0, 0,    0,   3,    0,   0),

            /* 61-70 */
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -1,  0, 0, 2, 0, 0, 0,   11,  24,   11,  -5),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -9, 17, 0, 0, 0, 0, 0,   -3,  -4,   -2,   1),
            new PlanetaryNutation( 0, 0, 0, 2, 0, -3,  5,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0,-1, 2, 0, 0, 0,    0,  -8,   -4,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 1,-2, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation( 1, 0,-2, 0, 0, 17,-16,  0,-2, 0, 0, 0, 0,    0,   5,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 1,-3, 0, 0, 0,    0,   3,    2,   0),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  0,  5, -6, 0, 0, 0, 0, 0,   -6,   4,    2,   3),
            new PlanetaryNutation( 0,-2, 2, 0, 0,  0,  9,-13, 0, 0, 0, 0, 0,   -3,  -5,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -1,  0, 0, 1, 0, 0, 0,   -5,   0,    0,   2),

            /* 71-80 */
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  0,  0, 0, 1, 0, 0, 0,    4,  24,   13,  -2),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0, 0, 1, 0, 0, 0,  -42,  20,    0,   0),
            new PlanetaryNutation( 0,-2, 2, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,  -10, 233,    0,   0),
            new PlanetaryNutation( 0,-1, 1, 1, 0,  5, -7,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   78, -18,    0,   0),
            new PlanetaryNutation( 2, 1,-3, 1, 0, -6,  7,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 2, 0,  0,  0,  0, 1, 0, 0, 0, 0,    0,  -3,   -1,   0),
            new PlanetaryNutation( 0,-1, 1, 1, 0,  0,  1,  0, 1, 0, 0, 0, 0,    0,  -4,   -2,   1),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 2, 0, 0,    0,  -8,   -4,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 1,    0,  -5,    3,   0),

            /* 81-90 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 2,   -7,   0,    0,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 2,  -14,   8,    3,   6),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 1,    0,   8,   -4,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -9, 15, 0, 0, 0, 0, 0,    0,  19,   10,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   45, -22,    0,   0),
            new PlanetaryNutation( 1,-1,-1, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 2, 0,-2, 0, 0,  2, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0,-5, 5, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation( 2, 0,-2, 1, 0,  0, -6,  8, 0, 0, 0, 0, 0,    3,   5,    3,  -2),
            new PlanetaryNutation( 2, 0,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   89, -16,   -9, -48),

            /* 91-100 */
            new PlanetaryNutation(-2, 1, 1, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation(-2, 1, 1, 1, 0,  0,  1,  0,-3, 0, 0, 0, 0,   -3,   7,    4,   2),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0, -349, -62,    0,   0),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  6, -8, 0, 0, 0, 0, 0,  -15,  22,    0,   0),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0,-1,-5, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,  -53,   0,    0,   0),
            new PlanetaryNutation(-1, 1, 1, 1, 0,-20, 20,  0, 0, 0, 0, 0, 0,    5,   0,    0,  -3),
            new PlanetaryNutation( 1, 0,-2, 0, 0, 20,-21,  0, 0, 0, 0, 0, 0,    0,  -8,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  8,-15, 0, 0, 0, 0, 0,   15,  -7,   -4,  -8),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0,-10, 15, 0, 0, 0, 0, 0,   -3,   0,    0,   1),

            /* 101-110 */
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,  -21, -78,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  0,  0, 1, 0, 0, 0, 0,   20, -70,  -37, -11),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,   6,    3,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0,-2, 4, 0, 0, 0,    5,   3,    2,  -2),
            new PlanetaryNutation( 2, 0,-2, 1, 0, -6,  8,  0, 0, 0, 0, 0, 0,  -17,  -4,   -2,   9),
            new PlanetaryNutation( 0,-2, 2, 1, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   6,    3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0,-1, 0, 0, 1,   32,  15,   -8,  17),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0,-1, 0, 0, 0,  174,  84,   45, -93),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 0,   11,  56,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0, 1, 0, 0, 0,  -66, -12,   -6,  35),

            /* 111-120 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 1,   47,   8,    4, -25),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 2,    0,   8,    4,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -9, 13, 0, 0, 0, 0, 0,   10, -22,  -12,  -5),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  7,-13, 0, 0, 0, 0, 0,   -3,   0,    0,   2),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,  -24,  12,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  9,-17, 0, 0, 0, 0, 0,    5,  -6,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -9, 17, 0, 0, 0, 0, 2,    3,   0,    0,  -2),
            new PlanetaryNutation( 1, 0,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    4,   3,    1,  -2),
            new PlanetaryNutation( 1, 0,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  29,   15,   0),
            new PlanetaryNutation( 0, 0, 0, 2, 0,  0, -1,  2, 0, 0, 0, 0, 0,   -5,  -4,   -2,   2),

            /* 121-130 */
            new PlanetaryNutation( 0,-1, 1, 1, 0,  0,  0,  2, 0, 0, 0, 0, 0,    8,  -3,   -1,  -5),
            new PlanetaryNutation( 0,-2, 2, 0, 1,  0, -2,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -5,  0, 2, 0, 0, 0, 0,   10,   0,    0,   0),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  0,  2,  0,-3, 1, 0, 0, 0,    3,   0,    0,  -2),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  3, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   3),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  8,-13,  0, 0, 0, 0, 0, 0,   46,  66,   35, -25),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  8,-12,  0, 0, 0, 0, 0, 0,  -14,   7,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0),
            new PlanetaryNutation(-1, 0, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,   -5,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 0, 1, 0, 18,-16,  0, 0, 0, 0, 0, 0,  -68, -34,  -18,  36),

            /* 131-140 */
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0,-1, 1, 0, 0, 0,    0,  14,    7,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  3, -7,  4, 0, 0, 0, 0, 0,   10,  -6,   -3,  -5),
            new PlanetaryNutation(-2, 1, 1, 1, 0,  0, -3,  7, 0, 0, 0, 0, 0,   -5,  -4,   -2,   3),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -1,  0,-2, 5, 0, 0, 0,   -3,   5,    2,   1),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  0,  0,-2, 5, 0, 0, 0,   76,  17,    9, -41),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,   84, 298,  159, -45),
            new PlanetaryNutation( 1, 0, 0, 1, 0,-10,  3,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2),
            new PlanetaryNutation(-1, 0, 0, 1, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,  -82, 292,  156,  44),

            /* 141-150 */
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  0,  0, 2,-5, 0, 0, 0,  -73,  17,    9,  39),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,   -9, -16,    0,   0),
            new PlanetaryNutation( 2,-1,-1, 1, 0,  0,  3, -7, 0, 0, 0, 0, 0,    3,   0,   -1,  -2),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0, 0,-5, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -3,  7, -4, 0, 0, 0, 0, 0,   -9,  -5,   -3,   5),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0, -439,   0,    0,   0),
            new PlanetaryNutation( 1, 0, 0, 1, 0,-18, 16,  0, 0, 0, 0, 0, 0,   57, -28,  -15, -30),
            new PlanetaryNutation(-2, 1, 1, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -6,   -3,   0),
            new PlanetaryNutation( 0, 1,-1, 2, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -8, 13,  0, 0, 0, 0, 0, 0,  -40,  57,   30,  21),

            /* 151-160 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 1,   23,   7,    3, -13),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0,  0, -2, 0, 0, 0, 0, 0,  273,  80,   43,-146),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 0, -449, 430,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,   -8, -47,  -25,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  2, 0, 0, 0, 0, 1,    6,  47,   25,  -3),
            new PlanetaryNutation(-1, 0, 1, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  23,   13,   0),
            new PlanetaryNutation(-1, 0, 1, 1, 0,  0,  3, -4, 0, 0, 0, 0, 0,   -3,   0,    0,   2),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0,-2, 0, 0, 0,    3,  -4,   -2,  -2),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0, 2, 0, 0, 0,  -48,-110,  -59,  26),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 1,   51, 114,   61, -27),

            /* 161-170 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 2, -133,   0,    0,  57),
            new PlanetaryNutation( 0, 1,-1, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -3,  5,  0, 0, 0, 0, 0, 0,  -21,  -6,   -3,  11),
            new PlanetaryNutation( 0, 1,-1, 2, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  -3,   -1,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -2,  4, 0, 0, 0, 0, 0,  -11, -21,  -11,   6),
            new PlanetaryNutation( 0, 2,-2, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,  -18,-436, -233,   9),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   35,  -7,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  5, -8,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  6, -8,  0, 0, 0, 0, 0, 0,   11,  -3,   -1,  -6),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -8, 15, 0, 0, 0, 0, 0,   -5,  -3,   -1,   3),

            /* 171-180 */
            new PlanetaryNutation(-2, 0, 2, 1, 0,  0,  2,  0,-3, 0, 0, 0, 0,  -53,  -9,   -5,  28),
            new PlanetaryNutation(-2, 0, 2, 1, 0,  0,  6, -8, 0, 0, 0, 0, 0,    0,   3,    2,   1),
            new PlanetaryNutation( 1, 0,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    4,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 3,-5, 0, 0, 0,    0,  -4,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0,-1, 0, 0, 0, 0,  -50, 194,  103,  27),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0,-1, 0, 0, 0, 1,  -13,  52,   28,   7),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 0,  -91, 248,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    6,  49,   26,  -3),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   -6, -47,  -25,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    0,   5,    3,   0),

            /* 181-190 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 2,   52,  23,   10, -23),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -1,  0, 0,-1, 0, 0, 0,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  0,  0, 0,-1, 0, 0, 0,    0,   5,    3,   0),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,   -4,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -7, 13, 0, 0, 0, 0, 2,   -4,   8,    3,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  7,-13, 0, 0, 0, 0, 0,   10,   0,    0,   0),
            new PlanetaryNutation( 2, 0,-2, 1, 0,  0, -5,  6, 0, 0, 0, 0, 0,    3,   0,    0,  -2),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -8, 11, 0, 0, 0, 0, 0,    0,   8,    4,   0),
            new PlanetaryNutation( 0, 2,-2, 1,-1,  0,  2,  0, 0, 0, 0, 0, 0,    0,   8,    4,   1),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,   -4,   0,    0,   0),

            /* 191-200 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2,-2, 0, 0, 0,   -4,   0,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 0, 3, 0, 0, 0,   -8,   4,    2,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 1,    8,  -4,   -2,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 2,    0,  15,    7,   0),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0, -138,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,  -7,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -7,   -3,   0),
            new PlanetaryNutation( 2, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   54,   0,    0, -29),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,  10,    4,   0),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0,  0, -2, 0, 0, 0, 0, 0,   -7,   0,    0,   3),

            /* 201-210 */
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  1, -2, 0, 0, 0, 0, 0,  -37,  35,   19,  20),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,    0,   4,    0,   0),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -4,   9,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -2,  0, 0, 2, 0, 0, 0,    8,   0,    0,  -4),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  3, -6,  0, 0, 0, 0, 0, 0,   -9, -14,   -8,   5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 1,   -3,  -9,   -5,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 0, -145,  47,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,  -10,  40,   21,   5),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 1,   11, -49,  -26,  -7),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,-2150,   0,    0, 932),

            /* 211-220 */
            new PlanetaryNutation( 0, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,  -12,   0,    0,   5),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,   85,   0,    0, -37),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 1,    4,   0,    0,  -2),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0,  1, -4, 0, 0, 0, 0, 0,    3,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 0,  -86, 153,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -6,   9,    5,   3),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    9, -13,   -7,  -5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -8,  12,    6,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 2,  -51,   0,    0,  22),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,  -11,-268, -116,   5),

            /* 221-230 */
            new PlanetaryNutation( 0, 2,-2, 2, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,  12,    5,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,    0,   7,    3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   31,   6,    3, -17),
            new PlanetaryNutation( 0, 1,-1, 1, 0, -5,  7,  0, 0, 0, 0, 0, 0,  140,  27,   14, -75),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   57,  11,    6, -30),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -8,  0, 0, 0, 0, 0, 0,  -14, -39,    0,   0),
            new PlanetaryNutation( 0, 1,-1, 2, 0,  0, -1,  0,-1, 0, 0, 0, 0,    0,  -6,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  0,  0,-1, 0, 0, 0, 0,    4,  15,    8,  -2),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    0,   4,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -2,  0, 1, 0, 0, 0, 0,   -3,   0,    0,   1),

            /* 231-240 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -6, 11, 0, 0, 0, 0, 2,    0,  11,    5,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,-11, 0, 0, 0, 0, 0,    9,   6,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0,-1,  0,  4,  0, 0, 0, 0, 0, 2,   -4,  10,    4,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 1,  0, -4,  0, 0, 0, 0, 0, 0,    5,   3,    0,   0),
            new PlanetaryNutation( 2, 0,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,   16,   0,    0,  -9),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -7,  9, 0, 0, 0, 0, 0,    0,   3,    2,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 4,-5, 0, 0, 2,    7,   0,    0,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 0,  -25,  22,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,   42, 223,  119, -22),

            /* 241-250 */
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -27,-143,  -77,  14),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,    9,  49,   26,  -5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 2,-1166,   0,    0, 505),
            new PlanetaryNutation( 0, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -5,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 5, 0, 0, 2,   -6,   0,    0,   3),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -8,   0,    1,   4),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,  117,   0,    0, -63),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  2, -4, 0, 0, 0, 0, 0,   -4,   8,    4,   2),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -4,  4, 0, 0, 0, 0, 0,    3,   0,    0,  -2),

            /* 251-260 */
            new PlanetaryNutation( 0, 1,-1, 2, 0, -5,  7,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -6, 0, 0, 0, 0, 0,    0,  31,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -5,   0,    1,   3),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -4,  6, 0, 0, 0, 0, 0,    4,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 2,  -24, -13,   -6,  10),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0, -32,  -17,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 2,    8,  12,    5,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 1,    3,   0,    0,  -1),

            /* 261-270 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -9, 0, 0, 0, 0, 0,    7,  13,    0,   0),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0,   -3,  16,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   50,   0,    0, -27),
            new PlanetaryNutation(-2, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -5,   -3,   0),
            new PlanetaryNutation( 0,-2, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 1,    0,   5,    3,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 2,   24,   5,    2, -11),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 2,    5, -11,   -5,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 1,   30,  -3,   -2, -16),
            new PlanetaryNutation( 0, 1,-1, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   18,   0,    0,  -9),

            /* 271-280 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 0,    8, 614,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 1,    3,  -3,   -1,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    6,  17,    9,  -3),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0, -1,  0, 3, 0, 0, 0, 0,   -3,  -9,   -5,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    0,   6,    3,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 2, -127,  21,    9,  55),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -8, 0, 0, 0, 0, 0,    3,   5,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -4,  8, 0, 0, 0, 0, 2,   -6, -10,   -4,   3),
            new PlanetaryNutation( 0,-2, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    5,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 2,   16,   9,    4,  -7),

            /* 281-290 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 1,    3,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -7, 0, 0, 0, 0, 0,    0,  22,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,  19,   10,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,    7,   0,    0,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -5, 10, 0, 0, 0, 0, 2,    0,  -5,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 4, 0, 0, 0, 2,   -9,   3,    1,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 2,   17,   0,    0,  -7),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 1,    0,  -3,   -2,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -5, 0, 0, 0, 0, 0,  -20,  34,    0,   0),

            /* 291-300 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 1,  -10,   0,    1,   5),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  1, -3,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 0,   22, -87,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 2,   -3,  -6,   -2,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 2,  -16,  -3,   -1,   7),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0),
            new PlanetaryNutation( 0,-2, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2, -3, 0, 0, 0, 0, 0,  -68,  39,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0, -4,  4,  0, 0, 0, 0, 0, 0,   27,   0,    0, -14),

            /* 301-310 */
            new PlanetaryNutation( 0,-1, 1, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1, -1, 0, 0, 0, 0, 0,  -25,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 1,  -12,  -3,   -2,   6),
            new PlanetaryNutation( 0, 1,-1, 1, 0, -4,  6,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 2,    3,  66,   29,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 2,  490,   0,    0,-213),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,  -22,  93,   49,  12),
            new PlanetaryNutation( 0, 1,-1, 1, 0, -4,  5,  0, 0, 0, 0, 0, 0,   -7,  28,   15,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,   -3,  13,    7,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -6,  0, 0, 0, 0, 0, 0,  -46,  14,    0,   0),

            /* 311-320 */
            new PlanetaryNutation(-2, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  1, 0, 0, 0, 0, 0,    2,   1,    0,   0),
            new PlanetaryNutation( 0,-1, 1, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,  -28,   0,    0,  15),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1, -3, 0, 0, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  3, 0, 0, 0, 0, 2,  -11,   0,    0,   5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -7, 12, 0, 0, 0, 0, 2,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 1,   25, 106,   57, -13),

            /* 321-330 */
            new PlanetaryNutation( 0, 1,-1, 1, 0, -1,  0,  0, 0, 0, 0, 0, 0,    5,  21,   11,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0, 1485,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 1,   -7, -32,  -17,   4),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  1, -2,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  5, 0, 0, 0, 0, 2,   -6,  -3,   -2,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 4, 0, 0, 0, 2,   30,  -6,   -2, -13),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0,-4, 0, 0, 0, 0,   -4,   4,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -1,  1,  0, 0, 0, 0, 0, 0,  -19,   0,    0,  10),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 2,    0,   4,    2,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 0,    0,   3,    0,   0),

            /* 331-340 */
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -3,  0, 3, 0, 0, 0, 0,    4,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  7, 0, 0, 0, 0, 2,    0,  -3,   -1,   0),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -5,  8, 0, 0, 0, 0, 2,    5,   3,    1,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -8, 0, 0, 0, 0, 0,    0,  11,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 2,  118,   0,    0, -52),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 1,    0,  -5,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,  -28,  36,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2, -4,  0, 0, 0, 0, 0, 0,    5,  -5,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 1,   14, -59,  -31,  -8),

            /* 341-350 */
            new PlanetaryNutation( 0, 1,-1, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   9,    5,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 2, -458,   0,    0, 198),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 2,    0, -45,  -20,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 1,    9,   0,    0,  -5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  6, -9,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -4,   -2,  -1),
            new PlanetaryNutation( 0, 2,-2, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   11,   0,    0,  -6),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -4,  6, 0, 0, 0, 0, 2,    6,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -6, 0, 0, 0, 0, 0,  -16,  23,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,   -2,   0),

            /* 351-360 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 2, 0, 0, 0, 2,   -5,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0, -166, 269,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   15,   0,    0,  -8),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  9,  0, 0, 0, 0, 0, 2,   10,   0,    0,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -78,  45,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 2,    0,  -5,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 1,    7,   0,    0,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,   -5, 328,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -2),

            /* 361-370 */
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 0,-3, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 1,-5, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 1,    0,  -4,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,-1223, -26,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 1,    0,   7,    3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0,-3, 5, 0, 0, 0,    3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -6,  20,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0, -368,   0,    0,   0),

            /* 371-380 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,  -75,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   11,   0,    0,  -6),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,    3,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 14,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,  -13, -30,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 0,   21,   3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 2,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    8, -27,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -8, 3, 0, 0, 0, 0,  -19, -11,    0,   0),

            /* 381-390 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  8,-3, 0, 0, 0, 2,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0,-2, 5, 0, 0, 2,    0,   5,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 1,-2, 0, 0, 0,   -1,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 0, 1, 0, 0, 2,  -14,   0,    0,   6),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,    6,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 2,  -74,   0,    0,  32),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 0, 2, 0, 0, 2,    0,  -3,   -1,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0, -5,  5,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2),

            /* 391-400 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,    8,  11,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 1,    0,   3,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 2, -262,   0,    0, 114),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 1,   -7,   0,    0,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 2,    0, -27,  -12,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  4, 0, 0, 0, 0, 2,  -19,  -8,   -4,   8),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 2,  202,   0,    0, -87),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 1,   -8,  35,   19,   5),
            new PlanetaryNutation( 0, 1,-1, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,   4,    2,   0),

            /* 401-410 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   16,  -5,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    5,   0,    0,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,  -3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0,-1,  0,  3,  0, 0, 0, 0, 0, 2,    1,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 2, 0, 0, 0, 2,  -35, -48,  -21,  15),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  6, 0, 0, 0, 0, 2,   -3,  -5,   -2,   1),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -6,  9, 0, 0, 0, 0, 2,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6, -9, 0, 0, 0, 0, 0,    0,  -5,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  2,  0, 0, 0, 0, 0, 1,   12,  55,   29,  -6),

            /* 411-420 */
            new PlanetaryNutation( 0, 1,-1, 1, 0, -2,  1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0, -598,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 1,   -3, -13,   -7,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  0, 3, 0, 0, 0, 2,   -5,  -7,   -3,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -5,  7, 0, 0, 0, 0, 2,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -7, 0, 0, 0, 0, 0,    5,  -7,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -5, 0, 0, 0, 0, 0,   16,  -6,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1, -3,  0, 0, 0, 0, 0, 0,    8,  -3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 1,    8, -31,  -16,  -4),

            /* 421-430 */
            new PlanetaryNutation( 0, 1,-1, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 2,  113,   0,    0, -49),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 2,    0, -24,  -10,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 1,    4,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -3, 0, 0, 0, 0, 0,   27,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  8,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 2,    0,  -4,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 1,    5,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  1, 0, 0, 0, 0, 2,  -13,   0,    0,   6),

            /* 431-440 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  3, 0, 0, 0, 0, 2,  -18, -10,   -4,   8),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,   -4, -28,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 2,   -5,   6,    3,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -9, 13,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  5, 0, 0, 0, 0, 2,   -5,  -9,   -4,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  0, 4, 0, 0, 0, 2,   17,   0,    0,  -7),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-4, 0, 0, 0, 0,   11,   4,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  7, 0, 0, 0, 0, 2,    0,  -6,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   83,  15,    0,   0),

            /* 441-450 */
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 2,    0,-114,  -49,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 2,  117,   0,    0, -51),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 1,   -5,  19,   10,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  9, 0, 0, 0, 0, 2,    0,  -3,   -1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,    3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 2,    0,  -6,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  393,   3,    0,   0),

            /* 451-460 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 1,   -4,  21,   11,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 2,   -6,   0,   -1,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5, 10,  0, 0, 0, 0, 0, 2,   -3,   8,    4,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,    8,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 2,   18, -29,  -13,  -8),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  3,  0, 0, 0, 0, 0, 1,    8,  34,   18,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   89,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 1,    3,  12,    6,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 2,   54, -15,   -7, -24),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0,-3, 0, 0, 0,    0,   3,    0,   0),

            /* 461-470 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -5, 13, 0, 0, 0, 0, 2,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 0,    0,  35,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 2, -154, -30,  -13,  67),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   15,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 1,    0,   4,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 0,    0,   9,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 2,   80, -71,  -31, -35),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0,-1, 0, 0, 2,    0, -20,   -9,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -6, 15, 0, 0, 0, 0, 2,   11,   5,    2,  -5),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 15,  0, 0, 0, 0, 0, 2,   61, -96,  -42, -27),

            /* 471-480 */
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  9, -4, 0, 0, 0, 0, 2,   14,   9,    4,  -6),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 2,-5, 0, 0, 2,  -11,  -6,   -3,   5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  8,-1,-5, 0, 0, 2,    0,  -3,   -1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6, -8, 3, 0, 0, 0, 2,  123,-415, -180, -53),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,    0,   0,    0, -35),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    7, -32,  -17,  -4),
            new PlanetaryNutation( 0, 1,-1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -9,   -5,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    0,  -4,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 2,  -89,   0,    0,  38),

            /* 481-490 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -6, 16,-4,-5, 0, 0, 2,    0, -86,  -19,  -6),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2,    0,   0,  -19,   6),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2, -123,-416, -180,  53),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6, -8, 1, 5, 0, 0, 2,    0,  -3,   -1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-2, 5, 0, 0, 2,   12,  -6,   -3,  -5),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -5,  4, 0, 0, 0, 0, 2,  -13,   9,    4,   6),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,    0, -15,   -7,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,  -62, -97,  -42,  27),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, 11,  0, 0, 0, 0, 0, 2,  -11,   5,    2,   5),

            /* 491-500 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 0, 1, 0, 0, 2,    0, -19,   -8,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -3,  0, 2, 0, 0, 0, 2,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   4,    2,   0),
            new PlanetaryNutation( 0, 1,-1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   4,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  1,  2, 0, 0, 0, 0, 2,  -85, -70,  -31,  37),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 1, 0, 0, 0, 2,  163, -12,   -5, -72),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  7,  0, 0, 0, 0, 0, 2,  -63, -16,   -7,  28),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  4, 0, 0, 0, 0, 2,  -21, -32,  -14,   9),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 2,    0,  -3,   -1,   0),

            /* 501-510 */
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   8,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 2,    3,  10,    4,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0, 2, 0, 0, 0, 2,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  6, 0, 0, 0, 0, 2,    0,  -7,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  7, -9, 0, 0, 0, 0, 2,    0,  -4,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 0,    6,  19,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 2,    5,-173,  -75,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6, -7, 0, 0, 0, 0, 2,    0,  -7,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -5, 0, 0, 0, 0, 2,    7, -12,   -5,  -3),

            /* 511-520 */
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 1,   -3,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 2,    3,  -4,   -2,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 2,   74,   0,    0, -32),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 1,   -3,  12,    6,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -3, 0, 0, 0, 0, 2,   26, -14,   -6, -11),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3, -1, 0, 0, 0, 0, 2,   19,   0,    0,  -8),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -4,  4,  0, 0, 0, 0, 0, 1,    6,  24,   13,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   83,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 1,    0, -10,   -5,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 2,   11,  -3,   -1,  -5),

            /* 521-530 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  1, 0, 0, 0, 0, 2,    3,   0,    1,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -3,  0, 5, 0, 0, 0, 2,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 1,    5, -23,  -12,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 2, -339,   0,    0, 147),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -9, 12,  0, 0, 0, 0, 0, 2,    0, -10,   -5,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0,-4, 0, 0, 0, 0,    5,   0,    0,   0),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  7, -8, 0, 0, 0, 0, 2,    0,  -4,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 0,   18,  -3,    0,   0),

            /* 531-540 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 2,    9, -11,   -5,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -2,  6,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6,  7,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  6, -7,  0, 0, 0, 0, 0, 0,    0,   9,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6, -6, 0, 0, 0, 0, 2,    6,  -9,   -4,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 0,   -4, -12,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 2,   67, -91,  -39, -29),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -4, 0, 0, 0, 0, 2,   30, -18,   -8, -13),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 0,    0,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 2,    0,-114,  -50,   0),

            /* 541-550 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,    0,   0,    0,  23),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,  517,  16,    7,-224),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0, 0,-2, 0, 0, 2,    0,  -7,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4, -2, 0, 0, 0, 0, 2,  143,  -3,   -1, -62),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0, 0,-1, 0, 0, 2,   29,   0,    0, -13),
            new PlanetaryNutation( 0, 2,-2, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -4,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 16,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0, 2,-5, 0, 0, 2,    5,  12,    5,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  7, -8, 3, 0, 0, 0, 2,  -25,   0,    0,  11),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -5, 16,-4,-5, 0, 0, 2,   -3,   0,    0,   1),

            /* 551-560 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0, -1,  8,-3, 0, 0, 0, 2,  -22,  12,    5,  10),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,   50,   0,    0, -22),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  2, 0, 0, 0, 0, 2,   -4,   4,    2,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  3,  0, 1, 0, 0, 0, 2,   -5, -11,   -5,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -3,  8,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -5,  5,  0, 0, 0, 0, 0, 1,    4,  17,    9,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 0,   59,   0,    0,   0),

            /* 561-570 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 1,    0,  -4,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 1,    4, -15,   -8,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 2,  370,  -8,    0,-160),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   0,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6, -5, 0, 0, 0, 0, 2,   -6,   3,    1,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  7, -8,  0, 0, 0, 0, 0, 0,    0,   6,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -3, 0, 0, 0, 0, 2,  -10,   0,    0,   4),

            /* 571-580 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -3,  0, 0, 0, 0, 0, 2,    0,   9,    4,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  1,  2,  0, 0, 0, 0, 0, 2,    4,  17,    7,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 2,   34,   0,    0, -15),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4,  0,-4, 0, 0, 0, 2,   -5,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4,  0,-3, 0, 0, 0, 2,  -37,  -7,   -3,  16),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -6,  6,  0, 0, 0, 0, 0, 1,    3,  13,    7,  -2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 0,   40,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4,  0,-2, 0, 0, 0, 2, -184,  -3,   -1,  80),

            /* 581-590 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6, -4, 0, 0, 0, 0, 2,   -3,   0,    0,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 1,    0, -10,   -6,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 2,   31,  -6,    0, -13),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4,  0,-1, 0, 0, 0, 2,   -3, -32,  -14,   1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4,  0, 0,-2, 0, 0, 2,   -7,   0,    0,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5, -2, 0, 0, 0, 0, 2,    0,  -8,   -4,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  4,  0, 0, 0, 0, 0, 0,    3,  -4,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  8, -9,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -4,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0),

            /* 591-600 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 2,   19, -23,  -10,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   0,    0, -10),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   3,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -7,  7,  0, 0, 0, 0, 0, 1,    0,   9,    5,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  7, -7,  0, 0, 0, 0, 0, 0,   28,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 1,    0,  -7,   -4,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 2,    8,  -4,    0,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   0,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5,  0,-4, 0, 0, 0, 2,   -3,   0,    0,   1),

            /* 601-610 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5,  0,-3, 0, 0, 0, 2,   -9,   0,    1,   4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  5,  0,-2, 0, 0, 0, 2,    3,  12,    5,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  3,  0,  0, 0, 0, 0, 0, 2,   17,  -3,   -1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -8,  8,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  8, -8,  0, 0, 0, 0, 0, 0,   19,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 1,    0,  -5,   -3,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 2,   14,  -3,    0,  -1),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,   -1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,    0,  -5),
            new PlanetaryNutation( 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0),

            /* 611-620 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  9, -9,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  6, -4,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    2,   9,    4,   3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    0,   0,    0,  -4),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    8,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   4,    2,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    6,   0,    0,  -3),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    5,   0,    0,  -2),

            /* 621-630 */
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 2,    3,   0,    0,  -1),
            new PlanetaryNutation( 1, 0,-2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 1, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0),
            new PlanetaryNutation( 1, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    7,   0,    0,   0),
            new PlanetaryNutation( 1, 0,-2, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    6,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0),
            new PlanetaryNutation( 1, 0,-2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0),
            new PlanetaryNutation(-2, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    5,   0,    0,   0),

            /* 631-640 */
            new PlanetaryNutation(-1, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0),
            new PlanetaryNutation( 1,-1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0),
            new PlanetaryNutation(-1, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   13,   0,    0,   0),
            new PlanetaryNutation(-2, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   21,  11,    0,   0),
            new PlanetaryNutation( 1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0),
            new PlanetaryNutation(-1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,  -5,   -2,   0),
            new PlanetaryNutation( 1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0),

            /* 641-650 */
            new PlanetaryNutation(-1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0),
            new PlanetaryNutation(-1, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2),
            new PlanetaryNutation( 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   20,  10,    0,   0),
            new PlanetaryNutation(-1, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  -34,   0,    0,   0),
            new PlanetaryNutation(-1, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,  -19,   0,    0,   0),
            new PlanetaryNutation( 1, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,    3,   0,    0,  -2),
            new PlanetaryNutation( 1, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1),
            new PlanetaryNutation( 1, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -6,   0,    0,   3),
            new PlanetaryNutation( 1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0),
            new PlanetaryNutation( 1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    3,   0,    0,   0),

            /* 651-660 */
            new PlanetaryNutation( 0, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0),
            new PlanetaryNutation( 0, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0),
            new PlanetaryNutation( 0, 2, 0, 2, 0, -2,  2,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 2, 0, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    6,   0,    0,  -3),
            new PlanetaryNutation( 0, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   3),
            new PlanetaryNutation( 0, 2, 0, 2, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0),
            new PlanetaryNutation( 0, 1, 1, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -3,   -2,   0),
            new PlanetaryNutation( 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  126, -63,  -27, -55),
            new PlanetaryNutation(-1, 2, 0, 2, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    1,   2),

            /* 661-670 */
            new PlanetaryNutation( 0, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,   -3,  28,   15,   2),
            new PlanetaryNutation( 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    5,   0,    1,  -2),
            new PlanetaryNutation( 0, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   9,    4,   1),
            new PlanetaryNutation( 0, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   9,    4,  -1),
            new PlanetaryNutation(-1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0, -126, -63,  -27,  55),
            new PlanetaryNutation( 2, 2,-2, 2, 0,  0, -2,  0, 3, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 1, 2, 0, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   21, -11,   -6, -11),
            new PlanetaryNutation( 0, 1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0),
            new PlanetaryNutation(-1, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -21, -11,   -6,  11),
            new PlanetaryNutation(-2, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   1),

            /* 671-680 */
            new PlanetaryNutation( 0, 2, 0, 2, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0),
            new PlanetaryNutation( 0, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    8,   0,    0,  -4),
            new PlanetaryNutation( 0, 2, 0, 2, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -6,   0,    0,   3),
            new PlanetaryNutation( 0, 2, 0, 2, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1),
            new PlanetaryNutation(-1, 2, 2, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 1, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1),
            new PlanetaryNutation(-1, 2, 2, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -5,   0,    0,   2),
            new PlanetaryNutation( 2, 2, 0, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   24, -12,   -5, -11),
            new PlanetaryNutation( 1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   3,    1,   0),
            new PlanetaryNutation( 1, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   3,    1,   0),

            /* 681-687 */
            new PlanetaryNutation( 1, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0),
            new PlanetaryNutation( 0, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -24, -12,   -5,  10),
            new PlanetaryNutation( 2, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    4,   0,   -1,  -2),
            new PlanetaryNutation(-1, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   13,   0,    0,  -6),
            new PlanetaryNutation(-1, 2, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    7,   0,    0,  -3),
            new PlanetaryNutation( 1, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1),
            new PlanetaryNutation( 0, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,    3,   0,    0,  -1)
        };

        /// <summary>
        /// Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
        /// with free core nutation omitted).
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation, luni-solar + planetary longitude</param>
        /// <param name="deps">nutation, luni-solar + planetary obliquity</param>
        public static void iauNut00a(double date1, double date2, out double dpsi, out double deps)
        {
            int i;
            double t, el, elp, f, d, om, arg, dp, de, sarg, carg,
                   al, af, ad, aom, alme, alve, alea, alma,
                   alju, alsa, alur, alne, apa, dpsils, depsls,
                   dpsipl, depspl;

            /* Units of 0.1 microarcsecond to radians */
            double U2R = DAS2R / 1e7;


            /* Number of terms in the luni-solar nutation model */
            //const int NLS = (int)(sizeof xls / sizeof xls[0]);
            int NLS = Nut00axls.Length;

            /* Number of terms in the planetary nutation model */
            //const int NPL = (int)(sizeof xpl / sizeof xpl[0]);
            int NPL = xpl.Length;

            /* ------------------------------------------------------------------ */

            /* Interval between fundamental date J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* ------------------- */
            /* LUNI-SOLAR NUTATION */
            /* ------------------- */

            /* Fundamental (Delaunay) arguments */

            /* Mean anomaly of the Moon (IERS 2003). */
            el = iauFal03(t);

            /* Mean anomaly of the Sun (MHB2000). */
            elp = ((1287104.79305 +
                     t * (129596581.0481 +
                     t * (-0.5532 +
                     t * (0.000136 +
                     t * (-0.00001149))))) % TURNAS) * DAS2R;

            /* Mean longitude of the Moon minus that of the ascending node */
            /* (IERS 2003. */
            f = iauFaf03(t);

            /* Mean elongation of the Moon from the Sun (MHB2000). */
            d = ((1072260.70369 +
                   t * (1602961601.2090 +
                   t * (-6.3706 +
                   t * (0.006593 +
                   t * (-0.00003169))))) % TURNAS) * DAS2R;

            /* Mean longitude of the ascending node of the Moon (IERS 2003). */
            om = iauFaom03(t);

            /* Initialize the nutation values. */
            dp = 0.0;
            de = 0.0;

            /* Summation of luni-solar nutation series (in reverse order). */
            for (i = NLS - 1; i >= 0; i--)
            {

                /* Argument and functions. */
                arg = (((double)Nut00axls[i].nl * el +
                           (double)Nut00axls[i].nlp * elp +
                           (double)Nut00axls[i].nf * f +
                           (double)Nut00axls[i].nd * d +
                           (double)Nut00axls[i].nom * om) % D2PI);
                sarg = Math.Sin(arg);
                carg = Math.Cos(arg);

                /* Term. */
                dp += (Nut00axls[i].sp + Nut00axls[i].spt * t) * sarg + Nut00axls[i].cp * carg;
                de += (Nut00axls[i].ce + Nut00axls[i].cet * t) * carg + Nut00axls[i].se * sarg;
            }

            /* Convert from 0.1 microarcsec units to radians. */
            dpsils = dp * U2R;
            depsls = de * U2R;

            /* ------------------ */
            /* PLANETARY NUTATION */
            /* ------------------ */

            /* n.b.  The MHB2000 code computes the luni-solar and planetary nutation */
            /* in different functions, using slightly different Delaunay */
            /* arguments in the two cases.  This behaviour is faithfully */
            /* reproduced here.  Use of the IERS 2003 expressions for both */
            /* cases leads to negligible changes, well below */
            /* 0.1 microarcsecond. */

            /* Mean anomaly of the Moon (MHB2000). */
            al = ((2.35555598 + 8328.6914269554 * t) % D2PI);

            /* Mean longitude of the Moon minus that of the ascending node */
            /*(MHB2000). */
            af = ((1.627905234 + 8433.466158131 * t) % D2PI);

            /* Mean elongation of the Moon from the Sun (MHB2000). */
            ad = ((5.198466741 + 7771.3771468121 * t) % D2PI);

            /* Mean longitude of the ascending node of the Moon (MHB2000). */
            aom = ((2.18243920 - 33.757045 * t) % D2PI);

            /* General accumulated precession in longitude (IERS 2003). */
            apa = iauFapa03(t);

            /* Planetary longitudes, Mercury through Uranus (IERS 2003). */
            alme = iauFame03(t);
            alve = iauFave03(t);
            alea = iauFae03(t);
            alma = iauFama03(t);
            alju = iauFaju03(t);
            alsa = iauFasa03(t);
            alur = iauFaur03(t);

            /* Neptune longitude (MHB2000). */
            alne = ((5.321159000 + 3.8127774000 * t) % D2PI);

            /* Initialize the nutation values. */
            dp = 0.0;
            de = 0.0;

            /* Summation of planetary nutation series (in reverse order). */
            for (i = NPL - 1; i >= 0; i--)
            {

                /* Argument and functions. */
                arg = (((double)xpl[i].nl * al +
                           (double)xpl[i].nf * af +
                           (double)xpl[i].nd * ad +
                           (double)xpl[i].nom * aom +
                           (double)xpl[i].nme * alme +
                           (double)xpl[i].nve * alve +
                           (double)xpl[i].nea * alea +
                           (double)xpl[i].nma * alma +
                           (double)xpl[i].nju * alju +
                           (double)xpl[i].nsa * alsa +
                           (double)xpl[i].nur * alur +
                           (double)xpl[i].nne * alne +
                           (double)xpl[i].npa * apa) % D2PI);
                sarg = Math.Sin(arg);
                carg = Math.Cos(arg);

                /* Term. */
                dp += (double)xpl[i].sp * sarg + (double)xpl[i].cp * carg;
                de += (double)xpl[i].se * sarg + (double)xpl[i].ce * carg;

            }

            /* Convert from 0.1 microarcsec units to radians. */
            dpsipl = dp * U2R;
            depspl = de * U2R;

            /* ------- */
            /* RESULTS */
            /* ------- */

            /* Add luni-solar and planetary components. */
            dpsi = dpsils + dpsipl;
            deps = depsls + depspl;

        }


        /* --------------------------------------------------- */
        /* Luni-solar nutation: argument and term coefficients */
        /* --------------------------------------------------- */

        /* The units for the sine and cosine coefficients are */
        /* 0.1 microarcsec and the same per Julian century    */

        static readonly LuniSolar[] Nut00bx =
        {
            /* 1-10 */
            new LuniSolar( 0, 0, 0, 0,1,
            -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0),
            new LuniSolar( 0, 0, 2,-2,2,
            -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0),
            new LuniSolar( 0, 0, 2, 0,2,-2276413.0,-234.0, 2796.0, 978459.0,-485.0,1374.0),
            new LuniSolar( 0, 0, 0, 0,2,2074554.0,  207.0, -698.0,-897492.0, 470.0,-291.0),
            new LuniSolar( 0, 1, 0, 0,0,1475877.0,-3633.0,11817.0, 73871.0,-184.0,-1924.0),
            new LuniSolar( 0, 1, 2,-2,2,-516821.0, 1226.0, -524.0, 224386.0,-677.0,-174.0),
            new LuniSolar( 1, 0, 0, 0,0, 711159.0,   73.0, -872.0,  -6750.0,   0.0, 358.0),
            new LuniSolar( 0, 0, 2, 0,1,-387298.0, -367.0,  380.0, 200728.0,  18.0, 318.0),
            new LuniSolar( 1, 0, 2, 0,2,-301461.0,  -36.0,  816.0, 129025.0, -63.0, 367.0),
            new LuniSolar( 0,-1, 2,-2,2, 215829.0, -494.0,  111.0, -95929.0, 299.0, 132.0),

            /* 11-20 */
            new LuniSolar( 0, 0, 2,-2,1, 128227.0,  137.0,  181.0, -68982.0,  -9.0,  39.0),
            new LuniSolar(-1, 0, 2, 0,2, 123457.0,   11.0,   19.0, -53311.0,  32.0,  -4.0),
            new LuniSolar(-1, 0, 0, 2,0, 156994.0,   10.0, -168.0,  -1235.0,   0.0,  82.0),
            new LuniSolar( 1, 0, 0, 0,1,  63110.0,   63.0,   27.0, -33228.0,   0.0,  -9.0),
            new LuniSolar(-1, 0, 0, 0,1, -57976.0,  -63.0, -189.0,  31429.0,   0.0, -75.0),
            new LuniSolar(-1, 0, 2, 2,2, -59641.0,  -11.0,  149.0,  25543.0, -11.0,  66.0),
            new LuniSolar( 1, 0, 2, 0,1, -51613.0,  -42.0,  129.0,  26366.0,   0.0,  78.0),
            new LuniSolar(-2, 0, 2, 0,1,  45893.0,   50.0,   31.0, -24236.0, -10.0,  20.0),
            new LuniSolar( 0, 0, 0, 2,0,  63384.0,   11.0, -150.0,  -1220.0,   0.0,  29.0),
            new LuniSolar( 0, 0, 2, 2,2, -38571.0,   -1.0,  158.0,  16452.0, -11.0,  68.0),

            /* 21-30 */
            new LuniSolar( 0,-2, 2,-2,2,  32481.0,    0.0,    0.0, -13870.0,   0.0,   0.0),
            new LuniSolar(-2, 0, 0, 2,0, -47722.0,    0.0,  -18.0,    477.0,   0.0, -25.0),
            new LuniSolar( 2, 0, 2, 0,2, -31046.0,   -1.0,  131.0,  13238.0, -11.0,  59.0),
            new LuniSolar( 1, 0, 2,-2,2,  28593.0,    0.0,   -1.0, -12338.0,  10.0,  -3.0),
            new LuniSolar(-1, 0, 2, 0,1,  20441.0,   21.0,   10.0, -10758.0,   0.0,  -3.0),
            new LuniSolar( 2, 0, 0, 0,0,  29243.0,    0.0,  -74.0,   -609.0,   0.0,  13.0),
            new LuniSolar( 0, 0, 2, 0,0,  25887.0,    0.0,  -66.0,   -550.0,   0.0,  11.0),
            new LuniSolar( 0, 1, 0, 0,1, -14053.0,  -25.0,   79.0,   8551.0,  -2.0, -45.0),
            new LuniSolar(-1, 0, 0, 2,1,  15164.0,   10.0,   11.0,  -8001.0,   0.0,  -1.0),
            new LuniSolar( 0, 2, 2,-2,2, -15794.0,   72.0,  -16.0,   6850.0, -42.0,  -5.0),

            /* 31-40 */
            new LuniSolar( 0, 0,-2, 2,0,  21783.0,    0.0,   13.0,   -167.0,   0.0,  13.0),
            new LuniSolar( 1, 0, 0,-2,1, -12873.0,  -10.0,  -37.0,   6953.0,   0.0, -14.0),
            new LuniSolar( 0,-1, 0, 0,1, -12654.0,   11.0,   63.0,   6415.0,   0.0,  26.0),
            new LuniSolar(-1, 0, 2, 2,1, -10204.0,    0.0,   25.0,   5222.0,   0.0,  15.0),
            new LuniSolar( 0, 2, 0, 0,0,  16707.0,  -85.0,  -10.0,    168.0,  -1.0,  10.0),
            new LuniSolar( 1, 0, 2, 2,2,  -7691.0,    0.0,   44.0,   3268.0,   0.0,  19.0),
            new LuniSolar(-2, 0, 2, 0,0, -11024.0,    0.0,  -14.0,    104.0,   0.0,   2.0),
            new LuniSolar( 0, 1, 2, 0,2,   7566.0,  -21.0,  -11.0,  -3250.0,   0.0,  -5.0),
            new LuniSolar( 0, 0, 2, 2,1,  -6637.0,  -11.0,   25.0,   3353.0,   0.0,  14.0),
            new LuniSolar( 0,-1, 2, 0,2,  -7141.0,   21.0,    8.0,   3070.0,   0.0,   4.0),

            /* 41-50 */
            new LuniSolar( 0, 0, 0, 2,1,  -6302.0,  -11.0,    2.0,   3272.0,   0.0,   4.0),
            new LuniSolar( 1, 0, 2,-2,1,   5800.0,   10.0,    2.0,  -3045.0,   0.0,  -1.0),
            new LuniSolar( 2, 0, 2,-2,2,   6443.0,    0.0,   -7.0,  -2768.0,   0.0,  -4.0),
            new LuniSolar(-2, 0, 0, 2,1,  -5774.0,  -11.0,  -15.0,   3041.0,   0.0,  -5.0),
            new LuniSolar( 2, 0, 2, 0,1,  -5350.0,    0.0,   21.0,   2695.0,   0.0,  12.0),
            new LuniSolar( 0,-1, 2,-2,1,  -4752.0,  -11.0,   -3.0,   2719.0,   0.0,  -3.0),
            new LuniSolar( 0, 0, 0,-2,1,  -4940.0,  -11.0,  -21.0,   2720.0,   0.0,  -9.0),
            new LuniSolar(-1,-1, 0, 2,0,   7350.0,    0.0,   -8.0,    -51.0,   0.0,   4.0),
            new LuniSolar( 2, 0, 0,-2,1,   4065.0,    0.0,    6.0,  -2206.0,   0.0,   1.0),
            new LuniSolar( 1, 0, 0, 2,0,   6579.0,    0.0,  -24.0,   -199.0,   0.0,   2.0),

            /* 51-60 */
            new LuniSolar( 0, 1, 2,-2,1,   3579.0,    0.0,    5.0,  -1900.0,   0.0,   1.0),
            new LuniSolar( 1,-1, 0, 0,0,   4725.0,    0.0,   -6.0,    -41.0,   0.0,   3.0),
            new LuniSolar(-2, 0, 2, 0,2,  -3075.0,    0.0,   -2.0,   1313.0,   0.0,  -1.0),
            new LuniSolar( 3, 0, 2, 0,2,  -2904.0,    0.0,   15.0,   1233.0,   0.0,   7.0),
            new LuniSolar( 0,-1, 0, 2,0,   4348.0,    0.0,  -10.0,    -81.0,   0.0,   2.0),
            new LuniSolar( 1,-1, 2, 0,2,  -2878.0,    0.0,    8.0,   1232.0,   0.0,   4.0),
            new LuniSolar( 0, 0, 0, 1,0,  -4230.0,    0.0,    5.0,    -20.0,   0.0,  -2.0),
            new LuniSolar(-1,-1, 2, 2,2,  -2819.0,    0.0,    7.0,   1207.0,   0.0,   3.0),
            new LuniSolar(-1, 0, 2, 0,0,  -4056.0,    0.0,    5.0,     40.0,   0.0,  -2.0),
            new LuniSolar( 0,-1, 2, 2,2,  -2647.0,    0.0,   11.0,   1129.0,   0.0,   5.0),

            /* 61-70 */
            new LuniSolar(-2, 0, 0, 0,1,  -2294.0,    0.0,  -10.0,   1266.0,   0.0,  -4.0),
            new LuniSolar( 1, 1, 2, 0,2,   2481.0,    0.0,   -7.0,  -1062.0,   0.0,  -3.0),
            new LuniSolar( 2, 0, 0, 0,1,   2179.0,    0.0,   -2.0,  -1129.0,   0.0,  -2.0),
            new LuniSolar(-1, 1, 0, 1,0,   3276.0,    0.0,    1.0,     -9.0,   0.0,   0.0),
            new LuniSolar( 1, 1, 0, 0,0,  -3389.0,    0.0,    5.0,     35.0,   0.0,  -2.0),
            new LuniSolar( 1, 0, 2, 0,0,   3339.0,    0.0,  -13.0,   -107.0,   0.0,   1.0),
            new LuniSolar(-1, 0, 2,-2,1,  -1987.0,    0.0,   -6.0,   1073.0,   0.0,  -2.0),
            new LuniSolar( 1, 0, 0, 0,2,  -1981.0,    0.0,    0.0,    854.0,   0.0,   0.0),
            new LuniSolar(-1, 0, 0, 1,0,   4026.0,    0.0, -353.0,   -553.0,   0.0,-139.0),
            new LuniSolar( 0, 0, 2, 1,2,   1660.0,    0.0,   -5.0,   -710.0,   0.0,  -2.0),

            /* 71-77 */
            new LuniSolar(-1, 0, 2, 4,2,  -1521.0,    0.0,    9.0,    647.0,   0.0,   4.0),
            new LuniSolar(-1, 1, 0, 1,1,   1314.0,    0.0,    0.0,   -700.0,   0.0,   0.0),
            new LuniSolar( 0,-2, 2,-2,1,  -1283.0,    0.0,    0.0,    672.0,   0.0,   0.0),
            new LuniSolar( 1, 0, 2, 2,1,  -1331.0,    0.0,    8.0,    663.0,   0.0,   4.0),
            new LuniSolar(-2, 0, 2, 2,2,   1383.0,    0.0,   -2.0,   -594.0,   0.0,  -2.0),
            new LuniSolar(-1, 0, 0, 0,2,   1405.0,    0.0,    4.0,   -610.0,   0.0,   2.0),
            new LuniSolar( 1, 1, 2,-2,2,   1290.0,    0.0,    0.0,   -556.0,   0.0,   0.0)
        };

        /// <summary>
        /// Nutation, IAU 2000B model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation, luni-solar + planetary longitude</param>
        /// <param name="deps">nutation, luni-solar + planetary obliquity</param>
        public static void iauNut00b(double date1, double date2, out double dpsi, out double deps)
        {
            double t, el, elp, f, d, om, arg, dp, de, sarg, carg,
                   dpsils, depsls, dpsipl, depspl;
            int i;

            /* Units of 0.1 microarcsecond to radians */
            double U2R = DAS2R / 1e7;

            /* ---------------------------------------- */
            /* Fixed offsets in lieu of planetary terms */
            /* ---------------------------------------- */

            double DPPLAN = -0.135 * DMAS2R;
            double DEPLAN = 0.388 * DMAS2R;

            /* Number of terms in the series */
            //const int NLS = (int)(sizeof x / sizeof x[0]);
            int NLS = Nut00bx.Length;


            /* ------------------------------------------------------------------ */

            /* Interval between fundamental epoch J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* --------------------*/
            /* LUNI-SOLAR NUTATION */
            /* --------------------*/

            /* Fundamental (Delaunay) arguments from Simon et al. (1994) */

            /* Mean anomaly of the Moon. */
            el = ((485868.249036 + (1717915923.2178) * t) % TURNAS) * DAS2R;

            /* Mean anomaly of the Sun. */
            elp = ((1287104.79305 + (129596581.0481) * t) % TURNAS) * DAS2R;

            /* Mean argument of the latitude of the Moon. */
            f = ((335779.526232 + (1739527262.8478) * t) % TURNAS) * DAS2R;

            /* Mean elongation of the Moon from the Sun. */
            d = ((1072260.70369 + (1602961601.2090) * t) % TURNAS) * DAS2R;

            /* Mean longitude of the ascending node of the Moon. */
            om = ((450160.398036 + (-6962890.5431) * t) % TURNAS) * DAS2R;

            /* Initialize the nutation values. */
            dp = 0.0;
            de = 0.0;

            /* Summation of luni-solar nutation series (smallest terms first). */
            for (i = NLS - 1; i >= 0; i--)
            {

                /* Argument and functions. */
                arg = (((double)Nut00bx[i].nl * el +
                            (double)Nut00bx[i].nlp * elp +
                            (double)Nut00bx[i].nf * f +
                            (double)Nut00bx[i].nd * d +
                            (double)Nut00bx[i].nom * om) % D2PI);
                sarg = Math.Sin(arg);
                carg = Math.Cos(arg);

                /* Term. */
                dp += (Nut00bx[i].sp + Nut00bx[i].spt * t) * sarg + Nut00bx[i].cp * carg;
                de += (Nut00bx[i].ce + Nut00bx[i].cet * t) * carg + Nut00bx[i].se * sarg;
            }

            /* Convert from 0.1 microarcsec units to radians. */
            dpsils = dp * U2R;
            depsls = de * U2R;

            /* ------------------------------*/
            /* IN LIEU OF PLANETARY NUTATION */
            /* ------------------------------*/

            /* Fixed offset to correct for missing terms in truncated series. */
            dpsipl = DPPLAN;
            depspl = DEPLAN;

            /* --------*/
            /* RESULTS */
            /* --------*/

            /* Add luni-solar and planetary components. */
            dpsi = dpsils + dpsipl;
            deps = depsls + depspl;
        }


        /// <summary>
        /// IAU 2000A nutation with adjustments to match the IAU 2006
        /// precession.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation, luni-solar + planetary longitude (radians)</param>
        /// <param name="deps">nutation, luni-solar + planetary obliquity (radians)</param>
        public static void iauNut06a(double date1, double date2, out double dpsi, out double deps)
        {
            double t, fj2, dp, de;


            /* Interval between fundamental date J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Factor correcting for secular variation of J2. */
            fj2 = -2.7774e-6 * t;

            /* Obtain IAU 2000A nutation. */
            iauNut00a(date1, date2, out dp, out de);

            /* Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5). */
            dpsi = dp + dp * (0.4697e-6 + fj2);
            deps = de + de * fj2;
        }


        readonly struct Nut80Coef
        {
            public readonly int nl, nlp, nf, nd, nom; /* coefficients of l,l',F,D,Om */
            public readonly double sp, spt;        /* longitude sine, 1 and t coefficients */
            public readonly double ce, cet;        /* obliquity cosine, 1 and t coefficients */

            public Nut80Coef(int nl, int nlp, int nf, int nd, int nom, double sp, double spt, double ce, double cet)
            {
                this.nl = nl;
                this.nlp = nlp;
                this.nf = nf;
                this.nd = nd;
                this.nom = nom;
                this.sp = sp;
                this.spt = spt;
                this.ce = ce;
                this.cet = cet;
            }
        }

        /* ------------------------------------------------ */
        /* Table of multiples of arguments and coefficients */
        /* ------------------------------------------------ */

        /* The units for the sine and cosine coefficients are 0.1 mas and */
        /* the same per Julian century */
        static readonly Nut80Coef[] x =
        {
            /* 1-10 */
            new Nut80Coef(  0,  0,  0,  0,  1, -171996.0, -174.2,  92025.0,    8.9 ),
            new Nut80Coef(  0,  0,  0,  0,  2,    2062.0,    0.2,   -895.0,    0.5 ),
            new Nut80Coef( -2,  0,  2,  0,  1,      46.0,    0.0,    -24.0,    0.0 ),
            new Nut80Coef(  2,  0, -2,  0,  0,      11.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef( -2,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  1, -1,  0, -1,  0,      -3.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0, -2,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  2,  0, -2,  0,  1,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0,  2, -2,  2,  -13187.0,   -1.6,   5736.0,   -3.1 ),
            new Nut80Coef(  0,  1,  0,  0,  0,    1426.0,   -3.4,     54.0,   -0.1 ),

            /* 11-20 */
            new Nut80Coef(  0,  1,  2, -2,  2,    -517.0,    1.2,    224.0,   -0.6 ),
            new Nut80Coef(  0, -1,  2, -2,  2,     217.0,   -0.5,    -95.0,    0.3 ),
            new Nut80Coef(  0,  0,  2, -2,  1,     129.0,    0.1,    -70.0,    0.0 ),
            new Nut80Coef(  2,  0,  0, -2,  0,      48.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  0,  0,  2, -2,  0,     -22.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  2,  0,  0,  0,      17.0,   -0.1,      0.0,    0.0 ),
            new Nut80Coef(  0,  1,  0,  0,  1,     -15.0,    0.0,      9.0,    0.0 ),
            new Nut80Coef(  0,  2,  2, -2,  2,     -16.0,    0.1,      7.0,    0.0 ),
            new Nut80Coef(  0, -1,  0,  0,  1,     -12.0,    0.0,      6.0,    0.0 ),
            new Nut80Coef( -2,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 ),

            /* 21-30 */
            new Nut80Coef(  0, -1,  2, -2,  1,      -5.0,    0.0,      3.0,    0.0 ),
            new Nut80Coef(  2,  0,  0, -2,  1,       4.0,    0.0,     -2.0,    0.0 ),
            new Nut80Coef(  0,  1,  2, -2,  1,       4.0,    0.0,     -2.0,    0.0 ),
            new Nut80Coef(  1,  0,  0, -1,  0,      -4.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  2,  1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0, -2,  2,  1,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  1, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  1,  0,  0,  2,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef( -1,  0,  0,  1,  1,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  1,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ),

            /* 31-40 */
            new Nut80Coef(  0,  0,  2,  0,  2,   -2274.0,   -0.2,    977.0,   -0.5 ),
            new Nut80Coef(  1,  0,  0,  0,  0,     712.0,    0.1,     -7.0,    0.0 ),
            new Nut80Coef(  0,  0,  2,  0,  1,    -386.0,   -0.4,    200.0,    0.0 ),
            new Nut80Coef(  1,  0,  2,  0,  2,    -301.0,    0.0,    129.0,   -0.1 ),
            new Nut80Coef(  1,  0,  0, -2,  0,    -158.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef( -1,  0,  2,  0,  2,     123.0,    0.0,    -53.0,    0.0 ),
            new Nut80Coef(  0,  0,  0,  2,  0,      63.0,    0.0,     -2.0,    0.0 ),
            new Nut80Coef(  1,  0,  0,  0,  1,      63.0,    0.1,    -33.0,    0.0 ),
            new Nut80Coef( -1,  0,  0,  0,  1,     -58.0,   -0.1,     32.0,    0.0 ),
            new Nut80Coef( -1,  0,  2,  2,  2,     -59.0,    0.0,     26.0,    0.0 ),

            /* 41-50 */
            new Nut80Coef(  1,  0,  2,  0,  1,     -51.0,    0.0,     27.0,    0.0 ),
            new Nut80Coef(  0,  0,  2,  2,  2,     -38.0,    0.0,     16.0,    0.0 ),
            new Nut80Coef(  2,  0,  0,  0,  0,      29.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef(  1,  0,  2, -2,  2,      29.0,    0.0,    -12.0,    0.0 ),
            new Nut80Coef(  2,  0,  2,  0,  2,     -31.0,    0.0,     13.0,    0.0 ),
            new Nut80Coef(  0,  0,  2,  0,  0,      26.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef( -1,  0,  2,  0,  1,      21.0,    0.0,    -10.0,    0.0 ),
            new Nut80Coef( -1,  0,  0,  2,  1,      16.0,    0.0,     -8.0,    0.0 ),
            new Nut80Coef(  1,  0,  0, -2,  1,     -13.0,    0.0,      7.0,    0.0 ),
            new Nut80Coef( -1,  0,  2,  2,  1,     -10.0,    0.0,      5.0,    0.0 ),

            /* 51-60 */
            new Nut80Coef(  1,  1,  0, -2,  0,      -7.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  1,  2,  0,  2,       7.0,    0.0,     -3.0,    0.0 ),
            new Nut80Coef(  0, -1,  2,  0,  2,      -7.0,    0.0,      3.0,    0.0 ),
            new Nut80Coef(  1,  0,  2,  2,  2,      -8.0,    0.0,      3.0,    0.0 ),
            new Nut80Coef(  1,  0,  0,  2,  0,       6.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  2,  0,  2, -2,  2,       6.0,    0.0,     -3.0,    0.0 ),
            new Nut80Coef(  0,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 ),
            new Nut80Coef(  0,  0,  2,  2,  1,      -7.0,    0.0,      3.0,    0.0 ),
            new Nut80Coef(  1,  0,  2, -2,  1,       6.0,    0.0,     -3.0,    0.0 ),
            new Nut80Coef(  0,  0,  0, -2,  1,      -5.0,    0.0,      3.0,    0.0 ),

            /* 61-70 */
            new Nut80Coef(  1, -1,  0,  0,  0,       5.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  2,  0,  2,  0,  1,      -5.0,    0.0,      3.0,    0.0 ),
            new Nut80Coef(  0,  1,  0, -2,  0,      -4.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  0, -2,  0,  0,       4.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0,  0,  1,  0,      -4.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  1,  0,  0,  0,      -3.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  0,  2,  0,  0,       3.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1, -1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef( -1, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef( -2,  0,  0,  0,  1,      -2.0,    0.0,      1.0,    0.0 ),

            /* 71-80 */
            new Nut80Coef(  3,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  0, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  1,  1,  2,  0,  2,       2.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef( -1,  0,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  2,  0,  0,  0,  1,       2.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef(  1,  0,  0,  0,  2,      -2.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  3,  0,  0,  0,  0,       2.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0,  2,  1,  2,       2.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef( -1,  0,  0,  0,  2,       1.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef(  1,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 ),

            /* 81-90 */
            new Nut80Coef( -2,  0,  2,  2,  2,       1.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef( -1,  0,  2,  4,  2,      -2.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef(  2,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  1,  2, -2,  2,       1.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef(  1,  0,  2,  2,  1,      -1.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef( -2,  0,  2,  4,  2,      -1.0,    0.0,      1.0,    0.0 ),
            new Nut80Coef( -1,  0,  4,  0,  2,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1, -1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  2,  0,  2, -2,  1,       1.0,    0.0,     -1.0,    0.0 ),
            new Nut80Coef(  2,  0,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0 ),

            /* 91-100 */
            new Nut80Coef(  1,  0,  0,  2,  1,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0,  4, -2,  2,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  3,  0,  2, -2,  2,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  0,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  1,  2,  0,  1,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef( -1, -1,  0,  2,  1,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0, -2,  0,  1,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0,  2, -1,  2,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  1,  0,  2,  0,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  0, -2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ),

            /* 101-106 */
            new Nut80Coef(  0, -1,  2,  0,  1,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  1,  0, -2,  1,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  1,  0, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  2,  0,  0,  2,  0,       1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  0,  2,  4,  2,      -1.0,    0.0,      0.0,    0.0 ),
            new Nut80Coef(  0,  1,  0,  1,  0,       1.0,    0.0,      0.0,    0.0 )
        };

        /// <summary>
        /// Nutation, IAU 1980 model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation in longitude (radians)</param>
        /// <param name="deps">nutation in obliquity (radians)</param>
        public static void iauNut80(double date1, double date2, out double dpsi, out double deps)
        {
            double t, el, elp, f, d, om, dp, de, arg, s, c;
            int j;

            /* Units of 0.1 milliarcsecond to radians */
            double U2R = DAS2R / 1e4;

            /* Number of terms in the series */
            //const int NT = (int)(sizeof x / sizeof x[0]);
            int NT = x.Length;

            /* ------------------------------------------------------------------ */

            /* Interval between fundamental epoch J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* --------------------- */
            /* Fundamental arguments */
            /* --------------------- */

            /* Mean longitude of Moon minus mean longitude of Moon's perigee. */
            el = iauAnpm(
                 (485866.733 + (715922.633 + (31.310 + 0.064 * t) * t) * t)
                 * DAS2R + ((1325.0 * t) % 1.0) * D2PI);

            /* Mean longitude of Sun minus mean longitude of Sun's perigee. */
            elp = iauAnpm(
                  (1287099.804 + (1292581.224 + (-0.577 - 0.012 * t) * t) * t)
                  * DAS2R + ((99.0 * t) % 1.0) * D2PI);

            /* Mean longitude of Moon minus mean longitude of Moon's node. */
            f = iauAnpm(
                (335778.877 + (295263.137 + (-13.257 + 0.011 * t) * t) * t)
                * DAS2R + ((1342.0 * t) % 1.0) * D2PI);

            /* Mean elongation of Moon from Sun. */
            d = iauAnpm(
                (1072261.307 + (1105601.328 + (-6.891 + 0.019 * t) * t) * t)
                * DAS2R + ((1236.0 * t) % 1.0) * D2PI);

            /* Longitude of the mean ascending node of the lunar orbit on the */
            /* ecliptic, measured from the mean equinox of date. */
            om = iauAnpm(
                 (450160.280 + (-482890.539 + (7.455 + 0.008 * t) * t) * t)
                 * DAS2R + ((-5.0 * t) % 1.0) * D2PI);

            /* --------------- */
            /* Nutation series */
            /* --------------- */

            /* Initialize nutation components. */
            dp = 0.0;
            de = 0.0;

            /* Sum the nutation terms, ending with the biggest. */
            for (j = NT - 1; j >= 0; j--)
            {

                /* Form argument for current term. */
                arg = (double)x[j].nl * el
                    + (double)x[j].nlp * elp
                    + (double)x[j].nf * f
                    + (double)x[j].nd * d
                    + (double)x[j].nom * om;

                /* Accumulate current nutation term. */
                s = x[j].sp + x[j].spt * t;
                c = x[j].ce + x[j].cet * t;
                if (s != 0.0) dp += s * Math.Sin(arg);
                if (c != 0.0) de += c * Math.Cos(arg);
            }

            /* Convert results from 0.1 mas units to radians. */
            dpsi = dp * U2R;
            deps = de * U2R;
        }


        /// <summary>
        /// Form the matrix of nutation for a given date, IAU 1980 model.
        /// </summary>
        /// <param name="date1">TT date</param>
        /// <param name="date2">TT date</param>
        /// <param name="rmatn">double[3][3]    nutation matrix</param>
        public static void iauNutm80(double date1, double date2, out double[][] rmatn)
        {
            double dpsi, deps, epsa;


            /* Nutation components and mean obliquity. */
            iauNut80(date1, date2, out dpsi, out deps);
            epsa = iauObl80(date1, date2);

            /* Build the rotation matrix. */
            iauNumat(epsa, dpsi, deps, out rmatn);
        }


        /// <summary>
        /// Mean obliquity of the ecliptic, IAU 2006 precession model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>obliquity of the ecliptic (radians)</returns>
        public static double iauObl06(double date1, double date2)
        {
            double t, eps0;


            /* Interval between fundamental date J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Mean obliquity. */
            eps0 = (84381.406 +
                   (-46.836769 +
                   (-0.0001831 +
                   (0.00200340 +
                   (-0.000000576 +
                   (-0.0000000434) * t) * t) * t) * t) * t) * DAS2R;

            return eps0;
        }


        /// <summary>
        /// Mean obliquity of the ecliptic, IAU 1980 model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>obliquity of the ecliptic (radians)</returns>
        public static double iauObl80(double date1, double date2)
        {
            double t, eps0;


            /* Interval between fundamental epoch J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Mean obliquity of date. */
            eps0 = DAS2R * (84381.448 +
                           (-46.8150 +
                           (-0.00059 +
                           (0.001813) * t) * t) * t);

            return eps0;
        }


        /// <summary>
        /// Precession angles, IAU 2006, equinox based.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="eps0">epsilon_0</param>
        /// <param name="psia">psi_A</param>
        /// <param name="oma">omega_A</param>
        /// <param name="bpa">P_A</param>
        /// <param name="bqa">Q_A</param>
        /// <param name="pia">pi_A</param>
        /// <param name="bpia">pi_A</param>
        /// <param name="epsa">obliquity epsilon_A</param>
        /// <param name="chia">chi_A</param>
        /// <param name="za">z_A</param>
        /// <param name="zetaa">zeta_A</param>
        /// <param name="thetaa">theta_A</param>
        /// <param name="pa">P_A</param>
        /// <param name="gam">F-W angle gamma_J2000</param>
        /// <param name="phi">F-W angle phi_J2000</param>
        /// <param name="psi">F-W angle psi_J2000</param>
        public static void iauP06e(double date1, double date2,
             out double eps0, out double psia, out double oma, out double bpa,
             out double bqa, out double pia, out double bpia,
             out double epsa, out double chia, out double za, out double zetaa,
             out double thetaa, out double pa,
             out double gam, out double phi, out double psi)
        {
            double t;


            /* Interval between fundamental date J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Obliquity at J2000.0. */

            eps0 = 84381.406 * DAS2R;

            /* Luni-solar precession. */

            psia = (5038.481507 +
                    (-1.0790069 +
                    (-0.00114045 +
                    (0.000132851 +
                    (-0.0000000951)
                    * t) * t) * t) * t) * t * DAS2R;

            /* Inclination of mean equator with respect to the J2000.0 ecliptic. */

            oma = eps0 + (-0.025754 +
                           (0.0512623 +
                           (-0.00772503 +
                           (-0.000000467 +
                           (0.0000003337)
                           * t) * t) * t) * t) * t * DAS2R;

            /* Ecliptic pole x, J2000.0 ecliptic triad. */

            bpa = (4.199094 +
                   (0.1939873 +
                   (-0.00022466 +
                   (-0.000000912 +
                   (0.0000000120)
                   * t) * t) * t) * t) * t * DAS2R;

            /* Ecliptic pole -y, J2000.0 ecliptic triad. */

            bqa = (-46.811015 +
                   (0.0510283 +
                   (0.00052413 +
                   (-0.000000646 +
                   (-0.0000000172)
                   * t) * t) * t) * t) * t * DAS2R;

            /* Angle between moving and J2000.0 ecliptics. */

            pia = (46.998973 +
                   (-0.0334926 +
                   (-0.00012559 +
                   (0.000000113 +
                   (-0.0000000022)
                   * t) * t) * t) * t) * t * DAS2R;

            /* Longitude of ascending node of the moving ecliptic. */

            bpia = (629546.7936 +
                    (-867.95758 +
                    (0.157992 +
                    (-0.0005371 +
                    (-0.00004797 +
                    (0.000000072)
                    * t) * t) * t) * t) * t) * DAS2R;

            /* Mean obliquity of the ecliptic. */

            epsa = iauObl06(date1, date2);

            /* Planetary precession. */

            chia = (10.556403 +
                    (-2.3814292 +
                    (-0.00121197 +
                    (0.000170663 +
                    (-0.0000000560)
                    * t) * t) * t) * t) * t * DAS2R;

            /* Equatorial precession: minus the third of the 323 Euler angles. */

            za = (-2.650545 +
                  (2306.077181 +
                  (1.0927348 +
                  (0.01826837 +
                  (-0.000028596 +
                  (-0.0000002904)
                  * t) * t) * t) * t) * t) * DAS2R;

            /* Equatorial precession: minus the first of the 323 Euler angles. */

            zetaa = (2.650545 +
                     (2306.083227 +
                     (0.2988499 +
                     (0.01801828 +
                     (-0.000005971 +
                     (-0.0000003173)
                     * t) * t) * t) * t) * t) * DAS2R;

            /* Equatorial precession: second of the 323 Euler angles. */

            thetaa = (2004.191903 +
                      (-0.4294934 +
                      (-0.04182264 +
                      (-0.000007089 +
                      (-0.0000001274)
                      * t) * t) * t) * t) * t * DAS2R;

            /* General precession. */

            pa = (5028.796195 +
                  (1.1054348 +
                  (0.00007964 +
                  (-0.000023857 +
                  (-0.0000000383)
                  * t) * t) * t) * t) * t * DAS2R;

            /* Fukushima-Williams angles for precession. */

            gam = (10.556403 +
                   (0.4932044 +
                   (-0.00031238 +
                   (-0.000002788 +
                   (0.0000000260)
                   * t) * t) * t) * t) * t * DAS2R;

            phi = eps0 + (-46.811015 +
                           (0.0511269 +
                           (0.00053289 +
                           (-0.000000440 +
                           (-0.0000000176)
                           * t) * t) * t) * t) * t * DAS2R;

            psi = (5038.481507 +
                   (1.5584176 +
                   (-0.00018522 +
                   (-0.000026452 +
                   (-0.0000000148)
                   * t) * t) * t) * t) * t * DAS2R;
        }


        /// <summary>
        /// This function forms three Euler angles which implement general
        /// precession from epoch J2000.0, using the IAU 2006 model.  Frame
        /// bias (the offset between ICRS and mean J2000.0) is included.
        /// <para>
        /// The precession-bias matrix is
        /// R_3(-z) x R_2(+theta) x R_3(-zeta).
        /// </para>
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="bzeta">1st rotation: radians cw around z</param>
        /// <param name="bz">3rd rotation: radians cw around z</param>
        /// <param name="btheta">2nd rotation: radians ccw around y</param>
        public static void iauPb06(double date1, double date2,
             out double bzeta, out double bz, out double btheta)
        {
            double y, x;
            double[][] r;

            /* Precession matrix via Fukushima-Williams angles. */
            iauPmat06(date1, date2, out r);

            /* Solve for z, choosing the +/- pi alternative. */
            y = r[1][2];
            x = -r[0][2];
            if (x < 0.0)
            {
                y = -y;
                x = -x;
            }
            bz = (x != 0.0 || y != 0.0) ? -Math.Atan2(y, x) : 0.0;

            /* Derotate it out of the matrix. */
            iauRz(bz, r);

            /* Solve for the remaining two angles. */
            y = r[0][2];
            x = r[2][2];
            btheta = (x != 0.0 || y != 0.0) ? -Math.Atan2(y, x) : 0.0;

            y = -r[1][0];
            x = r[1][1];
            bzeta = (x != 0.0 || y != 0.0) ? -Math.Atan2(y, x) : 0.0;
        }


        /// <summary>
        /// Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="gamb">F-W angle gamma_bar (radians)</param>
        /// <param name="phib">F-W angle phi_bar (radians)</param>
        /// <param name="psib">F-W angle psi_bar (radians)</param>
        /// <param name="epsa">F-W angle epsilon_A (radians)</param>
        public static void iauPfw06(double date1, double date2,
              out double gamb, out double phib, out double psib, out double epsa)
        {
            double t;


            /* Interval between fundamental date J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* P03 bias+precession angles. */
            gamb = (-0.052928 +
                    (10.556378 +
                    (0.4932044 +
                    (-0.00031238 +
                    (-0.000002788 +
                    (0.0000000260)
                    * t) * t) * t) * t) * t) * DAS2R;
            phib = (84381.412819 +
                    (-46.811016 +
                    (0.0511268 +
                    (0.00053289 +
                    (-0.000000440 +
                    (-0.0000000176)
                    * t) * t) * t) * t) * t) * DAS2R;
            psib = (-0.041775 +
                    (5038.481484 +
                    (1.5584175 +
                    (-0.00018522 +
                    (-0.000026452 +
                    (-0.0000000148)
                    * t) * t) * t) * t) * t) * DAS2R;
            epsa = iauObl06(date1, date2);
        }


        /// <summary>
        /// Precession matrix (including frame bias) from GCRS to a specified
        /// date, IAU 2000 model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rbp">double[3][3]    bias-precession matrix</param>
        public static void iauPmat00(double date1, double date2, out double[][] rbp)
        {
            double[][] rb;
            double[][] rp;

            /* Obtain the required matrix (discarding others). */
            iauBp00(date1, date2, out rb, out rp, out rbp);
        }


        /// <summary>
        /// Precession matrix (including frame bias) from GCRS to a specified
        /// date, IAU 2006 model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rbp">double[3][3]    bias-precession matrix</param>
        public static void iauPmat06(double date1, double date2, out double[][] rbp)
        {
            double gamb, phib, psib, epsa;

            /* Bias-precession Fukushima-Williams angles. */
            iauPfw06(date1, date2, out gamb, out phib, out psib, out epsa);

            /* Form the matrix. */
            iauFw2m(gamb, phib, psib, epsa, out rbp);
        }


        /// <summary>
        /// Precession matrix from J2000.0 to a specified date, IAU 1976 model.
        /// </summary>
        /// <param name="date1">ending date, TT</param>
        /// <param name="date2">ending date, TT</param>
        /// <param name="rmatp">double[3][3] precession matrix, J2000.0 -> date1+date2</param>
        public static void iauPmat76(double date1, double date2, out double[][] rmatp)
        {
            double zeta, z, theta;
            double[][] wmat;
            rmatp = new double[][] { new double[3], new double[3], new double[3] };

            /* Precession Euler angles, J2000.0 to specified date. */
            iauPrec76(DJ00, 0.0, date1, date2, out zeta, out z, out theta);

            /* Form the rotation matrix. */
            iauIr(out wmat);
            iauRz(-zeta, wmat);
            iauRy(theta, wmat);
            iauRz(-z, wmat);
            iauCr(wmat, rmatp);
        }


        /// <summary>
        /// Precession-nutation, IAU 2000 model:  a multi-purpose function,
        /// supporting classical (equinox-based) use directly and CIO-based
        /// use indirectly.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation longitude</param>
        /// <param name="deps">nutation obliquity</param>
        /// <param name="epsa">mean obliquity</param>
        /// <param name="rb">double[3][3]    frame bias matrix</param>
        /// <param name="rp">double[3][3]    precession matrix</param>
        /// <param name="rbp">double[3][3]    bias-precession matrix</param>
        /// <param name="rn">double[3][3]    nutation matrix</param>
        /// <param name="rbpn">double[3][3]    GCRS-to-true matrix</param>
        public static void iauPn00(double date1, double date2, double dpsi, double deps,
             out double epsa,
             out double[][] rb, out double[][] rp, out double[][] rbp,
             out double[][] rn, out double[][] rbpn)
        {
            double dpsipr, depspr;
            double[][] rbpw;
            double[][] rnw;

            rbp = new double[][] { new double[3], new double[3], new double[3] };
            rn = new double[][] { new double[3], new double[3], new double[3] };
            rbpn = new double[][] { new double[3], new double[3], new double[3] };

            /* IAU 2000 precession-rate adjustments. */
            iauPr00(date1, date2, out dpsipr, out depspr);

            /* Mean obliquity, consistent with IAU 2000 precession-nutation. */
            epsa = iauObl80(date1, date2) + depspr;

            /* Frame bias and precession matrices and their product. */
            iauBp00(date1, date2, out rb, out rp, out rbpw);
            iauCr(rbpw, rbp);

            /* Nutation matrix. */
            iauNumat(epsa, dpsi, deps, out rnw);
            iauCr(rnw, rn);

            /* Bias-precession-nutation matrix (classical). */
            iauRxr(rnw, rbpw, ref rbpn);
        }


        /// <summary>
        /// Precession-nutation, IAU 2000A model:  a multi-purpose function,
        /// supporting classical (equinox-based) use directly and CIO-based
        /// use indirectly.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation longitude</param>
        /// <param name="deps">nutation obliquity</param>
        /// <param name="epsa">mean obliquity</param>
        /// <param name="rb">double[3][3]    frame bias matrix</param>
        /// <param name="rp">double[3][3]    precession matrix</param>
        /// <param name="rbp">double[3][3]    bias-precession matrix</param>
        /// <param name="rn">double[3][3]    nutation matrix</param>
        /// <param name="rbpn">double[3][3]    GCRS-to-true matrix</param>
        public static void iauPn00a(double date1, double date2,
              out double dpsi, out double deps, out double epsa,
              out double[][] rb, out double[][] rp, out double[][] rbp,
              out double[][] rn, out double[][] rbpn)
        {
            /* Nutation. */
            iauNut00a(date1, date2, out dpsi, out deps);

            /* Remaining results. */
            iauPn00(date1, date2, dpsi, deps, out epsa, out rb, out rp, out rbp, out rn, out rbpn);
        }


        /// <summary>
        /// Precession-nutation, IAU 2000B model:  a multi-purpose function,
        /// supporting classical (equinox-based) use directly and CIO-based
        /// use indirectly.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation longitude</param>
        /// <param name="deps">nutation obliquity</param>
        /// <param name="epsa">mean obliquity</param>
        /// <param name="rb">double[3][3]    frame bias matrix</param>
        /// <param name="rp">double[3][3]    precession matrix</param>
        /// <param name="rbp">double[3][3]    bias-precession matrix</param>
        /// <param name="rn">double[3][3]    nutation matrix</param>
        /// <param name="rbpn">double[3][3]    GCRS-to-true matrix</param>
        public static void iauPn00b(double date1, double date2,
              out double dpsi, out double deps, out double epsa,
              out double[][] rb, out double[][] rp, out double[][] rbp,
              out double[][] rn, out double[][] rbpn)
        {
            /* Nutation. */
            iauNut00b(date1, date2, out dpsi, out deps);

            /* Remaining results. */
            iauPn00(date1, date2, dpsi, deps, out epsa, out rb, out rp, out rbp, out rn, out rbpn);
        }


        /// <summary>
        /// Precession-nutation, IAU 2006 model:  a multi-purpose function,
        /// supporting classical (equinox-based) use directly and CIO-based use
        /// indirectly.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation longitude</param>
        /// <param name="deps">nutation obliquity</param>
        /// <param name="epsa">mean obliquity</param>
        /// <param name="rb">double[3][3]    frame bias matrix</param>
        /// <param name="rp">double[3][3]    precession matrix</param>
        /// <param name="rbp">double[3][3]    bias-precession matrix</param>
        /// <param name="rn">double[3][3]    nutation matrix</param>
        /// <param name="rbpn">double[3][3]    GCRS-to-true matrix</param>
        public static void iauPn06(double date1, double date2, double dpsi, double deps,
             out double epsa,
             out double[][] rb, out double[][] rp, out double[][] rbp,
             out double[][] rn, out double[][] rbpn)
        {
            double gamb, phib, psib, eps;
            double[][] r1;
            double[][] r2;
            double[][] rt = new double[][] { new double[3], new double[3], new double[3] };

            rb = new double[][] { new double[3], new double[3], new double[3] };
            rp = new double[][] { new double[3], new double[3], new double[3] };
            rn = new double[][] { new double[3], new double[3], new double[3] };
            rbp = new double[][] { new double[3], new double[3], new double[3] };
            rbpn = new double[][] { new double[3], new double[3], new double[3] };

            /* Bias-precession Fukushima-Williams angles of J2000.0 = frame bias. */
            iauPfw06(DJM0, DJM00, out gamb, out phib, out psib, out eps);

            /* B matrix. */
            iauFw2m(gamb, phib, psib, eps, out r1);
            iauCr(r1, rb);

            /* Bias-precession Fukushima-Williams angles of date. */
            iauPfw06(date1, date2, out gamb, out phib, out psib, out eps);

            /* Bias-precession matrix. */
            iauFw2m(gamb, phib, psib, eps, out r2);
            iauCr(r2, rbp);

            /* Solve for precession matrix. */
            iauTr(r1, ref rt);
            iauRxr(r2, rt, ref rp);

            /* Equinox-based bias-precession-nutation matrix. */
            iauFw2m(gamb, phib, psib + dpsi, eps + deps, out r1);
            iauCr(r1, rbpn);

            /* Solve for nutation matrix. */
            iauTr(r2, ref rt);
            iauRxr(r1, rt, ref rn);

            /* Obliquity, mean of date. */
            epsa = eps;
        }


        /// <summary>
        /// Precession-nutation, IAU 2006/2000A models:  a multi-purpose function,
        /// supporting classical (equinox-based) use directly and CIO-based use
        /// indirectly.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsi">nutation longitude</param>
        /// <param name="deps">nutation obliquity</param>
        /// <param name="epsa">mean obliquity</param>
        /// <param name="rb">double[3][3]    frame bias matrix</param>
        /// <param name="rp">double[3][3]    precession matrix</param>
        /// <param name="rbp">double[3][3]    bias-precession matrix</param>
        /// <param name="rn">double[3][3]    nutation matrix</param>
        /// <param name="rbpn">double[3][3]    GCRS-to-true matrix</param>
        public static void iauPn06a(double date1, double date2,
              out double dpsi, out double deps, out double epsa,
             out double[][] rb, out double[][] rp, out double[][] rbp,
             out double[][] rn, out double[][] rbpn)
        {
            /* Nutation. */
            iauNut06a(date1, date2, out dpsi, out deps);

            /* Remaining results. */
            iauPn06(date1, date2, dpsi, deps, out epsa, out rb, out rp, out rbp, out rn, out rbpn);
        }


        /// <summary>
        /// Form the matrix of precession-nutation for a given date (including
        /// frame bias), equinox based, IAU 2000A model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rbpn">double[3][3] bias-precession-nutation matrix</param>
        public static void iauPnm00a(double date1, double date2, out double[][] rbpn)
        {
            double dpsi, deps, epsa;


            /* Obtain the required matrix (discarding other results). */
            iauPn00a(date1, date2, out dpsi, out deps, out epsa,
                out double[][] rb, out double[][] rp, out double[][] rbp, out double[][] rn, out rbpn);
        }


        /// <summary>
        /// Form the matrix of precession-nutation for a given date (including
        /// frame bias), equinox-based, IAU 2000B model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rbpn">double[3][3] bias-precession-nutation matrix</param>
        public static void iauPnm00b(double date1, double date2, out double[][] rbpn)
        {
            double dpsi, deps, epsa;


            /* Obtain the required matrix (discarding other results). */
            iauPn00b(date1, date2, out dpsi, out deps, out epsa,
                out double[][] rb, out double[][] rp, out double[][] rbp, out double[][] rn, out rbpn);
        }


        /// <summary>
        /// Form the matrix of precession-nutation for a given date (including
        /// frame bias), equinox based, IAU 2006 precession and IAU 2000A
        /// nutation models.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rbpn">double[3][3] bias-precession-nutation matrix</param>
        public static void iauPnm06a(double date1, double date2, out double[][] rbpn)
        {
            double gamb, phib, psib, epsa, dp, de;


            /* Fukushima-Williams angles for frame bias and precession. */
            iauPfw06(date1, date2, out gamb, out phib, out psib, out epsa);

            /* Nutation components. */
            iauNut06a(date1, date2, out dp, out de);

            /* Equinox based nutation x precession x bias matrix. */
            iauFw2m(gamb, phib, psib + dp, epsa + de, out rbpn);
        }


        /// <summary>
        /// Form the matrix of precession/nutation for a given date, IAU 1976
        /// precession model, IAU 1980 nutation model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="rmatpn">double[3][3]   combined precession/nutation matrix</param>
        public static void iauPnm80(double date1, double date2, out double[][] rmatpn)
        {
            double[][] rmatp, rmatn;

            rmatpn = new double[][] { new double[3], new double[3], new double[3] };

            /* Precession matrix, J2000.0 to date. */
            iauPmat76(date1, date2, out rmatp);

            /* Nutation matrix. */
            iauNutm80(date1, date2, out rmatn);

            /* Combine the matrices:  PN = N x P. */
            iauRxr(rmatn, rmatp, ref rmatpn);
        }


        /// <summary>
        /// Form the matrix of polar motion for a given date, IAU 2000.
        /// </summary>
        /// <param name="xp">coordinates of the pole (radians)</param>
        /// <param name="yp">coordinates of the pole (radians)</param>
        /// <param name="sp">the TIO locator s' (radians)</param>
        /// <param name="rpom">double[3][3]   polar-motion matrix</param>
        public static void iauPom00(double xp, double yp, double sp, out double[][] rpom)
        {
            /* Construct the matrix. */
            iauIr(out rpom);
            iauRz(sp, rpom);
            iauRy(-xp, rpom);
            iauRx(-yp, rpom);
        }


        /// <summary>
        /// Precession-rate part of the IAU 2000 precession-nutation models
        /// (part of MHB2000).
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="dpsipr">precession corrections in longitude</param>
        /// <param name="depspr">precession corrections in obliquity</param>
        public static void iauPr00(double date1, double date2, out double dpsipr, out double depspr)
        {
            double t;

            /* Precession and obliquity corrections (radians per century) */
            double PRECOR = -0.29965 * DAS2R, OBLCOR = -0.02524 * DAS2R;


            /* Interval between fundamental epoch J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Precession rate contributions with respect to IAU 1976/80. */
            dpsipr = PRECOR * t;
            depspr = OBLCOR * t;
        }


        /// <summary>
        /// IAU 1976 precession model.
        /// <para>R_3(-z) x R_2(+theta) x R_3(-zeta).</para>
        /// </summary>
        /// <param name="date01">TDB starting date</param>
        /// <param name="date02">TDB starting date</param>
        /// <param name="date11">TDB ending date</param>
        /// <param name="date12">TDB ending date</param>
        /// <param name="zeta">1st rotation: radians cw around z</param>
        /// <param name="z">3rd rotation: radians cw around z</param>
        /// <param name="theta">2nd rotation: radians ccw around y</param>
        public static void iauPrec76(double date01, double date02, double date11, double date12,
               out double zeta, out double z, out double theta)
        {
            double t0, t, tas2r, w;


            /* Interval between fundamental epoch J2000.0 and start date (JC). */
            t0 = ((date01 - DJ00) + date02) / DJC;

            /* Interval over which precession required (JC). */
            t = ((date11 - date01) + (date12 - date02)) / DJC;

            /* Euler angles. */
            tas2r = t * DAS2R;
            w = 2306.2181 + (1.39656 - 0.000139 * t0) * t0;

            zeta = (w + ((0.30188 - 0.000344 * t0) + 0.017998 * t) * t) * tas2r;

            z = (w + ((1.09468 + 0.000066 * t0) + 0.018203 * t) * t) * tas2r;

            theta = ((2004.3109 + (-0.85330 - 0.000217 * t0) * t0)
                   + ((-0.42665 - 0.000217 * t0) - 0.041833 * t) * t) * tas2r;
        }


        /* --------------------- */
        /* The series for s+XY/2 */
        /* --------------------- */

        readonly struct S00TERM
        {
            public readonly int[] nfa;      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
            public readonly double s, c;     /* sine and cosine coefficients */

            public S00TERM(int[] nfa, double s, double c)
            {
                this.nfa = nfa;
                this.s = s;
                this.c = c;
            }
        }

        /* Terms of order t^0 */
        static readonly S00TERM[] S00s0 = {
            /* 1-10 */
            new S00TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 ),
            new S00TERM( new int[] { 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 ),
            new S00TERM( new int[] { 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 ),
            new S00TERM( new int[] { 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 ),
            new S00TERM( new int[] { 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 ),
            new S00TERM( new int[] { 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 ),

            /* 11-20 */
            new S00TERM( new int[] { 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 ),
            new S00TERM( new int[] { 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 ),

            /* 21-30 */
            new S00TERM( new int[] { 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 ),
            new S00TERM( new int[] { 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 ),
            new S00TERM( new int[] { 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 ),

            /* 31-33 */
            new S00TERM( new int[] { 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 )
        };

        /* Terms of order t^1 */
        static readonly S00TERM[] S00s1 ={
            /* 1-3 */
            new S00TERM( new int[] { 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 ),
            new S00TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},     1.71e-6,  -0.03e-6 ),
            new S00TERM( new int[] { 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 )
        };

        /* Terms of order t^2 */
        static readonly S00TERM[] S00s2 ={
            /* 1-10 */
            new S00TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},   743.53e-6,  -0.17e-6 ),
            new S00TERM( new int[] { 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 ),
            new S00TERM( new int[] { 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 ),
            new S00TERM( new int[] { 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 ),
            new S00TERM( new int[] { 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 ),
            new S00TERM( new int[] { 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 ),

            /* 11-20 */
            new S00TERM( new int[] { 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 ),
            new S00TERM( new int[] { 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 ),
            new S00TERM( new int[] { 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 ),

            /* 21-25 */
            new S00TERM( new int[] { 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 ),
            new S00TERM( new int[] { 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 ),
            new S00TERM( new int[] { 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 ),
            new S00TERM( new int[] { 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 ),
            new S00TERM( new int[] { 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 )
        };

        /* Terms of order t^3 */
        static readonly S00TERM[] S00s3 ={
            /* 1-4 */
            new S00TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.51e-6 ),
            new S00TERM( new int[] { 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.39e-6 ),
            new S00TERM(new int[] { 0, 0, 2, 0, 2, 0, 0, 0 }, -0.01e-6, -0.24e-6 ),
            new S00TERM(new int[] { 0, 0, 0, 0, 2, 0, 0, 0 }, 0.00e-6, 0.22e-6 )
        };

        /* Terms of order t^4 */
        static readonly S00TERM[] S00s4 ={
            /* 1-1 */
            new S00TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 )
        };

        /// <summary>
        /// The CIO locator s, positioning the Celestial Intermediate Origin on
        /// the equator of the Celestial Intermediate Pole, given the CIP's X,Y
        /// coordinates.  Compatible with IAU 2000A precession-nutation.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="x">CIP coordinates</param>
        /// <param name="y">CIP coordinates</param>
        /// <returns>the CIO locator s in radians</returns>
        public static double iauS00(double date1, double date2, double x, double y)
        {
            /* Time since J2000.0, in Julian centuries */
            double t;

            /* Miscellaneous */
            int i, j;
            double a, w0, w1, w2, w3, w4, w5;

            /* Fundamental arguments */
            double[] fa = new double[8];

            /* Returned value */
            double s;


            /* Polynomial coefficients */
            double[] sp = {
                /* 1-6 */
                94.00e-6,
                3808.35e-6,
                -119.94e-6,
                -72574.09e-6,
                27.70e-6,
                15.61e-6
            };


            /* Number of terms in the series */
            //const int NS0 = (int)(sizeof s0 / sizeof(TERM));
            //const int NS1 = (int)(sizeof s1 / sizeof(TERM));
            //const int NS2 = (int)(sizeof s2 / sizeof(TERM));
            //const int NS3 = (int)(sizeof s3 / sizeof(TERM));
            //const int NS4 = (int)(sizeof s4 / sizeof(TERM));
            int NS0 = S00s0.Length;
            int NS1 = S00s1.Length;
            int NS2 = S00s2.Length;
            int NS3 = S00s3.Length;
            int NS4 = S00s4.Length;


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

            /* Evaluate s. */
            w0 = sp[0];
            w1 = sp[1];
            w2 = sp[2];
            w3 = sp[3];
            w4 = sp[4];
            w5 = sp[5];

            for (i = NS0 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S00s0[i].nfa[j] * fa[j];
                }
                w0 += S00s0[i].s * Math.Sin(a) + S00s0[i].c * Math.Cos(a);
            }

            for (i = NS1 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S00s1[i].nfa[j] * fa[j];
                }
                w1 += S00s1[i].s * Math.Sin(a) + S00s1[i].c * Math.Cos(a);
            }

            for (i = NS2 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S00s2[i].nfa[j] * fa[j];
                }
                w2 += S00s2[i].s * Math.Sin(a) + S00s2[i].c * Math.Cos(a);
            }

            for (i = NS3 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S00s3[i].nfa[j] * fa[j];
                }
                w3 += S00s3[i].s * Math.Sin(a) + S00s3[i].c * Math.Cos(a);
            }

            for (i = NS4 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S00s4[i].nfa[j] * fa[j];
                }
                w4 += S00s4[i].s * Math.Sin(a) + S00s4[i].c * Math.Cos(a);
            }

            s = (w0 +
                (w1 +
                (w2 +
                (w3 +
                (w4 +
                 w5 * t) * t) * t) * t) * t) * DAS2R - x * y / 2.0;

            return s;
        }


        /// <summary>
        /// The CIO locator s, positioning the Celestial Intermediate Origin on
        /// the equator of the Celestial Intermediate Pole, using the IAU 2000A
        /// precession-nutation model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>the CIO locator s in radians</returns>
        public static double iauS00a(double date1, double date2)
        {
            double[][] rbpn;
            double x, y, s;


            /* Bias-precession-nutation-matrix, IAU 2000A. */
            iauPnm00a(date1, date2, out rbpn);

            /* Extract the CIP coordinates. */
            iauBpn2xy(rbpn, out x, out y);

            /* Compute the CIO locator s, given the CIP coordinates. */
            s = iauS00(date1, date2, x, y);

            return s;
        }


        /// <summary>
        /// The CIO locator s, positioning the Celestial Intermediate Origin on
        /// the equator of the Celestial Intermediate Pole, using the IAU 2000B
        /// precession-nutation model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>the CIO locator s in radians</returns>
        public static double iauS00b(double date1, double date2)
        {
            double[][] rbpn;
            double x, y, s;


            /* Bias-precession-nutation-matrix, IAU 2000B. */
            iauPnm00b(date1, date2, out rbpn);

            /* Extract the CIP coordinates. */
            iauBpn2xy(rbpn, out x, out y);

            /* Compute the CIO locator s, given the CIP coordinates. */
            s = iauS00(date1, date2, x, y);

            return s;
        }


        /* --------------------- */
        /* The series for s+XY/2 */
        /* --------------------- */

        readonly struct S06TERM
        {
            public readonly int[] nfa;      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
            public readonly double s, c;     /* sine and cosine coefficients */

            public S06TERM(int[] nfa, double s, double c)
            {
                this.nfa = nfa;
                this.s = s;
                this.c = c;
            }
        }

        /* Terms of order t^0 */
        static readonly S06TERM[] S06s0 = {

            /* 1-10 */
            new S06TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 ),
            new S06TERM( new int[] { 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 ),
            new S06TERM( new int[] { 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 ),
            new S06TERM( new int[] { 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 ),
            new S06TERM( new int[] { 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 ),

            /* 11-20 */
            new S06TERM( new int[] { 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 ),

            /* 21-30 */
            new S06TERM( new int[] { 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 ),
            new S06TERM( new int[] { 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 ),

            /* 31-33 */
            new S06TERM( new int[] { 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 )
        };

        /* Terms of order t^1 */
        static readonly S06TERM[] S06s1 = {

            /* 1 - 3 */
            new S06TERM( new int[] { 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},     1.73e-6,  -0.03e-6 ),
            new S06TERM( new int[] { 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 )
        };

        /* Terms of order t^2 */
        static readonly S06TERM[] S06s2 = {

            /* 1-10 */
            new S06TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},   743.52e-6,  -0.17e-6 ),
            new S06TERM( new int[] { 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 ),
            new S06TERM( new int[] { 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 ),
            new S06TERM( new int[] { 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 ),

            /* 11-20 */
            new S06TERM( new int[] { 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 ),
            new S06TERM( new int[] { 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 ),
            new S06TERM( new int[] { 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 ),

            /* 21-25 */
            new S06TERM( new int[] { 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 ),
            new S06TERM( new int[] { 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 ),
            new S06TERM( new int[] { 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 ),
            new S06TERM( new int[] { 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 )
        };

        /* Terms of order t^3 */
        static readonly S06TERM[] S06s3 = {

            /* 1-4 */
            new S06TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.42e-6 ),
            new S06TERM( new int[] { 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.46e-6 ),
            new S06TERM( new int[] { 0,  0,  2,  0,  2,  0,  0,  0},    -0.01e-6,  -0.25e-6 ),
            new S06TERM( new int[] { 0,  0,  0,  0,  2,  0,  0,  0},     0.00e-6,   0.23e-6 )
        };

        /* Terms of order t^4 */
        static readonly S06TERM[] S06s4 = {

            /* 1-1 */
            new S06TERM( new int[] { 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 )
        };

        /// <summary>
        /// The CIO locator s, positioning the Celestial Intermediate Origin on
        /// the equator of the Celestial Intermediate Pole, given the CIP's X,Y
        /// coordinates.  Compatible with IAU 2006/2000A precession-nutation.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="x">CIP coordinates</param>
        /// <param name="y">CIP coordinates</param>
        /// <returns>the CIO locator s in radians</returns>
        public static double iauS06(double date1, double date2, double x, double y)
        {
            /* Time since J2000.0, in Julian centuries */
            double t;

            /* Miscellaneous */
            int i, j;
            double a, w0, w1, w2, w3, w4, w5;

            /* Fundamental arguments */
            double[] fa = new double[8];

            /* Returned value */
            double s;


            /* Polynomial coefficients */
            double[] sp = {
                /* 1-6 */
                94.00e-6,
                3808.65e-6,
                -122.68e-6,
                -72574.11e-6,
                27.98e-6,
                15.62e-6
            };


            /* Number of terms in the series */
            //static const int NS0 = (int)(sizeof s0 / sizeof(TERM));
            //static const int NS1 = (int)(sizeof s1 / sizeof(TERM));
            //static const int NS2 = (int)(sizeof s2 / sizeof(TERM));
            //static const int NS3 = (int)(sizeof s3 / sizeof(TERM));
            //static const int NS4 = (int)(sizeof s4 / sizeof(TERM));
            int NS0 = S06s0.Length;
            int NS1 = S06s1.Length;
            int NS2 = S06s2.Length;
            int NS3 = S06s3.Length;
            int NS4 = S06s4.Length;


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

            /* Evaluate s. */
            w0 = sp[0];
            w1 = sp[1];
            w2 = sp[2];
            w3 = sp[3];
            w4 = sp[4];
            w5 = sp[5];

            for (i = NS0 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S06s0[i].nfa[j] * fa[j];
                }
                w0 += S06s0[i].s * Math.Sin(a) + S06s0[i].c * Math.Cos(a);
            }

            for (i = NS1 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S06s1[i].nfa[j] * fa[j];
                }
                w1 += S06s1[i].s * Math.Sin(a) + S06s1[i].c * Math.Cos(a);
            }

            for (i = NS2 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S06s2[i].nfa[j] * fa[j];
                }
                w2 += S06s2[i].s * Math.Sin(a) + S06s2[i].c * Math.Cos(a);
            }

            for (i = NS3 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S06s3[i].nfa[j] * fa[j];
                }
                w3 += S06s3[i].s * Math.Sin(a) + S06s3[i].c * Math.Cos(a);
            }

            for (i = NS4 - 1; i >= 0; i--)
            {
                a = 0.0;
                for (j = 0; j < 8; j++)
                {
                    a += (double)S06s4[i].nfa[j] * fa[j];
                }
                w4 += S06s4[i].s * Math.Sin(a) + S06s4[i].c * Math.Cos(a);
            }

            s = (w0 +
                (w1 +
                (w2 +
                (w3 +
                (w4 +
                 w5 * t) * t) * t) * t) * t) * DAS2R - x * y / 2.0;

            return s;
        }


        /// <summary>
        /// The CIO locator s, positioning the Celestial Intermediate Origin on
        /// the equator of the Celestial Intermediate Pole, using the IAU 2006
        /// precession and IAU 2000A nutation models.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>the CIO locator s in radians</returns>
        public static double iauS06a(double date1, double date2)
        {
            double[][] rnpb;
            double x, y, s;


            /* Bias-precession-nutation-matrix, IAU 20006/2000A. */
            iauPnm06a(date1, date2, out rnpb);

            /* Extract the CIP coordinates. */
            iauBpn2xy(rnpb, out x, out y);

            /* Compute the CIO locator s, given the CIP coordinates. */
            s = iauS06(date1, date2, x, y);

            return s;
        }


        /// <summary>
        /// The TIO locator s', positioning the Terrestrial Intermediate Origin
        /// on the equator of the Celestial Intermediate Pole.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <returns>the TIO locator s' in radians</returns>
        public static double iauSp00(double date1, double date2)
        {
            double t, sp;


            /* Interval between fundamental epoch J2000.0 and current date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Approximate s'. */
            sp = -47e-6 * t * DAS2R;

            return sp;
        }


        /* Polynomial coefficients (arcsec, X then Y). */
        static readonly double[][] Xy06xyp = {
                new double[] {
                    -0.016617,
                    2004.191898,
                    -0.4297829,
                    -0.19861834,
                    0.000007578,
                    0.0000059285
                },
                new double[] {
                    -0.006951,
                    -0.025896,
                    -22.4072747,
                    0.00190059,
                    0.001112526,
                    0.0000001358
                }
            };

        /* Fundamental-argument multipliers:  luni-solar terms */
        static readonly int[][] Xy06mfals = {
            /* 1-10 */
            new int[] {  0,   0,   0,   0,   1 },
            new int[] {  0,   0,   2,  -2,   2 },
            new int[] {  0,   0,   2,   0,   2 },
            new int[] {  0,   0,   0,   0,   2 },
            new int[] {  0,   1,   0,   0,   0 },
            new int[] {  0,   1,   2,  -2,   2 },
            new int[] {  1,   0,   0,   0,   0 },
            new int[] {  0,   0,   2,   0,   1 },
            new int[] {  1,   0,   2,   0,   2 },
            new int[] {  0,   1,  -2,   2,  -2 },

            /* 11-20 */
            new int[] {  0,   0,   2,  -2,   1 },
            new int[] {  1,   0,  -2,   0,  -2 },
            new int[] {  1,   0,   0,  -2,   0 },
            new int[] {  1,   0,   0,   0,   1 },
            new int[] {  1,   0,   0,   0,  -1 },
            new int[] {  1,   0,  -2,  -2,  -2 },
            new int[] {  1,   0,   2,   0,   1 },
            new int[] {  2,   0,  -2,   0,  -1 },
            new int[] {  0,   0,   0,   2,   0 },
            new int[] {  0,   0,   2,   2,   2 },

            /* 21-30 */
            new int[] {  2,   0,   0,  -2,   0 },
            new int[] {  0,   2,  -2,   2,  -2 },
            new int[] {  2,   0,   2,   0,   2 },
            new int[] {  1,   0,   2,  -2,   2 },
            new int[] {  1,   0,  -2,   0,  -1 },
            new int[] {  2,   0,   0,   0,   0 },
            new int[] {  0,   0,   2,   0,   0 },
            new int[] {  0,   1,   0,   0,   1 },
            new int[] {  1,   0,   0,  -2,  -1 },
            new int[] {  0,   2,   2,  -2,   2 },

            /* 31-40 */
            new int[] {  0,   0,   2,  -2,   0 },
            new int[] {  1,   0,   0,  -2,   1 },
            new int[] {  0,   1,   0,   0,  -1 },
            new int[] {  0,   2,   0,   0,   0 },
            new int[] {  1,   0,  -2,  -2,  -1 },
            new int[] {  1,   0,   2,   2,   2 },
            new int[] {  0,   1,   2,   0,   2 },
            new int[] {  2,   0,  -2,   0,   0 },
            new int[] {  0,   0,   2,   2,   1 },
            new int[] {  0,   1,  -2,   0,  -2 },

            /* 41-50 */
            new int[] {  0,   0,   0,   2,   1 },
            new int[] {  1,   0,   2,  -2,   1 },
            new int[] {  2,   0,   0,  -2,  -1 },
            new int[] {  2,   0,   2,  -2,   2 },
            new int[] {  2,   0,   2,   0,   1 },
            new int[] {  0,   0,   0,   2,  -1 },
            new int[] {  0,   1,  -2,   2,  -1 },
            new int[] {  1,   1,   0,  -2,   0 },
            new int[] {  2,   0,   0,  -2,   1 },
            new int[] {  1,   0,   0,   2,   0 },

            /* 51-60 */
            new int[] {  0,   1,   2,  -2,   1 },
            new int[] {  1,  -1,   0,   0,   0 },
            new int[] {  0,   1,  -1,   1,  -1 },
            new int[] {  2,   0,  -2,   0,  -2 },
            new int[] {  0,   1,   0,  -2,   0 },
            new int[] {  1,   0,   0,  -1,   0 },
            new int[] {  3,   0,   2,   0,   2 },
            new int[] {  0,   0,   0,   1,   0 },
            new int[] {  1,  -1,   2,   0,   2 },
            new int[] {  1,   1,  -2,  -2,  -2 },

            /* 61-70 */
            new int[] {  1,   0,  -2,   0,   0 },
            new int[] {  2,   0,   0,   0,  -1 },
            new int[] {  0,   1,  -2,  -2,  -2 },
            new int[] {  1,   1,   2,   0,   2 },
            new int[] {  2,   0,   0,   0,   1 },
            new int[] {  1,   1,   0,   0,   0 },
            new int[] {  1,   0,  -2,   2,  -1 },
            new int[] {  1,   0,   2,   0,   0 },
            new int[] {  1,  -1,   0,  -1,   0 },
            new int[] {  1,   0,   0,   0,   2 },

            /* 71-80 */
            new int[] {  1,   0,  -1,   0,  -1 },
            new int[] {  0,   0,   2,   1,   2 },
            new int[] {  1,   0,  -2,  -4,  -2 },
            new int[] {  1,  -1,   0,  -1,  -1 },
            new int[] {  1,   0,   2,   2,   1 },
            new int[] {  0,   2,  -2,   2,  -1 },
            new int[] {  1,   0,   0,   0,  -2 },
            new int[] {  2,   0,  -2,  -2,  -2 },
            new int[] {  1,   1,   2,  -2,   2 },
            new int[] {  2,   0,  -2,  -4,  -2 },

            /* 81-90 */
            new int[] {  1,   0,  -4,   0,  -2 },
            new int[] {  2,   0,   2,  -2,   1 },
            new int[] {  1,   0,   0,  -1,  -1 },
            new int[] {  2,   0,   2,   2,   2 },
            new int[] {  3,   0,   0,   0,   0 },
            new int[] {  1,   0,   0,   2,   1 },
            new int[] {  0,   0,   2,  -2,  -1 },
            new int[] {  3,   0,   2,  -2,   2 },
            new int[] {  0,   0,   4,  -2,   2 },
            new int[] {  1,   0,   0,  -4,   0 },

            /* 91-100 */
            new int[] {  0,   1,   2,   0,   1 },
            new int[] {  2,   0,   0,  -4,   0 },
            new int[] {  1,   1,   0,  -2,  -1 },
            new int[] {  2,   0,  -2,   0,   1 },
            new int[] {  0,   0,   2,   0,  -1 },
            new int[] {  0,   1,  -2,   0,  -1 },
            new int[] {  0,   1,   0,   0,   2 },
            new int[] {  0,   0,   2,  -1,   2 },
            new int[] {  0,   0,   2,   4,   2 },
            new int[] {  2,   1,   0,  -2,   0 },

            /* 101-110 */
            new int[] {  1,   1,   0,  -2,   1 },
            new int[] {  1,  -1,   0,  -2,   0 },
            new int[] {  1,  -1,   0,  -1,  -2 },
            new int[] {  1,  -1,   0,   0,   1 },
            new int[] {  0,   1,  -2,   2,   0 },
            new int[] {  0,   1,   0,   0,  -2 },
            new int[] {  1,  -1,   2,   2,   2 },
            new int[] {  1,   0,   0,   2,  -1 },
            new int[] {  1,  -1,  -2,  -2,  -2 },
            new int[] {  3,   0,   2,   0,   1 },

            /* 111-120 */
            new int[] {  0,   1,   2,   2,   2 },
            new int[] {  1,   0,   2,  -2,   0 },
            new int[] {  1,   1,  -2,  -2,  -1 },
            new int[] {  1,   0,   2,  -4,   1 },
            new int[] {  0,   1,  -2,  -2,  -1 },
            new int[] {  2,  -1,   2,   0,   2 },
            new int[] {  0,   0,   0,   2,   2 },
            new int[] {  1,  -1,   2,   0,   1 },
            new int[] {  1,  -1,  -2,   0,  -2 },
            new int[] {  0,   1,   0,   2,   0 },

            /* 121-130 */
            new int[] {  0,   1,   2,  -2,   0 },
            new int[] {  0,   0,   0,   1,   1 },
            new int[] {  1,   0,  -2,  -2,   0 },
            new int[] {  0,   3,   2,  -2,   2 },
            new int[] {  2,   1,   2,   0,   2 },
            new int[] {  1,   1,   0,   0,   1 },
            new int[] {  2,   0,   0,   2,   0 },
            new int[] {  1,   1,   2,   0,   1 },
            new int[] {  1,   0,   0,  -2,  -2 },
            new int[] {  1,   0,  -2,   2,   0 },

            /* 131-140 */
            new int[] {  1,   0,  -1,   0,  -2 },
            new int[] {  0,   1,   0,  -2,   1 },
            new int[] {  0,   1,   0,   1,   0 },
            new int[] {  0,   0,   0,   1,  -1 },
            new int[] {  1,   0,  -2,   2,  -2 },
            new int[] {  1,  -1,   0,   0,  -1 },
            new int[] {  0,   0,   0,   4,   0 },
            new int[] {  1,  -1,   0,   2,   0 },
            new int[] {  1,   0,   2,   1,   2 },
            new int[] {  1,   0,   2,  -1,   2 },

            /* 141-150 */
            new int[] {  0,   0,   2,   1,   1 },
            new int[] {  1,   0,   0,  -2,   2 },
            new int[] {  1,   0,  -2,   0,   1 },
            new int[] {  1,   0,  -2,  -4,  -1 },
            new int[] {  0,   0,   2,   2,   0 },
            new int[] {  1,   1,   2,  -2,   1 },
            new int[] {  1,   0,  -2,   1,  -1 },
            new int[] {  0,   0,   1,   0,   1 },
            new int[] {  2,   0,  -2,  -2,  -1 },
            new int[] {  4,   0,   2,   0,   2 },

            /* 151-160 */
            new int[] {  2,  -1,   0,   0,   0 },
            new int[] {  2,   1,   2,  -2,   2 },
            new int[] {  0,   1,   2,   1,   2 },
            new int[] {  1,   0,   4,  -2,   2 },
            new int[] {  1,   1,   0,   0,  -1 },
            new int[] {  2,   0,   2,   0,   0 },
            new int[] {  2,   0,  -2,  -4,  -1 },
            new int[] {  1,   0,  -1,   0,   0 },
            new int[] {  1,   0,   0,   1,   0 },
            new int[] {  0,   1,   0,   2,   1 },

            /* 161-170 */
            new int[] {  1,   0,  -4,   0,  -1 },
            new int[] {  1,   0,   0,  -4,  -1 },
            new int[] {  2,   0,   2,   2,   1 },
            new int[] {  2,   1,   0,   0,   0 },
            new int[] {  0,   0,   2,  -3,   2 },
            new int[] {  1,   2,   0,  -2,   0 },
            new int[] {  0,   3,   0,   0,   0 },
            new int[] {  0,   0,   4,   0,   2 },
            new int[] {  0,   0,   2,  -4,   1 },
            new int[] {  2,   0,   0,  -2,  -2 },

            /* 171-180 */
            new int[] {  1,   1,  -2,  -4,  -2 },
            new int[] {  0,   1,   0,  -2,  -1 },
            new int[] {  0,   0,   0,   4,   1 },
            new int[] {  3,   0,   2,  -2,   1 },
            new int[] {  1,   0,   2,   4,   2 },
            new int[] {  1,   1,  -2,   0,  -2 },
            new int[] {  0,   0,   4,  -2,   1 },
            new int[] {  2,  -2,   0,  -2,   0 },
            new int[] {  2,   1,   0,  -2,  -1 },
            new int[] {  0,   2,   0,  -2,   0 },

            /* 181-190 */
            new int[] {  1,   0,   0,  -1,   1 },
            new int[] {  1,   1,   2,   2,   2 },
            new int[] {  3,   0,   0,   0,  -1 },
            new int[] {  2,   0,   0,  -4,  -1 },
            new int[] {  3,   0,   2,   2,   2 },
            new int[] {  0,   0,   2,   4,   1 },
            new int[] {  0,   2,  -2,  -2,  -2 },
            new int[] {  1,  -1,   0,  -2,  -1 },
            new int[] {  0,   0,   2,  -1,   1 },
            new int[] {  2,   0,   0,   2,   1 },

            /* 191-200 */
            new int[] {  1,  -1,  -2,   2,  -1 },
            new int[] {  0,   0,   0,   2,  -2 },
            new int[] {  2,   0,   0,  -4,   1 },
            new int[] {  1,   0,   0,  -4,   1 },
            new int[] {  2,   0,   2,  -4,   1 },
            new int[] {  4,   0,   2,  -2,   2 },
            new int[] {  2,   1,  -2,   0,  -1 },
            new int[] {  2,   1,  -2,  -4,  -2 },
            new int[] {  3,   0,   0,  -4,   0 },
            new int[] {  1,  -1,   2,   2,   1 },

            /* 201-210 */
            new int[] {  1,  -1,  -2,   0,  -1 },
            new int[] {  0,   2,   0,   0,   1 },
            new int[] {  1,   2,  -2,  -2,  -2 },
            new int[] {  1,   1,   0,  -4,   0 },
            new int[] {  2,   0,   0,  -2,   2 },
            new int[] {  0,   2,   2,  -2,   1 },
            new int[] {  1,   0,   2,   0,  -1 },
            new int[] {  2,   1,   0,  -2,   1 },
            new int[] {  2,  -1,  -2,   0,  -1 },
            new int[] {  1,  -1,  -2,  -2,  -1 },

            /* 211-220 */
            new int[] {  0,   1,  -2,   1,  -2 },
            new int[] {  1,   0,  -4,   2,  -2 },
            new int[] {  0,   1,   2,   2,   1 },
            new int[] {  3,   0,   0,   0,   1 },
            new int[] {  2,  -1,   2,   2,   2 },
            new int[] {  0,   1,  -2,  -4,  -2 },
            new int[] {  1,   0,  -2,  -3,  -2 },
            new int[] {  2,   0,   0,   0,   2 },
            new int[] {  1,  -1,   0,  -2,  -2 },
            new int[] {  2,   0,  -2,   2,  -1 },

            /* 221-230 */
            new int[] {  0,   2,  -2,   0,  -2 },
            new int[] {  3,   0,  -2,   0,  -1 },
            new int[] {  2,  -1,   2,   0,   1 },
            new int[] {  1,   0,  -2,  -1,  -2 },
            new int[] {  0,   0,   2,   0,   3 },
            new int[] {  2,   0,  -4,   0,  -2 },
            new int[] {  2,   1,   0,  -4,   0 },
            new int[] {  1,   1,  -2,   1,  -1 },
            new int[] {  0,   2,   2,   0,   2 },
            new int[] {  1,  -1,   2,  -2,   2 },

            /* 231-240 */
            new int[] {  1,  -1,   0,  -2,   1 },
            new int[] {  2,   1,   2,   0,   1 },
            new int[] {  1,   0,   2,  -4,   2 },
            new int[] {  1,   1,  -2,   0,  -1 },
            new int[] {  1,   1,   0,   2,   0 },
            new int[] {  1,   0,   0,  -3,   0 },
            new int[] {  2,   0,   2,  -1,   2 },
            new int[] {  0,   2,   0,   0,  -1 },
            new int[] {  2,  -1,   0,  -2,   0 },
            new int[] {  4,   0,   0,   0,   0 },

            /* 241-250 */
            new int[] {  2,   1,  -2,  -2,  -2 },
            new int[] {  0,   2,  -2,   2,   0 },
            new int[] {  1,   0,   2,   1,   1 },
            new int[] {  1,   0,  -1,   0,  -3 },
            new int[] {  3,  -1,   2,   0,   2 },
            new int[] {  2,   0,   2,  -2,   0 },
            new int[] {  1,  -2,   0,   0,   0 },
            new int[] {  2,   0,   0,   0,  -2 },
            new int[] {  1,   0,   0,   4,   0 },
            new int[] {  0,   1,   0,   1,   1 },

            /* 251-260 */
            new int[] {  1,   0,   2,   2,   0 },
            new int[] {  0,   1,   0,   2,  -1 },
            new int[] {  0,   1,   0,   1,  -1 },
            new int[] {  0,   0,   2,  -2,   3 },
            new int[] {  3,   1,   2,   0,   2 },
            new int[] {  1,   1,   2,   1,   2 },
            new int[] {  1,   1,  -2,   2,  -1 },
            new int[] {  2,  -1,   2,  -2,   2 },
            new int[] {  1,  -2,   2,   0,   2 },
            new int[] {  1,   0,   2,  -4,   0 },

            /* 261-270 */
            new int[] {  0,   0,   1,   0,   0 },
            new int[] {  1,   0,   2,  -3,   1 },
            new int[] {  1,  -2,   0,  -2,   0 },
            new int[] {  2,   0,   0,   2,  -1 },
            new int[] {  1,   1,   2,  -4,   1 },
            new int[] {  4,   0,   2,   0,   1 },
            new int[] {  0,   1,   2,   1,   1 },
            new int[] {  1,   2,   2,  -2,   2 },
            new int[] {  2,   0,   2,   1,   2 },
            new int[] {  2,   1,   2,  -2,   1 },

            /* 271-280 */
            new int[] {  1,   0,   2,  -1,   1 },
            new int[] {  1,   0,   4,  -2,   1 },
            new int[] {  1,  -1,   2,  -2,   1 },
            new int[] {  0,   1,   0,  -4,   0 },
            new int[] {  3,   0,  -2,  -2,  -2 },
            new int[] {  0,   0,   4,  -4,   2 },
            new int[] {  2,   0,  -4,  -2,  -2 },
            new int[] {  2,  -2,   0,  -2,  -1 },
            new int[] {  1,   0,   2,  -2,  -1 },
            new int[] {  2,   0,  -2,  -6,  -2 },

            /* 281-290 */
            new int[] {  1,   0,  -2,   1,  -2 },
            new int[] {  1,   0,  -2,   2,   1 },
            new int[] {  1,  -1,   0,   2,  -1 },
            new int[] {  1,   0,  -2,   1,   0 },
            new int[] {  2,  -1,   0,  -2,   1 },
            new int[] {  1,  -1,   0,   2,   1 },
            new int[] {  2,   0,  -2,  -2,   0 },
            new int[] {  1,   0,   2,  -3,   2 },
            new int[] {  0,   0,   0,   4,  -1 },
            new int[] {  2,  -1,   0,   0,   1 },

            /* 291-300 */
            new int[] {  2,   0,   4,  -2,   2 },
            new int[] {  0,   0,   2,   3,   2 },
            new int[] {  0,   1,   4,  -2,   2 },
            new int[] {  0,   1,  -2,   2,   1 },
            new int[] {  1,   1,   0,   2,   1 },
            new int[] {  1,   0,   0,   4,   1 },
            new int[] {  0,   0,   4,   0,   1 },
            new int[] {  2,   0,   0,  -3,   0 },
            new int[] {  1,   0,   0,  -1,  -2 },
            new int[] {  1,  -2,  -2,  -2,  -2 },

            /* 301-310 */
            new int[] {  3,   0,   0,   2,   0 },
            new int[] {  2,   0,   2,  -4,   2 },
            new int[] {  1,   1,  -2,  -4,  -1 },
            new int[] {  1,   0,  -2,  -6,  -2 },
            new int[] {  2,  -1,   0,   0,  -1 },
            new int[] {  2,  -1,   0,   2,   0 },
            new int[] {  0,   1,   2,  -2,  -1 },
            new int[] {  1,   1,   0,   1,   0 },
            new int[] {  1,   2,   0,  -2,  -1 },
            new int[] {  1,   0,   0,   1,  -1 },

            /* 311-320 */
            new int[] {  0,   0,   1,   0,   2 },
            new int[] {  3,   1,   2,  -2,   2 },
            new int[] {  1,   0,  -4,  -2,  -2 },
            new int[] {  1,   0,   2,   4,   1 },
            new int[] {  1,  -2,   2,   2,   2 },
            new int[] {  1,  -1,  -2,  -4,  -2 },
            new int[] {  0,   0,   2,  -4,   2 },
            new int[] {  0,   0,   2,  -3,   1 },
            new int[] {  2,   1,  -2,   0,   0 },
            new int[] {  3,   0,  -2,  -2,  -1 },

            /* 321-330 */
            new int[] {  2,   0,   2,   4,   2 },
            new int[] {  0,   0,   0,   0,   3 },
            new int[] {  2,  -1,  -2,  -2,  -2 },
            new int[] {  2,   0,   0,  -1,   0 },
            new int[] {  3,   0,   2,  -4,   2 },
            new int[] {  2,   1,   2,   2,   2 },
            new int[] {  0,   0,   3,   0,   3 },
            new int[] {  1,   1,   2,   2,   1 },
            new int[] {  2,   1,   0,   0,  -1 },
            new int[] {  1,   2,   0,  -2,   1 },

            /* 331-340 */
            new int[] {  3,   0,   2,   2,   1 },
            new int[] {  1,  -1,  -2,   2,  -2 },
            new int[] {  1,   1,   0,  -1,   0 },
            new int[] {  1,   2,   0,   0,   0 },
            new int[] {  1,   0,   4,   0,   2 },
            new int[] {  1,  -1,   2,   4,   2 },
            new int[] {  2,   1,   0,   0,   1 },
            new int[] {  1,   0,   0,   2,   2 },
            new int[] {  1,  -1,  -2,   2,   0 },
            new int[] {  0,   2,  -2,  -2,  -1 },

            /* 341-350 */
            new int[] {  2,   0,  -2,   0,   2 },
            new int[] {  5,   0,   2,   0,   2 },
            new int[] {  3,   0,  -2,  -6,  -2 },
            new int[] {  1,  -1,   2,  -1,   2 },
            new int[] {  3,   0,   0,  -4,  -1 },
            new int[] {  1,   0,   0,   1,   1 },
            new int[] {  1,   0,  -4,   2,  -1 },
            new int[] {  0,   1,   2,  -4,   1 },
            new int[] {  1,   2,   2,   0,   2 },
            new int[] {  0,   1,   0,  -2,  -2 },

            /* 351-360 */
            new int[] {  0,   0,   2,  -1,   0 },
            new int[] {  1,   0,   1,   0,   1 },
            new int[] {  0,   2,   0,  -2,   1 },
            new int[] {  3,   0,   2,   0,   0 },
            new int[] {  1,   1,  -2,   1,   0 },
            new int[] {  2,   1,  -2,  -4,  -1 },
            new int[] {  3,  -1,   0,   0,   0 },
            new int[] {  2,  -1,  -2,   0,   0 },
            new int[] {  4,   0,   2,  -2,   1 },
            new int[] {  2,   0,  -2,   2,   0 },

            /* 361-370 */
            new int[] {  1,   1,   2,  -2,   0 },
            new int[] {  1,   0,  -2,   4,  -1 },
            new int[] {  1,   0,  -2,  -2,   1 },
            new int[] {  2,   0,   2,  -4,   0 },
            new int[] {  1,   1,   0,  -2,  -2 },
            new int[] {  1,   1,  -2,  -2,   0 },
            new int[] {  1,   0,   1,  -2,   1 },
            new int[] {  2,  -1,  -2,  -4,  -2 },
            new int[] {  3,   0,  -2,   0,  -2 },
            new int[] {  0,   1,  -2,  -2,   0 },

            /* 371-380 */
            new int[] {  3,   0,   0,  -2,  -1 },
            new int[] {  1,   0,  -2,  -3,  -1 },
            new int[] {  0,   1,   0,  -4,  -1 },
            new int[] {  1,  -2,   2,  -2,   1 },
            new int[] {  0,   1,  -2,   1,  -1 },
            new int[] {  1,  -1,   0,   0,   2 },
            new int[] {  2,   0,   0,   1,   0 },
            new int[] {  1,  -2,   0,   2,   0 },
            new int[] {  1,   2,  -2,  -2,  -1 },
            new int[] {  0,   0,   4,  -4,   1 },

            /* 381-390 */
            new int[] {  0,   1,   2,   4,   2 },
            new int[] {  0,   1,  -4,   2,  -2 },
            new int[] {  3,   0,  -2,   0,   0 },
            new int[] {  2,  -1,   2,   2,   1 },
            new int[] {  0,   1,  -2,  -4,  -1 },
            new int[] {  4,   0,   2,   2,   2 },
            new int[] {  2,   0,  -2,  -3,  -2 },
            new int[] {  2,   0,   0,  -6,   0 },
            new int[] {  1,   0,   2,   0,   3 },
            new int[] {  3,   1,   0,   0,   0 },

            /* 391-400 */
            new int[] {  3,   0,   0,  -4,   1 },
            new int[] {  1,  -1,   2,   0,   0 },
            new int[] {  1,  -1,   0,  -4,   0 },
            new int[] {  2,   0,  -2,   2,  -2 },
            new int[] {  1,   1,   0,  -2,   2 },
            new int[] {  4,   0,   0,  -2,   0 },
            new int[] {  2,   2,   0,  -2,   0 },
            new int[] {  0,   1,   2,   0,   0 },
            new int[] {  1,   1,   0,  -4,   1 },
            new int[] {  1,   0,   0,  -4,  -2 },

            /* 401-410 */
            new int[] {  0,   0,   0,   1,   2 },
            new int[] {  3,   0,   0,   2,   1 },
            new int[] {  1,   1,   0,  -4,  -1 },
            new int[] {  0,   0,   2,   2,  -1 },
            new int[] {  1,   1,   2,   0,   0 },
            new int[] {  1,  -1,   2,  -4,   1 },
            new int[] {  1,   1,   0,   0,   2 },
            new int[] {  0,   0,   2,   6,   2 },
            new int[] {  4,   0,  -2,  -2,  -1 },
            new int[] {  2,   1,   0,  -4,  -1 },

            /* 411-420 */
            new int[] {  0,   0,   0,   3,   1 },
            new int[] {  1,  -1,  -2,   0,   0 },
            new int[] {  0,   0,   2,   1,   0 },
            new int[] {  1,   0,   0,   2,  -2 },
            new int[] {  3,  -1,   2,   2,   2 },
            new int[] {  3,  -1,   2,  -2,   2 },
            new int[] {  1,   0,   0,  -1,   2 },
            new int[] {  1,  -2,   2,  -2,   2 },
            new int[] {  0,   1,   0,   2,   2 },
            new int[] {  0,   1,  -2,  -1,  -2 },

            /* 421-430 */
            new int[] {  1,   1,  -2,   0,   0 },
            new int[] {  0,   2,   2,  -2,   0 },
            new int[] {  3,  -1,  -2,  -1,  -2 },
            new int[] {  1,   0,   0,  -6,   0 },
            new int[] {  1,   0,  -2,  -4,   0 },
            new int[] {  2,   1,   0,  -4,   1 },
            new int[] {  2,   0,   2,   0,  -1 },
            new int[] {  2,   0,  -4,   0,  -1 },
            new int[] {  0,   0,   3,   0,   2 },
            new int[] {  2,   1,  -2,  -2,  -1 },

            /* 431-440 */
            new int[] {  1,  -2,   0,   0,   1 },
            new int[] {  2,  -1,   0,  -4,   0 },
            new int[] {  0,   0,   0,   3,   0 },
            new int[] {  5,   0,   2,  -2,   2 },
            new int[] {  1,   2,  -2,  -4,  -2 },
            new int[] {  1,   0,   4,  -4,   2 },
            new int[] {  0,   0,   4,  -1,   2 },
            new int[] {  3,   1,   0,  -4,   0 },
            new int[] {  3,   0,   0,  -6,   0 },
            new int[] {  2,   0,   0,   2,   2 },

            /* 441-450 */
            new int[] {  2,  -2,   2,   0,   2 },
            new int[] {  1,   0,   0,  -3,   1 },
            new int[] {  1,  -2,  -2,   0,  -2 },
            new int[] {  1,  -1,  -2,  -3,  -2 },
            new int[] {  0,   0,   2,  -2,  -2 },
            new int[] {  2,   0,  -2,  -4,   0 },
            new int[] {  1,   0,  -4,   0,   0 },
            new int[] {  0,   1,   0,  -1,   0 },
            new int[] {  4,   0,   0,   0,  -1 },
            new int[] {  3,   0,   2,  -1,   2 },

            /* 451-460 */
            new int[] {  3,  -1,   2,   0,   1 },
            new int[] {  2,   0,   2,  -1,   1 },
            new int[] {  1,   2,   2,  -2,   1 },
            new int[] {  1,   1,   0,   2,  -1 },
            new int[] {  0,   2,   2,   0,   1 },
            new int[] {  3,   1,   2,   0,   1 },
            new int[] {  1,   1,   2,   1,   1 },
            new int[] {  1,   1,   0,  -1,   1 },
            new int[] {  1,  -2,   0,  -2,  -1 },
            new int[] {  4,   0,   0,  -4,   0 },

            /* 461-470 */
            new int[] {  2,   1,   0,   2,   0 },
            new int[] {  1,  -1,   0,   4,   0 },
            new int[] {  0,   1,   0,  -2,   2 },
            new int[] {  0,   0,   2,   0,  -2 },
            new int[] {  1,   0,  -1,   0,   1 },
            new int[] {  3,   0,   2,  -2,   0 },
            new int[] {  2,   0,   2,   2,   0 },
            new int[] {  1,   2,   0,  -4,   0 },
            new int[] {  1,  -1,   0,  -3,   0 },
            new int[] {  0,   1,   0,   4,   0 },

            /* 471 - 480 */
            new int[] {  0,   1,  -2,   0,   0 },
            new int[] {  2,   2,   2,  -2,   2 },
            new int[] {  0,   0,   0,   1,  -2 },
            new int[] {  0,   2,  -2,   0,  -1 },
            new int[] {  4,   0,   2,  -4,   2 },
            new int[] {  2,   0,  -4,   2,  -2 },
            new int[] {  2,  -1,  -2,   0,  -2 },
            new int[] {  1,   1,   4,  -2,   2 },
            new int[] {  1,   1,   2,  -4,   2 },
            new int[] {  1,   0,   2,   3,   2 },

            /* 481-490 */
            new int[] {  1,   0,   0,   4,  -1 },
            new int[] {  0,   0,   0,   4,   2 },
            new int[] {  2,   0,   0,   4,   0 },
            new int[] {  1,   1,  -2,   2,   0 },
            new int[] {  2,   1,   2,   1,   2 },
            new int[] {  2,   1,   2,  -4,   1 },
            new int[] {  2,   0,   2,   1,   1 },
            new int[] {  2,   0,  -4,  -2,  -1 },
            new int[] {  2,   0,  -2,  -6,  -1 },
            new int[] {  2,  -1,   2,  -1,   2 },

            /* 491-500 */
            new int[] {  1,  -2,   2,   0,   1 },
            new int[] {  1,  -2,   0,  -2,   1 },
            new int[] {  1,  -1,   0,  -4,  -1 },
            new int[] {  0,   2,   2,   2,   2 },
            new int[] {  0,   2,  -2,  -4,  -2 },
            new int[] {  0,   1,   2,   3,   2 },
            new int[] {  0,   1,   0,  -4,   1 },
            new int[] {  3,   0,   0,  -2,   1 },
            new int[] {  2,   1,  -2,   0,   1 },
            new int[] {  2,   0,   4,  -2,   1 },

            /* 501-510 */
            new int[] {  2,   0,   0,  -3,  -1 },
            new int[] {  2,  -2,   0,  -2,   1 },
            new int[] {  2,  -1,   2,  -2,   1 },
            new int[] {  1,   0,   0,  -6,  -1 },
            new int[] {  1,  -2,   0,   0,  -1 },
            new int[] {  1,  -2,  -2,  -2,  -1 },
            new int[] {  0,   1,   4,  -2,   1 },
            new int[] {  0,   0,   2,   3,   1 },
            new int[] {  2,  -1,   0,  -1,   0 },
            new int[] {  1,   3,   0,  -2,   0 },

            /* 511-520 */
            new int[] {  0,   3,   0,  -2,   0 },
            new int[] {  2,  -2,   2,  -2,   2 },
            new int[] {  0,   0,   4,  -2,   0 },
            new int[] {  4,  -1,   2,   0,   2 },
            new int[] {  2,   2,  -2,  -4,  -2 },
            new int[] {  4,   1,   2,   0,   2 },
            new int[] {  4,  -1,  -2,  -2,  -2 },
            new int[] {  2,   1,   0,  -2,  -2 },
            new int[] {  2,   1,  -2,  -6,  -2 },
            new int[] {  2,   0,   0,  -1,   1 },

            /* 521-530 */
            new int[] {  2,  -1,  -2,   2,  -1 },
            new int[] {  1,   1,  -2,   2,  -2 },
            new int[] {  1,   1,  -2,  -3,  -2 },
            new int[] {  1,   0,   3,   0,   3 },
            new int[] {  1,   0,  -2,   1,   1 },
            new int[] {  1,   0,  -2,   0,   2 },
            new int[] {  1,  -1,   2,   1,   2 },
            new int[] {  1,  -1,   0,   0,  -2 },
            new int[] {  1,  -1,  -4,   2,  -2 },
            new int[] {  0,   3,  -2,  -2,  -2 },

            /* 531-540 */
            new int[] {  0,   1,   0,   4,   1 },
            new int[] {  0,   0,   4,   2,   2 },
            new int[] {  3,   0,  -2,  -2,   0 },
            new int[] {  2,  -2,   0,   0,   0 },
            new int[] {  1,   1,   2,  -4,   0 },
            new int[] {  1,   1,   0,  -3,   0 },
            new int[] {  1,   0,   2,  -3,   0 },
            new int[] {  1,  -1,   2,  -2,   0 },
            new int[] {  0,   2,   0,   2,   0 },
            new int[] {  0,   0,   2,   4,   0 },

            /* 541-550 */
            new int[] {  1,   0,   1,   0,   0 },
            new int[] {  3,   1,   2,  -2,   1 },
            new int[] {  3,   0,   4,  -2,   2 },
            new int[] {  3,   0,   2,   1,   2 },
            new int[] {  3,   0,   0,   2,  -1 },
            new int[] {  3,   0,   0,   0,   2 },
            new int[] {  3,   0,  -2,   2,  -1 },
            new int[] {  2,   0,   4,  -4,   2 },
            new int[] {  2,   0,   2,  -3,   2 },
            new int[] {  2,   0,   0,   4,   1 },

            /* 551-560 */
            new int[] {  2,   0,   0,  -3,   1 },
            new int[] {  2,   0,  -4,   2,  -1 },
            new int[] {  2,   0,  -2,  -2,   1 },
            new int[] {  2,  -2,   2,   2,   2 },
            new int[] {  2,  -2,   0,  -2,  -2 },
            new int[] {  2,  -1,   0,   2,   1 },
            new int[] {  2,  -1,   0,   2,  -1 },
            new int[] {  1,   1,   2,   4,   2 },
            new int[] {  1,   1,   0,   1,   1 },
            new int[] {  1,   1,   0,   1,  -1 },

            /* 561-570 */
            new int[] {  1,   1,  -2,  -6,  -2 },
            new int[] {  1,   0,   0,  -3,  -1 },
            new int[] {  1,   0,  -4,  -2,  -1 },
            new int[] {  1,   0,  -2,  -6,  -1 },
            new int[] {  1,  -2,   2,   2,   1 },
            new int[] {  1,  -2,  -2,   2,  -1 },
            new int[] {  1,  -1,  -2,  -4,  -1 },
            new int[] {  0,   2,   0,   0,   2 },
            new int[] {  0,   1,   2,  -4,   2 },
            new int[] {  0,   1,  -2,   4,  -1 },

            /* 571-580 */
            new int[] {  5,   0,   0,   0,   0 },
            new int[] {  3,   0,   0,  -3,   0 },
            new int[] {  2,   2,   0,  -4,   0 },
            new int[] {  1,  -1,   2,   2,   0 },
            new int[] {  0,   1,   0,   3,   0 },
            new int[] {  4,   0,  -2,   0,  -1 },
            new int[] {  3,   0,  -2,  -6,  -1 },
            new int[] {  3,   0,  -2,  -1,  -1 },
            new int[] {  2,   1,   2,   2,   1 },
            new int[] {  2,   1,   0,   2,   1 },

            /* 581-590 */
            new int[] {  2,   0,   2,   4,   1 },
            new int[] {  2,   0,   2,  -6,   1 },
            new int[] {  2,   0,   2,  -2,  -1 },
            new int[] {  2,   0,   0,  -6,  -1 },
            new int[] {  2,  -1,  -2,  -2,  -1 },
            new int[] {  1,   2,   2,   0,   1 },
            new int[] {  1,   2,   0,   0,   1 },
            new int[] {  1,   0,   4,   0,   1 },
            new int[] {  1,   0,   2,  -6,   1 },
            new int[] {  1,   0,   2,  -4,  -1 },

            /* 591-600 */
            new int[] {  1,   0,  -1,  -2,  -1 },
            new int[] {  1,  -1,   2,   4,   1 },
            new int[] {  1,  -1,   2,  -3,   1 },
            new int[] {  1,  -1,   0,   4,   1 },
            new int[] {  1,  -1,  -2,   1,  -1 },
            new int[] {  0,   1,   2,  -2,   3 },
            new int[] {  3,   0,   0,  -2,   0 },
            new int[] {  1,   0,   1,  -2,   0 },
            new int[] {  0,   2,   0,  -4,   0 },
            new int[] {  0,   0,   2,  -4,   0 },

            /* 601-610 */
            new int[] {  0,   0,   1,  -1,   0 },
            new int[] {  0,   0,   0,   6,   0 },
            new int[] {  0,   2,   0,   0,  -2 },
            new int[] {  0,   1,  -2,   2,  -3 },
            new int[] {  4,   0,   0,   2,   0 },
            new int[] {  3,   0,   0,  -1,   0 },
            new int[] {  3,  -1,   0,   2,   0 },
            new int[] {  2,   1,   0,   1,   0 },
            new int[] {  2,   1,   0,  -6,   0 },
            new int[] {  2,  -1,   2,   0,   0 },

            /* 611-620 */
            new int[] {  1,   0,   2,  -1,   0 },
            new int[] {  1,  -1,   0,   1,   0 },
            new int[] {  1,  -1,  -2,  -2,   0 },
            new int[] {  0,   1,   2,   2,   0 },
            new int[] {  0,   0,   2,  -3,   0 },
            new int[] {  2,   2,   0,  -2,  -1 },
            new int[] {  2,  -1,  -2,   0,   1 },
            new int[] {  1,   2,   2,  -4,   1 },
            new int[] {  0,   1,   4,  -4,   2 },
            new int[] {  0,   0,   0,   3,   2 },

            /* 621-630 */
            new int[] {  5,   0,   2,   0,   1 },
            new int[] {  4,   1,   2,  -2,   2 },
            new int[] {  4,   0,  -2,  -2,   0 },
            new int[] {  3,   1,   2,   2,   2 },
            new int[] {  3,   1,   0,  -2,   0 },
            new int[] {  3,   1,  -2,  -6,  -2 },
            new int[] {  3,   0,   0,   0,  -2 },
            new int[] {  3,   0,  -2,  -4,  -2 },
            new int[] {  3,  -1,   0,  -3,   0 },
            new int[] {  3,  -1,   0,  -2,   0 },

            /* 631-640 */
            new int[] {  2,   1,   2,   0,   0 },
            new int[] {  2,   1,   2,  -4,   2 },
            new int[] {  2,   1,   2,  -2,   0 },
            new int[] {  2,   1,   0,  -3,   0 },
            new int[] {  2,   1,  -2,   0,  -2 },
            new int[] {  2,   0,   0,  -4,   2 },
            new int[] {  2,   0,   0,  -4,  -2 },
            new int[] {  2,   0,  -2,  -5,  -2 },
            new int[] {  2,  -1,   2,   4,   2 },
            new int[] {  2,  -1,   0,  -2,   2 },

            /* 641-650 */
            new int[] {  1,   3,  -2,  -2,  -2 },
            new int[] {  1,   1,   0,   0,  -2 },
            new int[] {  1,   1,   0,  -6,   0 },
            new int[] {  1,   1,  -2,   1,  -2 },
            new int[] {  1,   1,  -2,  -1,  -2 },
            new int[] {  1,   0,   2,   1,   0 },
            new int[] {  1,   0,   0,   3,   0 },
            new int[] {  1,   0,   0,  -4,   2 },
            new int[] {  1,   0,  -2,   4,  -2 },
            new int[] {  1,  -2,   0,  -1,   0 },

            /* 651-NFLS */
            new int[] {  0,   1,  -4,   2,  -1 },
            new int[] {  1,   0,  -2,   0,  -3 },
            new int[] {  0,   0,   4,  -4,   4 }
        };


        /* Fundamental-argument multipliers:  planetary terms */
        static readonly int[][] Xy06mfapl = {
            /* 1-10 */
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0 },

            /* 11-20 */
            new int[] {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -8,  3,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0 },

            /* 21-30 */
            new int[] {  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2, -1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },

            /* 31-40 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0 },

            /* 41-50 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  0,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2 },
            new int[] {  1,  0,  0,  0,  0,  0,-18, 16,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2 },

            /* 51-60 */
            new int[] {  0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  0,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2 },
            new int[] {  1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1 },
            new int[] {  1,  0, -2,  0, -2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0 },

            /* 61-70 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-11,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0 },

            /* 71-80 */
            new int[] {  0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0 },

            /* 81-90 */
            new int[] {  2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0, -1 },
            new int[] {  2,  0,  0, -2,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0 },
            new int[] {  1,  0,  0, -1,  0,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2 },
            new int[] {  1,  0,  0,  0, -1,  0,-18, 16,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0 },

            /* 91-100 */
            new int[] {  1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2 },
            new int[] {  1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,-16,  4,  5,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0, -2 },

            /* 101-110 */
            new int[] {  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1 },
            new int[] {  2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2, -1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2 },

            /* 111-120 */
            new int[] {  0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0, -6,  8,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0 },

            /* 121-130 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-10,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0, -2 },
            new int[] {  1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },

            /* 131-140 */
            new int[] {  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0 },

            /* 141-150 */
            new int[] {  1,  0,  0, -1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  9,-11,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2 },

            /* 151-160 */
            new int[] {  1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0, -1,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0, -2 },
            new int[] {  0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0 },

            /* 161-170 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0 },

            /* 171-180 */
            new int[] {  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0 },
            new int[] {  1,  0,  0, -1, -1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1 },
            new int[] {  1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  1,  0, -2,  0, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0 },

            /* 181-190 */
            new int[] {  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -8,  3,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1 },

            /* 191-200 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  7,-10,  0,  0,  0,  0,  0, -2 },
            new int[] {  1,  0,  0, -2,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0 },

            /* 201-210 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -4,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2 },
            new int[] {  2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },

            /* 211-220 */
            new int[] {  0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2 },
            new int[] {  2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0 },

            /* 221-230 */
            new int[] {  2,  0,  0, -2,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -9,  4,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0, -2 },

            /* 231-240 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  7,-11,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0, -1,  1,  0,  0,  0 },
            new int[] {  2,  0,  0,  0,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0, -1 },

            /* 241-250 */
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2 },
            new int[] {  0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0 },

            /* 251-260 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0, -2, -2, -2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1 },

            /* 261-270 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2 },
            new int[] {  2,  0,  0, -2, -1,  0,  0, -2,  0,  0,  5,  0,  0,  0 },
            new int[] {  2,  0,  0, -2, -1,  0, -6,  8,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },

            /* 271-280 */
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0, -2,  5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,-15,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2 },
            new int[] {  1,  0,  0, -1,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2 },

            /* 281-290 */
            new int[] {  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0, -2 },

            /* 291-300 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  9,-12,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  1, -1,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0, -1 },

            /* 301-310 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  1, -1, -1,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -5,  0,  0,  0,  0, -2 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  3, -1,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2 },

            /* 311-320 */
            new int[] {  0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1 },
            new int[] {  0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2 },
            new int[] {  0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },

            /* 321-330 */
            new int[] {  0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0 },
            new int[] {  2,  0, -1, -1, -1,  0,  0, -1,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  5,-10,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0 },

            /* 331-340 */
            new int[] {  0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -6,  0,  0,  0,  0,  0, -2 },

            /* 341-350 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-10,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1 },

            /* 351-360 */
            new int[] {  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2, -1,  0,  0, -5,  6,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,  0, -2 },
            new int[] {  2,  0, -1, -1, -1,  0,  0,  3, -7,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0 },

            /* 361-370 */
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  4, -3,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2 },
            new int[] {  1,  0, -2, -2, -2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1 },

            /* 371-380 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -6,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2 },

            /* 381-390 */
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0,  2,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2, -2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  2,  0, -1, -1, -1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
            new int[] {  1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  1,  0, -1,  1, -1,  0,-18, 17,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2, -1,  0, -5,  6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0 },

            /* 391-400 */
            new int[] {  0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-16,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2 },
            new int[] {  0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0 },
            new int[] {  2,  0, -1, -1, -2,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2 },
            new int[] {  2,  0,  0, -2, -1,  0,  0, -2,  0,  4, -5,  0,  0,  0 },

            /* 401-410 */
            new int[] {  2,  0,  0, -2, -1,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  2,  0, -1, -1, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -1, -1,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
            new int[] {  1,  0, -1, -1, -1,  0, 20,-20,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0 },

            /* 411-420 */
            new int[] {  0,  0,  0,  0,  0,  0,  9,-11,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0 },

            /* 421-430 */
            new int[] {  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -6,  0,  0,  0,  0, -2 },
            new int[] {  1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-12,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0 },

            /* 431-440 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -7,  0,  0,  0,  0, -2 },
            new int[] {  1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0 },
            new int[] {  1,  0, -2,  0, -2,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -4,  0,  0,  0,  0,  0, -2 },
            new int[] {  1,  0, -2, -2, -2,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  1,  0, -1,  1, -1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },

            /* 441-450 */
            new int[] {  0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5,-10,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0, -5,  0,  0,  0, -2 },

            /* 451-460 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -5,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1 },
            new int[] {  1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1, -1,  0,  0, -1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0, -2 },

            /* 461-470 */
            new int[] {  1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 },

            /* 471-480 */
            new int[] {  0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  7,-10,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -8,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2 },

            /* 481-490 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0, -1 },
            new int[] {  2,  0,  0, -2, -1,  0,  0, -6,  8,  0,  0,  0,  0,  0 },
            new int[] {  2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0, -1 },

            /* 491-500 */
            new int[] {  0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2, -1,  0,  0, -2,  0,  3, -1,  0,  0,  0 },
            new int[] {  1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2, -1,  0,  0, -2,  0,  2,  0,  0,  0,  0 },

            /* 501-510 */
            new int[] {  1,  0,  0, -1, -1,  0,  0, -3,  4,  0,  0,  0,  0,  0 },
            new int[] {  1,  0, -1,  0, -1,  0, -3,  5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0 },

            /* 511-520 */
            new int[] {  0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  7,-11,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1 },

            /* 521-530 */
            new int[] {  0,  0,  0,  0,  0,  0,  1, -4,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -4,  8, -3,  0,  0,  0,  0 },

            /* 531-540 */
            new int[] {  2,  0,  0, -2,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -1,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  0,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0 },

            /* 541-550 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  2,  2,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0,  0, -4,  8, -3,  0,  0,  0,  0 },

            /* 551-560 */
            new int[] {  1,  0,  0, -2,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0,  1,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0, -4,  5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  2,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0 },

            /* 561-570 */
            new int[] {  0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0 },
            new int[] {  2,  0, -2, -2, -2,  0,  0, -2,  0,  2,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0, -1,  0,-10,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0 },

            /* 571-580 */
            new int[] {  0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0 },
            new int[] {  2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  5, -5,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  1,  5,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  0,  5,  0,  0,  0 },

            /* 581-590 */
            new int[] {  2,  0,  0, -2,  0,  0,  0, -2,  0,  0,  2,  0,  0,  0 },
            new int[] {  2,  0,  0, -2,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  2,  0, -2,  0, -2,  0,  0,  5, -9,  0,  0,  0,  0,  0 },
            new int[] {  2,  0, -1, -1,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0 },
            new int[] {  1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 },
            new int[] {  1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0 },
            new int[] {  1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 },

            /* 591-600 */
            new int[] {  1,  0,  0,  0,  0,  0,  0, -2,  0,  3,  0,  0,  0,  0 },
            new int[] {  1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
            new int[] {  1,  0, -2, -2, -2,  0,  0,  1,  0, -1,  0,  0,  0,  0 },
            new int[] {  1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  2, -2,  0, -1,  0,  2,  0,  0,  0,  0,  0,  0 },

            /* 601-610 */
            new int[] {  0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0,  0, -1,  0,  0, -2,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  1, -1, -1,  0, -5,  7,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 },

            /* 611-620 */
            new int[] {  0,  0,  0,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-14,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  8,-11,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0, -2 },

            /* 621-630 */
            new int[] {  0,  0,  0,  0,  0,  0,  5, -6, -4,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  7,-12,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0, -2 },

            /* 631-640 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  6,-10,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0, -1 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5,-16,  4,  5,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  5,-13,  0,  0,  0,  0, -2 },

            /* 641-650 */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3,  0, -5,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -9,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  3, -7,  0,  0,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  2, -8,  1,  5,  0,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0 },

            /* 651-NFPL */
            new int[] {  0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -6,  3,  0, -2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2 },
            new int[] {  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }
        };

        /* Pointers into amplitudes array, one pointer per frequency */
        static readonly int[] Xy06nc = {

            /* 1-100 */
            1,    21,    37,    51,    65,    79,    91,   103,   115,   127,
            139,   151,   163,   172,   184,   196,   207,   219,   231,   240,
            252,   261,   273,   285,   297,   309,   318,   327,   339,   351,
            363,   372,   384,   396,   405,   415,   423,   435,   444,   452,
            460,   467,   474,   482,   490,   498,   506,   513,   521,   528,
            536,   543,   551,   559,   566,   574,   582,   590,   597,   605,
            613,   620,   628,   636,   644,   651,   658,   666,   674,   680,
            687,   695,   702,   710,   717,   725,   732,   739,   746,   753,
            760,   767,   774,   782,   790,   798,   805,   812,   819,   826,
            833,   840,   846,   853,   860,   867,   874,   881,   888,   895,

            /* 101-200 */
            901,   908,   914,   921,   928,   934,   941,   948,   955,   962,
            969,   976,   982,   989,   996,  1003,  1010,  1017,  1024,  1031,
            1037,  1043,  1050,  1057,  1064,  1071,  1078,  1084,  1091,  1098,
            1104,  1112,  1118,  1124,  1131,  1138,  1145,  1151,  1157,  1164,
            1171,  1178,  1185,  1192,  1199,  1205,  1212,  1218,  1226,  1232,
            1239,  1245,  1252,  1259,  1266,  1272,  1278,  1284,  1292,  1298,
            1304,  1310,  1316,  1323,  1329,  1335,  1341,  1347,  1353,  1359,
            1365,  1371,  1377,  1383,  1389,  1396,  1402,  1408,  1414,  1420,
            1426,  1434,  1440,  1446,  1452,  1459,  1465,  1471,  1477,  1482,
            1488,  1493,  1499,  1504,  1509,  1514,  1520,  1527,  1532,  1538,

            /* 201-300 */
            1543,  1548,  1553,  1558,  1564,  1569,  1574,  1579,  1584,  1589,
            1594,  1596,  1598,  1600,  1602,  1605,  1608,  1610,  1612,  1617,
            1619,  1623,  1625,  1627,  1629,  1632,  1634,  1640,  1642,  1644,
            1646,  1648,  1650,  1652,  1654,  1658,  1660,  1662,  1664,  1668,
            1670,  1672,  1673,  1675,  1679,  1681,  1683,  1684,  1686,  1688,
            1690,  1693,  1695,  1697,  1701,  1703,  1705,  1707,  1709,  1711,
            1712,  1715,  1717,  1721,  1723,  1725,  1727,  1729,  1731,  1733,
            1735,  1737,  1739,  1741,  1743,  1745,  1747,  1749,  1751,  1753,
            1755,  1757,  1759,  1761,  1762,  1764,  1766,  1768,  1769,  1771,
            1773,  1775,  1777,  1779,  1781,  1783,  1785,  1787,  1788,  1790,

            /* 301-400 */
            1792,  1794,  1796,  1798,  1800,  1802,  1804,  1806,  1807,  1809,
            1811,  1815,  1817,  1819,  1821,  1823,  1825,  1827,  1829,  1831,
            1833,  1835,  1837,  1839,  1840,  1842,  1844,  1848,  1850,  1852,
            1854,  1856,  1858,  1859,  1860,  1862,  1864,  1866,  1868,  1869,
            1871,  1873,  1875,  1877,  1879,  1881,  1883,  1885,  1887,  1889,
            1891,  1892,  1896,  1898,  1900,  1901,  1903,  1905,  1907,  1909,
            1910,  1911,  1913,  1915,  1919,  1921,  1923,  1927,  1929,  1931,
            1933,  1935,  1937,  1939,  1943,  1945,  1947,  1948,  1949,  1951,
            1953,  1955,  1957,  1958,  1960,  1962,  1964,  1966,  1968,  1970,
            1971,  1973,  1974,  1975,  1977,  1979,  1980,  1981,  1982,  1984,

            /* 401-500 */
            1986,  1988,  1990,  1992,  1994,  1995,  1997,  1999,  2001,  2003,
            2005,  2007,  2008,  2009,  2011,  2013,  2015,  2017,  2019,  2021,
            2023,  2024,  2025,  2027,  2029,  2031,  2033,  2035,  2037,  2041,
            2043,  2045,  2046,  2047,  2049,  2051,  2053,  2055,  2056,  2057,
            2059,  2061,  2063,  2065,  2067,  2069,  2070,  2071,  2072,  2074,
            2076,  2078,  2080,  2082,  2084,  2086,  2088,  2090,  2092,  2094,
            2095,  2096,  2097,  2099,  2101,  2105,  2106,  2107,  2108,  2109,
            2110,  2111,  2113,  2115,  2119,  2121,  2123,  2125,  2127,  2129,
            2131,  2133,  2135,  2136,  2137,  2139,  2141,  2143,  2145,  2147,
            2149,  2151,  2153,  2155,  2157,  2159,  2161,  2163,  2165,  2167,

            /* 501-600 */
            2169,  2171,  2173,  2175,  2177,  2179,  2181,  2183,  2185,  2186,
            2187,  2188,  2192,  2193,  2195,  2197,  2199,  2201,  2203,  2205,
            2207,  2209,  2211,  2213,  2217,  2219,  2221,  2223,  2225,  2227,
            2229,  2231,  2233,  2234,  2235,  2236,  2237,  2238,  2239,  2240,
            2241,  2244,  2246,  2248,  2250,  2252,  2254,  2256,  2258,  2260,
            2262,  2264,  2266,  2268,  2270,  2272,  2274,  2276,  2278,  2280,
            2282,  2284,  2286,  2288,  2290,  2292,  2294,  2296,  2298,  2300,
            2302,  2303,  2304,  2305,  2306,  2307,  2309,  2311,  2313,  2315,
            2317,  2319,  2321,  2323,  2325,  2327,  2329,  2331,  2333,  2335,
            2337,  2341,  2343,  2345,  2347,  2349,  2351,  2352,  2355,  2356,

            /* 601-700 */
            2357,  2358,  2359,  2361,  2363,  2364,  2365,  2366,  2367,  2368,
            2369,  2370,  2371,  2372,  2373,  2374,  2376,  2378,  2380,  2382,
            2384,  2385,  2386,  2387,  2388,  2389,  2390,  2391,  2392,  2393,
            2394,  2395,  2396,  2397,  2398,  2399,  2400,  2401,  2402,  2403,
            2404,  2405,  2406,  2407,  2408,  2409,  2410,  2411,  2412,  2413,
            2414,  2415,  2417,  2418,  2430,  2438,  2445,  2453,  2460,  2468,
            2474,  2480,  2488,  2496,  2504,  2512,  2520,  2527,  2535,  2543,
            2550,  2558,  2566,  2574,  2580,  2588,  2596,  2604,  2612,  2619,
            2627,  2634,  2642,  2648,  2656,  2664,  2671,  2679,  2685,  2693,
            2701,  2709,  2717,  2725,  2733,  2739,  2747,  2753,  2761,  2769,

            /* 701-800 */
            2777,  2785,  2793,  2801,  2809,  2817,  2825,  2833,  2841,  2848,
            2856,  2864,  2872,  2878,  2884,  2892,  2898,  2906,  2914,  2922,
            2930,  2938,  2944,  2952,  2958,  2966,  2974,  2982,  2988,  2996,
            3001,  3009,  3017,  3025,  3032,  3039,  3045,  3052,  3059,  3067,
            3069,  3076,  3083,  3090,  3098,  3105,  3109,  3111,  3113,  3120,
            3124,  3128,  3132,  3136,  3140,  3144,  3146,  3150,  3158,  3161,
            3165,  3166,  3168,  3172,  3176,  3180,  3182,  3185,  3189,  3193,
            3194,  3197,  3200,  3204,  3208,  3212,  3216,  3219,  3221,  3222,
            3226,  3230,  3234,  3238,  3242,  3243,  3247,  3251,  3254,  3258,
            3262,  3266,  3270,  3274,  3275,  3279,  3283,  3287,  3289,  3293,

            /* 801-900 */
            3296,  3300,  3303,  3307,  3311,  3315,  3319,  3321,  3324,  3327,
            3330,  3334,  3338,  3340,  3342,  3346,  3350,  3354,  3358,  3361,
            3365,  3369,  3373,  3377,  3381,  3385,  3389,  3393,  3394,  3398,
            3402,  3406,  3410,  3413,  3417,  3421,  3425,  3429,  3433,  3435,
            3439,  3443,  3446,  3450,  3453,  3457,  3458,  3461,  3464,  3468,
            3472,  3476,  3478,  3481,  3485,  3489,  3493,  3497,  3501,  3505,
            3507,  3511,  3514,  3517,  3521,  3524,  3525,  3527,  3529,  3533,
            3536,  3540,  3541,  3545,  3548,  3551,  3555,  3559,  3563,  3567,
            3569,  3570,  3574,  3576,  3578,  3582,  3586,  3590,  3593,  3596,
            3600,  3604,  3608,  3612,  3616,  3620,  3623,  3626,  3630,  3632,

            /* 901-1000 */
            3636,  3640,  3643,  3646,  3648,  3652,  3656,  3660,  3664,  3667,
            3669,  3671,  3675,  3679,  3683,  3687,  3689,  3693,  3694,  3695,
            3699,  3703,  3705,  3707,  3710,  3713,  3717,  3721,  3725,  3729,
            3733,  3736,  3740,  3744,  3748,  3752,  3754,  3757,  3759,  3763,
            3767,  3770,  3773,  3777,  3779,  3783,  3786,  3790,  3794,  3798,
            3801,  3805,  3809,  3813,  3817,  3821,  3825,  3827,  3831,  3835,
            3836,  3837,  3840,  3844,  3848,  3852,  3856,  3859,  3863,  3867,
            3869,  3871,  3875,  3879,  3883,  3887,  3890,  3894,  3898,  3901,
            3905,  3909,  3913,  3917,  3921,  3922,  3923,  3924,  3926,  3930,
            3932,  3936,  3938,  3940,  3944,  3948,  3952,  3956,  3959,  3963,

            /* 1001-1100 */
            3965,  3969,  3973,  3977,  3979,  3981,  3982,  3986,  3989,  3993,
            3997,  4001,  4004,  4006,  4009,  4012,  4016,  4020,  4024,  4026,
            4028,  4032,  4036,  4040,  4044,  4046,  4050,  4054,  4058,  4060,
            4062,  4063,  4064,  4068,  4071,  4075,  4077,  4081,  4083,  4087,
            4089,  4091,  4095,  4099,  4101,  4103,  4105,  4107,  4111,  4115,
            4119,  4123,  4127,  4129,  4131,  4135,  4139,  4141,  4143,  4145,
            4149,  4153,  4157,  4161,  4165,  4169,  4173,  4177,  4180,  4183,
            4187,  4191,  4195,  4198,  4201,  4205,  4209,  4212,  4213,  4216,
            4217,  4221,  4223,  4226,  4230,  4234,  4236,  4240,  4244,  4248,
            4252,  4256,  4258,  4262,  4264,  4266,  4268,  4270,  4272,  4276,

            /* 1101-1200 */
            4279,  4283,  4285,  4287,  4289,  4293,  4295,  4299,  4300,  4301,
            4305,  4309,  4313,  4317,  4319,  4323,  4325,  4329,  4331,  4333,
            4335,  4337,  4341,  4345,  4349,  4351,  4353,  4357,  4361,  4365,
            4367,  4369,  4373,  4377,  4381,  4383,  4387,  4389,  4391,  4395,
            4399,  4403,  4407,  4411,  4413,  4414,  4415,  4418,  4419,  4421,
            4423,  4427,  4429,  4431,  4433,  4435,  4437,  4439,  4443,  4446,
            4450,  4452,  4456,  4458,  4460,  4462,  4466,  4469,  4473,  4477,
            4481,  4483,  4487,  4489,  4491,  4493,  4497,  4499,  4501,  4504,
            4506,  4510,  4513,  4514,  4515,  4518,  4521,  4522,  4525,  4526,
            4527,  4530,  4533,  4534,  4537,  4541,  4542,  4543,  4544,  4545,

            /* 1201-1300 */
            4546,  4547,  4550,  4553,  4554,  4555,  4558,  4561,  4564,  4567,
            4568,  4571,  4574,  4575,  4578,  4581,  4582,  4585,  4586,  4588,
            4590,  4592,  4596,  4598,  4602,  4604,  4608,  4612,  4613,  4616,
            4619,  4622,  4623,  4624,  4625,  4626,  4629,  4632,  4633,  4636,
            4639,  4640,  4641,  4642,  4643,  4644,  4645,  4648,  4649,  4650,
            4651,  4652,  4653,  4656,  4657,  4660,  4661,  4664,  4667,  4670,
            4671,  4674,  4675,  4676,  4677,  4678,  4681,  4682,  4683,  4684,
            4687,  4688,  4689,  4692,  4693,  4696,  4697,  4700,  4701,  4702,
            4703,  4704,  4707,  4708,  4711,  4712,  4715,  4716,  4717,  4718,
            4719,  4720,  4721,  4722,  4723,  4726,  4729,  4730,  4733,  4736,

            /* 1301-(NFLS+NFPL) */
            4737,  4740,  4741,  4742,  4745,  4746,  4749,  4752,  4753
        };

        /* Amplitude coefficients (microarcsec);  indexed using the nc array. */
        static readonly double[] Xy06a = {

       /* 1-105 */
             -6844318.44,     9205236.26,1328.67,1538.18,      205833.11,
               153041.79,       -3309.73, 853.32,2037.98,       -2301.27,
           81.46, 120.56, -20.39, -15.22,   1.73,  -1.61,  -0.10,   0.11,
           -0.02,  -0.02,     -523908.04,      573033.42,-544.75,-458.66,
                12814.01,       11714.49, 198.97,-290.91, 155.74,-143.27,
           -2.75,  -1.03,  -1.27,  -1.16,   0.00,  -0.01,      -90552.22,
                97846.69, 111.23, 137.41,2187.91,2024.68,  41.44, -51.26,
           26.92, -24.46,  -0.46,  -0.28,  -0.22,  -0.20,       82168.76,
               -89618.24, -27.64, -29.05,       -2004.36,       -1837.32,
          -36.07,  48.00, -24.43,  22.41,   0.47,   0.24,   0.20,   0.18,
                58707.02,7387.02, 470.05,-192.40, 164.33,       -1312.21,
         -179.73, -28.93, -17.36,  -1.83,  -0.50,   3.57,   0.00,   0.13,
               -20557.78,       22438.42, -20.84, -17.40, 501.82, 459.68,
           59.20, -67.30,   6.08,  -5.61,  -1.36,  -1.19,       28288.28,
         -674.99, -34.69,  35.80, -15.07,-632.54, -11.19,   0.78,  -8.41,
            0.17,   0.01,   0.07,      -15406.85,       20069.50,  15.12,

       /* 106-219 */
           31.80, 448.76, 344.50,  -5.77,   1.41,   4.59,  -5.02,   0.17,
            0.24,      -11991.74,       12902.66,  32.46,  36.70, 288.49,
          268.14,   5.70,  -7.06,   3.57,  -3.23,  -0.06,  -0.04,
                -8584.95,       -9592.72,   4.42, -13.20,-214.50, 192.06,
           23.87,  29.83,   2.54,   2.40,   0.60,  -0.48,5095.50,
                -6918.22,   7.19,   3.92,-154.91,-113.94,   2.86,  -1.04,
           -1.52,   1.73,  -0.07,  -0.10,       -4910.93,       -5331.13,
            0.76,   0.40,-119.21, 109.81,   2.16,   3.20,   1.46,   1.33,
            0.04,  -0.02,       -6245.02,-123.48,  -6.68,  -8.20,  -2.76,
          139.64,   2.71,   0.15,   1.86,2511.85,       -3323.89,   1.07,
           -0.90, -74.33, -56.17,   1.16,  -0.01,  -0.75,   0.83,  -0.02,
           -0.04,2307.58,3143.98,  -7.52,   7.50,  70.31, -51.60,   1.46,
            0.16,  -0.69,  -0.79,   0.02,  -0.05,2372.58,2554.51,   5.93,
           -6.60,  57.12, -53.05,  -0.96,  -1.24,  -0.71,  -0.64,  -0.01,
                -2053.16,2636.13,   5.13,   7.80,  58.94,  45.91,  -0.42,
           -0.12,   0.61,  -0.66,   0.02,   0.03,       -1825.49,

       /* 220-339 */
                -2423.59,   1.23,  -2.00, -54.19,  40.82,  -1.07,  -1.02,
            0.54,   0.61,  -0.04,   0.04,2521.07,-122.28,  -5.97,   2.90,
           -2.73, -56.37,  -0.82,   0.13,  -0.75,       -1534.09,1645.01,
            6.29,   6.80,  36.78,  34.30,   0.92,  -1.25,   0.46,  -0.41,
           -0.02,  -0.01,1898.27,  47.70,  -0.72,   2.50,   1.07, -42.45,
           -0.94,   0.02,  -0.56,       -1292.02,       -1387.00,   0.00,
            0.00, -31.01,  28.89,   0.68,   0.00,   0.38,   0.35,  -0.01,
           -0.01,       -1234.96,1323.81,   5.21,   5.90,  29.60,  27.61,
            0.74,  -1.22,   0.37,  -0.33,  -0.02,  -0.01,1137.48,
                -1233.89,  -0.04,  -0.30, -27.59, -25.43,  -0.61,   1.00,
           -0.34,   0.31,   0.01,   0.01,-813.13,       -1075.60,   0.40,
            0.30, -24.05,  18.18,  -0.40,  -0.01,   0.24,   0.27,  -0.01,
            0.01,1163.22, -60.90,  -2.94,   1.30,  -1.36, -26.01,  -0.58,
            0.07,  -0.35,1029.70, -55.55,  -2.63,   1.10,  -1.25, -23.02,
           -0.52,   0.06,  -0.31,-556.26, 852.85,   3.16,  -4.48,  19.06,
           12.44,  -0.81,  -0.27,   0.17,  -0.21,   0.00,   0.02,-603.52,

       /* 340-467 */
         -800.34,   0.44,   0.10, -17.90,  13.49,  -0.08,  -0.01,   0.18,
            0.20,  -0.01,   0.01,-628.24, 684.99,  -0.64,  -0.50,  15.32,
           14.05,   3.18,  -4.19,   0.19,  -0.17,  -0.09,  -0.07,-866.48,
          -16.26,   0.52,  -1.30,  -0.36,  19.37,   0.43,  -0.01,   0.26,
         -512.37, 695.54,  -1.47,  -1.40,  15.55,  11.46,  -0.16,   0.03,
            0.15,  -0.17,   0.01,   0.01, 506.65, 643.75,   2.54,  -2.62,
           14.40, -11.33,  -0.77,  -0.06,  -0.15,  -0.16,   0.00,   0.01,
          664.57,  16.81,  -0.40,   1.00,   0.38, -14.86,  -3.71,  -0.09,
           -0.20, 405.91, 522.11,   0.99,  -1.50,  11.67,  -9.08,  -0.25,
           -0.02,  -0.12,  -0.13,-305.78, 326.60,   1.75,   1.90,   7.30,
            6.84,   0.20,  -0.04, 300.99,-325.03,  -0.44,  -0.50,  -7.27,
           -6.73,  -1.01,   0.01,   0.00,   0.08,   0.00,   0.02, 438.51,
           10.47,  -0.56,  -0.20,   0.24,  -9.81,  -0.24,   0.01,  -0.13,
         -264.02, 335.24,   0.99,   1.40,   7.49,   5.90,  -0.27,  -0.02,
          284.09, 307.03,   0.32,  -0.40,   6.87,  -6.35,  -0.99,  -0.01,
         -250.54, 327.11,   0.08,   0.40,   7.31,   5.60,  -0.30, 230.72,

       /* 468-595 */
         -304.46,   0.08,  -0.10,  -6.81,  -5.16,   0.27, 229.78, 304.17,
           -0.60,   0.50,   6.80,  -5.14,   0.33,   0.01, 256.30,-276.81,
           -0.28,  -0.40,  -6.19,  -5.73,  -0.14,   0.01,-212.82, 269.45,
            0.84,   1.20,   6.02,   4.76,   0.14,  -0.02, 196.64, 272.05,
           -0.84,   0.90,   6.08,  -4.40,   0.35,   0.02, 188.95, 272.22,
           -0.12,   0.30,   6.09,  -4.22,   0.34,-292.37,  -5.10,  -0.32,
           -0.40,  -0.11,   6.54,   0.14,   0.01, 161.79,-220.67,   0.24,
            0.10,  -4.93,  -3.62,  -0.08, 261.54, -19.94,  -0.95,   0.20,
           -0.45,  -5.85,  -0.13,   0.02, 142.16,-190.79,   0.20,   0.10,
           -4.27,  -3.18,  -0.07, 187.95,  -4.11,  -0.24,   0.30,  -0.09,
           -4.20,  -0.09,   0.01,   0.00,   0.00, -79.08, 167.90,   0.04,
            0.00,   3.75,   1.77, 121.98, 131.04,  -0.08,   0.10,   2.93,
           -2.73,  -0.06,-172.95,  -8.11,  -0.40,  -0.20,  -0.18,   3.87,
            0.09,   0.01,-160.15, -55.30, -14.04,  13.90,  -1.23,   3.58,
            0.40,   0.31,-115.40, 123.20,   0.60,   0.70,   2.75,   2.58,
            0.08,  -0.01,-168.26,  -2.00,   0.20,  -0.20,  -0.04,   3.76,

       /* 596-723 */
            0.08,-114.49, 123.20,   0.32,   0.40,   2.75,   2.56,   0.07,
           -0.01, 112.14, 120.70,   0.28,  -0.30,   2.70,  -2.51,  -0.07,
           -0.01, 161.34,   4.03,   0.20,   0.20,   0.09,  -3.61,  -0.08,
           91.31, 126.64,  -0.40,   0.40,   2.83,  -2.04,  -0.04,   0.01,
          105.29, 112.90,   0.44,  -0.50,   2.52,  -2.35,  -0.07,  -0.01,
           98.69,-106.20,  -0.28,  -0.30,  -2.37,  -2.21,  -0.06,   0.01,
           86.74,-112.94,  -0.08,  -0.20,  -2.53,  -1.94,  -0.05,-134.81,
            3.51,   0.20,  -0.20,   0.08,   3.01,   0.07,  79.03, 107.31,
           -0.24,   0.20,   2.40,  -1.77,  -0.04,   0.01, 132.81, -10.77,
           -0.52,   0.10,  -0.24,  -2.97,  -0.07,   0.01,-130.31,  -0.90,
            0.04,   0.00,   0.00,   2.91, -78.56,  85.32,   0.00,   0.00,
            1.91,   1.76,   0.04,   0.00,   0.00, -41.53,  89.10,   0.02,
            0.00,   1.99,   0.93,  66.03, -71.00,  -0.20,  -0.20,  -1.59,
           -1.48,  -0.04,  60.50,  64.70,   0.36,  -0.40,   1.45,  -1.35,
           -0.04,  -0.01, -52.27, -70.01,   0.00,   0.00,  -1.57,   1.17,
            0.03, -52.95,  66.29,   0.32,   0.40,   1.48,   1.18,   0.04,

       /* 724-851 */
           -0.01,  51.02,  67.25,   0.00,   0.00,   1.50,  -1.14,  -0.03,
          -55.66, -60.92,   0.16,  -0.20,  -1.36,   1.24,   0.03, -54.81,
          -59.20,  -0.08,   0.20,  -1.32,   1.23,   0.03,  51.32, -55.60,
            0.00,   0.00,  -1.24,  -1.15,  -0.03,  48.29,  51.80,   0.20,
           -0.20,   1.16,  -1.08,  -0.03, -45.59, -49.00,  -0.12,   0.10,
           -1.10,   1.02,   0.03,  40.54, -52.69,  -0.04,  -0.10,  -1.18,
           -0.91,  -0.02, -40.58, -49.51,  -1.00,   1.00,  -1.11,   0.91,
            0.04,   0.02, -43.76,  46.50,   0.36,   0.40,   1.04,   0.98,
            0.03,  -0.01,  62.65,  -5.00,  -0.24,   0.00,  -0.11,  -1.40,
           -0.03,   0.01, -38.57,  49.59,   0.08,   0.10,   1.11,   0.86,
            0.02, -33.22, -44.04,   0.08,  -0.10,  -0.98,   0.74,   0.02,
           37.15, -39.90,  -0.12,  -0.10,  -0.89,  -0.83,  -0.02,  36.68,
          -39.50,  -0.04,  -0.10,  -0.88,  -0.82,  -0.02, -53.22,  -3.91,
           -0.20,   0.00,  -0.09,   1.19,   0.03,  32.43, -42.19,  -0.04,
           -0.10,  -0.94,  -0.73,  -0.02, -51.00,  -2.30,  -0.12,  -0.10,
            0.00,   1.14, -29.53, -39.11,   0.04,   0.00,  -0.87,   0.66,

       /* 852-979 */
            0.02,  28.50, -38.92,  -0.08,  -0.10,  -0.87,  -0.64,  -0.02,
           26.54,  36.95,  -0.12,   0.10,   0.83,  -0.59,  -0.01,  26.54,
           34.59,   0.04,  -0.10,   0.77,  -0.59,  -0.02,  28.35, -32.55,
           -0.16,   0.20,  -0.73,  -0.63,  -0.01, -28.00,  30.40,   0.00,
            0.00,   0.68,   0.63,   0.01, -27.61,  29.40,   0.20,   0.20,
            0.66,   0.62,   0.02,  40.33,   0.40,  -0.04,   0.10,   0.00,
           -0.90, -23.28,  31.61,  -0.08,  -0.10,   0.71,   0.52,   0.01,
           37.75,   0.80,   0.04,   0.10,   0.00,  -0.84,  23.66,  25.80,
            0.00,   0.00,   0.58,  -0.53,  -0.01,  21.01, -27.91,   0.00,
            0.00,  -0.62,  -0.47,  -0.01, -34.81,   2.89,   0.04,   0.00,
            0.00,   0.78, -23.49, -25.31,   0.00,   0.00,  -0.57,   0.53,
            0.01, -23.47,  25.20,   0.16,   0.20,   0.56,   0.52,   0.02,
           19.58,  27.50,  -0.12,   0.10,   0.62,  -0.44,  -0.01, -22.67,
          -24.40,  -0.08,   0.10,  -0.55,   0.51,   0.01, -19.97,  25.00,
            0.12,   0.20,   0.56,   0.45,   0.01,  21.28, -22.80,  -0.08,
           -0.10,  -0.51,  -0.48,  -0.01, -30.47,   0.91,   0.04,   0.00,

       /* 980-1107 */
            0.00,   0.68,  18.58,  24.00,   0.04,  -0.10,   0.54,  -0.42,
           -0.01, -18.02,  24.40,  -0.04,  -0.10,   0.55,   0.40,   0.01,
           17.74,  22.50,   0.08,  -0.10,   0.50,  -0.40,  -0.01, -19.41,
           20.70,   0.08,   0.10,   0.46,   0.43,   0.01, -18.64,  20.11,
            0.00,   0.00,   0.45,   0.42,   0.01, -16.75,  21.60,   0.04,
            0.10,   0.48,   0.37,   0.01, -18.42, -20.00,   0.00,   0.00,
           -0.45,   0.41,   0.01, -26.77,   1.41,   0.08,   0.00,   0.00,
            0.60, -26.17,  -0.19,   0.00,   0.00,   0.00,   0.59, -15.52,
           20.51,   0.00,   0.00,   0.46,   0.35,   0.01, -25.42,  -1.91,
           -0.08,   0.00,  -0.04,   0.57,   0.45, -17.42,  18.10,   0.00,
            0.00,   0.40,   0.39,   0.01,  16.39, -17.60,  -0.08,  -0.10,
           -0.39,  -0.37,  -0.01, -14.37,  18.91,   0.00,   0.00,   0.42,
            0.32,   0.01,  23.39,  -2.40,  -0.12,   0.00,   0.00,  -0.52,
           14.32, -18.50,  -0.04,  -0.10,  -0.41,  -0.32,  -0.01,  15.69,
           17.08,   0.00,   0.00,   0.38,  -0.35,  -0.01, -22.99,   0.50,
            0.04,   0.00,   0.00,   0.51,   0.00,   0.00,  14.47, -17.60,

       /* 1108-1235 */
           -0.01,   0.00,  -0.39,  -0.32, -13.33,  18.40,  -0.04,  -0.10,
            0.41,   0.30,  22.47,  -0.60,  -0.04,   0.00,   0.00,  -0.50,
          -12.78, -17.41,   0.04,   0.00,  -0.39,   0.29,   0.01, -14.10,
          -15.31,   0.04,   0.00,  -0.34,   0.32,   0.01,  11.98,  16.21,
           -0.04,   0.00,   0.36,  -0.27,  -0.01,  19.65,  -1.90,  -0.08,
            0.00,   0.00,  -0.44,  19.61,  -1.50,  -0.08,   0.00,   0.00,
           -0.44,  13.41, -14.30,  -0.04,  -0.10,  -0.32,  -0.30,  -0.01,
          -13.29,  14.40,   0.00,   0.00,   0.32,   0.30,   0.01,  11.14,
          -14.40,  -0.04,   0.00,  -0.32,  -0.25,  -0.01,  12.24, -13.38,
            0.04,   0.00,  -0.30,  -0.27,  -0.01,  10.07, -13.81,   0.04,
            0.00,  -0.31,  -0.23,  -0.01,  10.46,  13.10,   0.08,  -0.10,
            0.29,  -0.23,  -0.01,  16.55,  -1.71,  -0.08,   0.00,   0.00,
           -0.37,   9.75, -12.80,   0.00,   0.00,  -0.29,  -0.22,  -0.01,
            9.11,  12.80,   0.00,   0.00,   0.29,  -0.20,   0.00,   0.00,
           -6.44, -13.80,   0.00,   0.00,  -0.31,   0.14,  -9.19, -12.00,
            0.00,   0.00,  -0.27,   0.21, -10.30,  10.90,   0.08,   0.10,

       /* 1236-1363 */
            0.24,   0.23,   0.01,  14.92,  -0.80,  -0.04,   0.00,   0.00,
           -0.33,  10.02, -10.80,   0.00,   0.00,  -0.24,  -0.22,  -0.01,
           -9.75,  10.40,   0.04,   0.00,   0.23,   0.22,   0.01,   9.67,
          -10.40,  -0.04,   0.00,  -0.23,  -0.22,  -0.01,  -8.28, -11.20,
            0.04,   0.00,  -0.25,   0.19,  13.32,  -1.41,  -0.08,   0.00,
            0.00,  -0.30,   8.27,  10.50,   0.04,   0.00,   0.23,  -0.19,
            0.00,   0.00,  13.13,   0.00,   0.00,   0.00,   0.00,  -0.29,
          -12.93,   0.70,   0.04,   0.00,   0.00,   0.29,   7.91, -10.20,
            0.00,   0.00,  -0.23,  -0.18,  -7.84, -10.00,  -0.04,   0.00,
           -0.22,   0.18,   7.44,   9.60,   0.00,   0.00,   0.21,  -0.17,
           -7.64,   9.40,   0.08,   0.10,   0.21,   0.17,   0.01, -11.38,
            0.60,   0.04,   0.00,   0.00,   0.25,  -7.48,   8.30,   0.00,
            0.00,   0.19,   0.17, -10.98,  -0.20,   0.00,   0.00,   0.00,
            0.25,  10.98,   0.20,   0.00,   0.00,   0.00,  -0.25,   7.40,
           -7.90,  -0.04,   0.00,  -0.18,  -0.17,  -6.09,   8.40,  -0.04,
            0.00,   0.19,   0.14,  -6.94,  -7.49,   0.00,   0.00,  -0.17,

       /* 1364-1491 */
            0.16,   6.92,   7.50,   0.04,   0.00,   0.17,  -0.15,   6.20,
            8.09,   0.00,   0.00,   0.18,  -0.14,  -6.12,   7.80,   0.04,
            0.00,   0.17,   0.14,   5.85,  -7.50,   0.00,   0.00,  -0.17,
           -0.13,  -6.48,   6.90,   0.08,   0.10,   0.15,   0.14,   0.01,
            6.32,   6.90,   0.00,   0.00,   0.15,  -0.14,   5.61,  -7.20,
            0.00,   0.00,  -0.16,  -0.13,   9.07,   0.00,   0.00,   0.00,
            0.00,  -0.20,   5.25,   6.90,   0.00,   0.00,   0.15,  -0.12,
           -8.47,  -0.40,   0.00,   0.00,   0.00,   0.19,   6.32,  -5.39,
           -1.11,   1.10,  -0.12,  -0.14,   0.02,   0.02,   5.73,  -6.10,
           -0.04,   0.00,  -0.14,  -0.13,   4.70,   6.60,  -0.04,   0.00,
            0.15,  -0.11,  -4.90,  -6.40,   0.00,   0.00,  -0.14,   0.11,
           -5.33,   5.60,   0.04,   0.10,   0.13,   0.12,   0.01,  -4.81,
            6.00,   0.04,   0.00,   0.13,   0.11,   5.13,   5.50,   0.04,
            0.00,   0.12,  -0.11,   4.50,   5.90,   0.00,   0.00,   0.13,
           -0.10,  -4.22,   6.10,   0.00,   0.00,   0.14,  -4.53,   5.70,
            0.00,   0.00,   0.13,   0.10,   4.18,   5.70,   0.00,   0.00,

       /* 1492-1619 */
            0.13,  -4.75,  -5.19,   0.00,   0.00,  -0.12,   0.11,  -4.06,
            5.60,   0.00,   0.00,   0.13,  -3.98,   5.60,  -0.04,   0.00,
            0.13,   4.02,  -5.40,   0.00,   0.00,  -0.12,   4.49,  -4.90,
           -0.04,   0.00,  -0.11,  -0.10,  -3.62,  -5.40,  -0.16,   0.20,
           -0.12,   0.00,   0.01,   4.38,   4.80,   0.00,   0.00,   0.11,
           -6.40,  -0.10,   0.00,   0.00,   0.00,   0.14,  -3.98,   5.00,
            0.04,   0.00,   0.11,  -3.82,  -5.00,   0.00,   0.00,  -0.11,
           -3.71,   5.07,   0.00,   0.00,   0.11,   4.14,   4.40,   0.00,
            0.00,   0.10,  -6.01,  -0.50,  -0.04,   0.00,   0.00,   0.13,
           -4.04,   4.39,   0.00,   0.00,   0.10,   3.45,  -4.72,   0.00,
            0.00,  -0.11,   3.31,   4.71,   0.00,   0.00,   0.11,   3.26,
           -4.50,   0.00,   0.00,  -0.10,  -3.26,  -4.50,   0.00,   0.00,
           -0.10,  -3.34,  -4.40,   0.00,   0.00,  -0.10,  -3.74,  -4.00,
            3.70,   4.00,   3.34,  -4.30,   3.30,  -4.30,  -3.66,   3.90,
            0.04,   3.66,   3.90,   0.04,  -3.62,  -3.90,  -3.61,   3.90,
           -0.20,   5.30,   0.00,   0.00,   0.12,   3.06,   4.30,   3.30,

       /* 1620-1747 */
            4.00,   0.40,   0.20,   3.10,   4.10,  -3.06,   3.90,  -3.30,
           -3.60,  -3.30,   3.36,   0.01,   3.14,   3.40,  -4.57,  -0.20,
            0.00,   0.00,   0.00,   0.10,  -2.70,  -3.60,   2.94,  -3.20,
           -2.90,   3.20,   2.47,  -3.40,   2.55,  -3.30,   2.80,  -3.08,
            2.51,   3.30,  -4.10,   0.30,  -0.12,  -0.10,   4.10,   0.20,
           -2.74,   3.00,   2.46,   3.23,  -3.66,   1.20,  -0.20,   0.20,
            3.74,  -0.40,  -2.51,  -2.80,  -3.74,   2.27,  -2.90,   0.00,
            0.00,  -2.50,   2.70,  -2.51,   2.60,  -3.50,   0.20,   3.38,
           -2.22,  -2.50,   3.26,  -0.40,   1.95,  -2.60,   3.22,  -0.40,
           -0.04,  -1.79,  -2.60,   1.91,   2.50,   0.74,   3.05,  -0.04,
            0.08,   2.11,  -2.30,  -2.11,   2.20,  -1.87,  -2.40,   2.03,
           -2.20,  -2.03,   2.20,   2.98,   0.00,   0.00,   2.98,  -1.71,
            2.40,   2.94,  -0.10,  -0.12,   0.10,   1.67,   2.40,  -1.79,
            2.30,  -1.79,   2.20,  -1.67,   2.20,   1.79,  -2.00,   1.87,
           -1.90,   1.63,  -2.10,  -1.59,   2.10,   1.55,  -2.10,  -1.55,
            2.10,  -2.59,  -0.20,  -1.75,  -1.90,  -1.75,   1.90,  -1.83,

       /* 1748-1875 */
           -1.80,   1.51,   2.00,  -1.51,  -2.00,   1.71,   1.80,   1.31,
            2.10,  -1.43,   2.00,   1.43,   2.00,  -2.43,  -1.51,   1.90,
           -1.47,   1.90,   2.39,   0.20,  -2.39,   1.39,   1.90,   1.39,
           -1.80,   1.47,  -1.60,   1.47,  -1.60,   1.43,  -1.50,  -1.31,
            1.60,   1.27,  -1.60,  -1.27,   1.60,   1.27,  -1.60,   2.03,
            1.35,   1.50,  -1.39,  -1.40,   1.95,  -0.20,  -1.27,   1.49,
            1.19,   1.50,   1.27,   1.40,   1.15,   1.50,   1.87,  -0.10,
           -1.12,  -1.50,   1.87,  -1.11,  -1.50,  -1.11,  -1.50,   0.00,
            0.00,   1.19,   1.40,   1.27,  -1.30,  -1.27,  -1.30,  -1.15,
            1.40,  -1.23,   1.30,  -1.23,  -1.30,   1.22,  -1.29,   1.07,
           -1.40,   1.75,  -0.20,  -1.03,  -1.40,  -1.07,   1.20,  -1.03,
            1.15,   1.07,   1.10,   1.51,  -1.03,   1.10,   1.03,  -1.10,
            0.00,   0.00,  -1.03,  -1.10,   0.91,  -1.20,  -0.88,  -1.20,
           -0.88,   1.20,  -0.95,   1.10,  -0.95,  -1.10,   1.43,  -1.39,
            0.95,  -1.00,  -0.95,   1.00,  -0.80,   1.10,   0.91,  -1.00,
           -1.35,   0.88,   1.00,  -0.83,   1.00,  -0.91,   0.90,   0.91,

       /* 1876-2003 */
            0.90,   0.88,  -0.90,  -0.76,  -1.00,  -0.76,   1.00,   0.76,
            1.00,  -0.72,   1.00,   0.84,  -0.90,   0.84,   0.90,   1.23,
            0.00,   0.00,  -0.52,  -1.10,  -0.68,   1.00,   1.19,  -0.20,
            1.19,   0.76,   0.90,   1.15,  -0.10,   1.15,  -0.10,   0.72,
           -0.90,  -1.15,  -1.15,   0.68,   0.90,  -0.68,   0.90,  -1.11,
            0.00,   0.00,   0.20,   0.79,   0.80,  -1.11,  -0.10,   0.00,
            0.00,  -0.48,  -1.00,  -0.76,  -0.80,  -0.72,  -0.80,  -1.07,
           -0.10,   0.64,   0.80,  -0.64,  -0.80,   0.64,   0.80,   0.40,
            0.60,   0.52,  -0.50,  -0.60,  -0.80,  -0.71,   0.70,  -0.99,
            0.99,   0.56,   0.80,  -0.56,   0.80,   0.68,  -0.70,   0.68,
            0.70,  -0.95,  -0.64,   0.70,   0.64,   0.70,  -0.60,   0.70,
           -0.60,  -0.70,  -0.91,  -0.10,  -0.51,   0.76,  -0.91,  -0.56,
            0.70,   0.88,   0.88,  -0.63,  -0.60,   0.55,  -0.60,  -0.80,
            0.80,  -0.80,  -0.52,   0.60,   0.52,   0.60,   0.52,  -0.60,
           -0.48,   0.60,   0.48,   0.60,   0.48,   0.60,  -0.76,   0.44,
           -0.60,   0.52,  -0.50,  -0.52,   0.50,   0.40,   0.60,  -0.40,

       /* 2004-2131 */
           -0.60,   0.40,  -0.60,   0.72,  -0.72,  -0.51,  -0.50,  -0.48,
            0.50,   0.48,  -0.50,  -0.48,   0.50,  -0.48,   0.50,   0.48,
           -0.50,  -0.48,  -0.50,  -0.68,  -0.68,   0.44,   0.50,  -0.64,
           -0.10,  -0.64,  -0.10,  -0.40,   0.50,   0.40,   0.50,   0.40,
            0.50,   0.00,   0.00,  -0.40,  -0.50,  -0.36,  -0.50,   0.36,
           -0.50,   0.60,  -0.60,   0.40,  -0.40,   0.40,   0.40,  -0.40,
            0.40,  -0.40,   0.40,  -0.56,  -0.56,   0.36,  -0.40,  -0.36,
            0.40,   0.36,  -0.40,  -0.36,  -0.40,   0.36,   0.40,   0.36,
            0.40,  -0.52,   0.52,   0.52,   0.32,   0.40,  -0.32,   0.40,
           -0.32,   0.40,  -0.32,   0.40,   0.32,  -0.40,  -0.32,  -0.40,
            0.32,  -0.40,   0.28,  -0.40,  -0.28,   0.40,   0.28,  -0.40,
            0.28,   0.40,   0.48,  -0.48,   0.48,   0.36,  -0.30,  -0.36,
           -0.30,   0.00,   0.00,   0.20,   0.40,  -0.44,   0.44,  -0.44,
           -0.44,  -0.44,  -0.44,   0.32,  -0.30,   0.32,   0.30,   0.24,
            0.30,  -0.12,  -0.10,  -0.28,   0.30,   0.28,   0.30,   0.28,
            0.30,   0.28,  -0.30,   0.28,  -0.30,   0.28,  -0.30,   0.28,

       /* 2132-2259 */
            0.30,  -0.28,   0.30,   0.40,   0.40,  -0.24,   0.30,   0.24,
           -0.30,   0.24,  -0.30,  -0.24,  -0.30,   0.24,   0.30,   0.24,
           -0.30,  -0.24,   0.30,   0.24,  -0.30,  -0.24,  -0.30,   0.24,
           -0.30,   0.24,   0.30,  -0.24,   0.30,  -0.24,   0.30,   0.20,
           -0.30,   0.20,  -0.30,   0.20,  -0.30,   0.20,   0.30,   0.20,
           -0.30,   0.20,  -0.30,   0.20,   0.30,   0.20,   0.30,  -0.20,
           -0.30,   0.20,  -0.30,   0.20,  -0.30,  -0.36,  -0.36,  -0.36,
           -0.04,   0.30,   0.12,  -0.10,  -0.32,  -0.24,   0.20,   0.24,
            0.20,   0.20,  -0.20,  -0.20,  -0.20,  -0.20,  -0.20,   0.20,
            0.20,   0.20,  -0.20,   0.20,   0.20,   0.20,   0.20,  -0.20,
           -0.20,   0.00,   0.00,  -0.20,  -0.20,  -0.20,   0.20,  -0.20,
            0.20,   0.20,  -0.20,  -0.20,  -0.20,   0.20,   0.20,   0.20,
            0.20,   0.20,  -0.20,   0.20,  -0.20,   0.28,   0.28,   0.28,
            0.28,   0.28,   0.28,  -0.28,   0.28,   0.12,   0.00,   0.24,
            0.16,  -0.20,   0.16,  -0.20,   0.16,  -0.20,   0.16,   0.20,
           -0.16,   0.20,   0.16,   0.20,  -0.16,   0.20,  -0.16,   0.20,

       /* 2260-2387 */
           -0.16,   0.20,   0.16,  -0.20,   0.16,   0.20,   0.16,  -0.20,
           -0.16,   0.20,  -0.16,  -0.20,  -0.16,   0.20,   0.16,   0.20,
            0.16,  -0.20,   0.16,  -0.20,   0.16,   0.20,   0.16,   0.20,
            0.16,   0.20,  -0.16,  -0.20,   0.16,   0.20,  -0.16,   0.20,
            0.16,   0.20,  -0.16,  -0.20,   0.16,  -0.20,   0.16,  -0.20,
           -0.16,  -0.20,   0.24,  -0.24,  -0.24,   0.24,   0.24,   0.12,
            0.20,   0.12,   0.20,  -0.12,  -0.20,   0.12,  -0.20,   0.12,
           -0.20,  -0.12,   0.20,  -0.12,   0.20,  -0.12,  -0.20,   0.12,
            0.20,   0.12,   0.20,   0.12,  -0.20,  -0.12,   0.20,   0.12,
           -0.20,  -0.12,   0.20,   0.12,   0.20,   0.00,   0.00,  -0.12,
            0.20,  -0.12,   0.20,   0.12,  -0.20,  -0.12,   0.20,   0.12,
            0.20,   0.00,  -0.21,  -0.20,   0.00,   0.00,   0.20,  -0.20,
           -0.20,  -0.20,   0.20,  -0.16,  -0.10,   0.00,   0.17,   0.16,
            0.16,   0.16,   0.16,  -0.16,   0.16,   0.16,  -0.16,   0.16,
           -0.16,   0.16,   0.12,   0.10,   0.12,  -0.10,  -0.12,   0.10,
           -0.12,   0.10,   0.12,  -0.10,  -0.12,   0.12,  -0.12,   0.12,

       /* 2388-2515 */
           -0.12,   0.12,  -0.12,  -0.12,  -0.12,  -0.12,  -0.12,  -0.12,
           -0.12,   0.12,   0.12,   0.12,   0.12,  -0.12,  -0.12,   0.12,
            0.12,   0.12,  -0.12,   0.12,  -0.12,  -0.12,  -0.12,   0.12,
           -0.12,  -0.12,   0.12,   0.00,   0.11,   0.11,-122.67, 164.70,
          203.78, 273.50,   3.58,   2.74,   6.18,  -4.56,   0.00,  -0.04,
            0.00,  -0.07,  57.44, -77.10,  95.82, 128.60,  -1.77,  -1.28,
            2.85,  -2.14,  82.14,  89.50,   0.00,   0.00,   2.00,  -1.84,
           -0.04,  47.73, -64.10,  23.79,  31.90,  -1.45,  -1.07,   0.69,
           -0.53, -46.38,  50.50,   0.00,   0.00,   1.13,   1.04,   0.02,
          -18.38,   0.00,  63.80,   0.00,   0.00,   0.41,   0.00,  -1.43,
           59.07,   0.00,   0.00,   0.00,   0.00,  -1.32,  57.28,   0.00,
            0.00,   0.00,   0.00,  -1.28, -48.65,   0.00,  -1.15,   0.00,
            0.00,   1.09,   0.00,   0.03, -18.30,  24.60, -17.30, -23.20,
            0.56,   0.41,  -0.51,   0.39, -16.91,  26.90,   8.43,  13.30,
            0.60,   0.38,   0.31,  -0.19,   1.23,  -1.70, -19.13, -25.70,
           -0.03,  -0.03,  -0.58,   0.43,  -0.72,   0.90, -17.34, -23.30,

       /* 2516-2643 */
            0.03,   0.02,  -0.52,   0.39, -19.49, -21.30,   0.00,   0.00,
           -0.48,   0.44,   0.01,  20.57, -20.10,   0.64,   0.70,  -0.45,
           -0.46,   0.00,  -0.01,   4.89,   5.90, -16.55,  19.90,   0.14,
           -0.11,   0.44,   0.37,  18.22,  19.80,   0.00,   0.00,   0.44,
           -0.41,  -0.01,   4.89,  -5.30, -16.51, -18.00,  -0.11,  -0.11,
           -0.41,   0.37, -17.86,   0.00,  17.10,   0.00,   0.00,   0.40,
            0.00,  -0.38,   0.32,   0.00,  24.42,   0.00,   0.00,  -0.01,
            0.00,  -0.55, -23.79,   0.00,   0.00,   0.00,   0.00,   0.53,
           14.72, -16.00,  -0.32,   0.00,  -0.36,  -0.33,  -0.01,   0.01,
            3.34,  -4.50,  11.86,  15.90,  -0.11,  -0.07,   0.35,  -0.27,
           -3.26,   4.40,  11.62,  15.60,   0.09,   0.07,   0.35,  -0.26,
          -19.53,   0.00,   5.09,   0.00,   0.00,   0.44,   0.00,  -0.11,
          -13.48,  14.70,   0.00,   0.00,   0.33,   0.30,   0.01,  10.86,
          -14.60,   3.18,   4.30,  -0.33,  -0.24,   0.09,  -0.07, -11.30,
          -15.10,   0.00,   0.00,  -0.34,   0.25,   0.01,   2.03,  -2.70,
           10.82,  14.50,  -0.07,  -0.05,   0.32,  -0.24,  17.46,   0.00,

       /* 2644-2771 */
            0.00,   0.00,   0.00,  -0.39,  16.43,   0.00,   0.52,   0.00,
            0.00,  -0.37,   0.00,  -0.01,   9.35,   0.00,  13.29,   0.00,
            0.00,  -0.21,   0.00,  -0.30, -10.42,  11.40,   0.00,   0.00,
            0.25,   0.23,   0.01,   0.44,   0.50, -10.38,  11.30,   0.02,
           -0.01,   0.25,   0.23, -14.64,   0.00,   0.00,   0.00,   0.00,
            0.33,   0.56,   0.80,  -8.67,  11.70,   0.02,  -0.01,   0.26,
            0.19,  13.88,   0.00,  -2.47,   0.00,   0.00,  -0.31,   0.00,
            0.06,  -1.99,   2.70,   7.72,  10.30,   0.06,   0.04,   0.23,
           -0.17,  -0.20,   0.00,  13.05,   0.00,   0.00,   0.00,   0.00,
           -0.29,   6.92,  -9.30,   3.34,   4.50,  -0.21,  -0.15,   0.10,
           -0.07,  -6.60,   0.00,  10.70,   0.00,   0.00,   0.15,   0.00,
           -0.24,  -8.04,  -8.70,   0.00,   0.00,  -0.19,   0.18, -10.58,
            0.00,  -3.10,   0.00,   0.00,   0.24,   0.00,   0.07,  -7.32,
            8.00,  -0.12,  -0.10,   0.18,   0.16,   1.63,   1.70,   6.96,
           -7.60,   0.03,  -0.04,  -0.17,  -0.16,  -3.62,   0.00,   9.86,
            0.00,   0.00,   0.08,   0.00,  -0.22,   0.20,  -0.20,  -6.88,

       /* 2772-2899 */
           -7.50,   0.00,   0.00,  -0.17,   0.15,  -8.99,   0.00,   4.02,
            0.00,   0.00,   0.20,   0.00,  -0.09,  -1.07,   1.40,  -5.69,
           -7.70,   0.03,   0.02,  -0.17,   0.13,   6.48,  -7.20,  -0.48,
           -0.50,  -0.16,  -0.14,  -0.01,   0.01,   5.57,  -7.50,   1.07,
            1.40,  -0.17,  -0.12,   0.03,  -0.02,   8.71,   0.00,   3.54,
            0.00,   0.00,  -0.19,   0.00,  -0.08,   0.40,   0.00,   9.27,
            0.00,   0.00,  -0.01,   0.00,  -0.21,  -6.13,   6.70,  -1.19,
           -1.30,   0.15,   0.14,  -0.03,   0.03,   5.21,  -5.70,  -2.51,
           -2.60,  -0.13,  -0.12,  -0.06,   0.06,   5.69,  -6.20,  -0.12,
           -0.10,  -0.14,  -0.13,  -0.01,   2.03,  -2.70,   4.53,   6.10,
           -0.06,  -0.05,   0.14,  -0.10,   5.01,   5.50,  -2.51,   2.70,
            0.12,  -0.11,   0.06,   0.06,  -1.91,   2.60,  -4.38,  -5.90,
            0.06,   0.04,  -0.13,   0.10,   4.65,  -6.30,   0.00,   0.00,
           -0.14,  -0.10,  -5.29,   5.70,   0.00,   0.00,   0.13,   0.12,
           -2.23,  -4.00,  -4.65,   4.20,  -0.09,   0.05,   0.10,   0.10,
           -4.53,   6.10,   0.00,   0.00,   0.14,   0.10,   2.47,   2.70,

       /* 2900-3027 */
           -4.46,   4.90,   0.06,  -0.06,   0.11,   0.10,  -5.05,   5.50,
            0.84,   0.90,   0.12,   0.11,   0.02,  -0.02,   4.97,  -5.40,
           -1.71,   0.00,  -0.12,  -0.11,   0.00,   0.04,  -0.99,  -1.30,
            4.22,  -5.70,  -0.03,   0.02,  -0.13,  -0.09,   0.99,   1.40,
            4.22,  -5.60,   0.03,  -0.02,  -0.13,  -0.09,  -4.69,  -5.20,
            0.00,   0.00,  -0.12,   0.10,  -3.42,   0.00,   6.09,   0.00,
            0.00,   0.08,   0.00,  -0.14,  -4.65,  -5.10,   0.00,   0.00,
           -0.11,   0.10,   0.00,   0.00,  -4.53,  -5.00,   0.00,   0.00,
           -0.11,   0.10,  -2.43,  -2.70,  -3.82,   4.20,  -0.06,   0.05,
            0.10,   0.09,   0.00,   0.00,  -4.53,   4.90,   0.00,   0.00,
            0.11,   0.10,  -4.49,  -4.90,   0.00,   0.00,  -0.11,   0.10,
            2.67,  -2.90,  -3.62,  -3.90,  -0.06,  -0.06,  -0.09,   0.08,
            3.94,  -5.30,   0.00,   0.00,  -0.12,  -3.38,   3.70,  -2.78,
           -3.10,   0.08,   0.08,  -0.07,   0.06,   3.18,  -3.50,  -2.82,
           -3.10,  -0.08,  -0.07,  -0.07,   0.06,  -5.77,   0.00,   1.87,
            0.00,   0.00,   0.13,   0.00,  -0.04,   3.54,  -4.80,  -0.64,

       /* 3028-3155 */
           -0.90,  -0.11,   0.00,  -0.02,  -3.50,  -4.70,   0.68,  -0.90,
           -0.11,   0.00,  -0.02,   5.49,   0.00,   0.00,   0.00,   0.00,
           -0.12,   1.83,  -2.50,   2.63,   3.50,  -0.06,   0.00,   0.08,
            3.02,  -4.10,   0.68,   0.90,  -0.09,   0.00,   0.02,   0.00,
            0.00,   5.21,   0.00,   0.00,   0.00,   0.00,  -0.12,  -3.54,
            3.80,   2.70,   3.60,  -1.35,   1.80,   0.08,   0.00,   0.04,
           -2.90,   3.90,   0.68,   0.90,   0.09,   0.00,   0.02,   0.80,
           -1.10,  -2.78,  -3.70,  -0.02,   0.00,  -0.08,   4.10,   0.00,
           -2.39,   0.00,   0.00,  -0.09,   0.00,   0.05,  -1.59,   2.10,
            2.27,   3.00,   0.05,   0.00,   0.07,  -2.63,   3.50,  -0.48,
           -0.60,  -2.94,  -3.20,  -2.94,   3.20,   2.27,  -3.00,  -1.11,
           -1.50,  -0.07,   0.00,  -0.03,  -0.56,  -0.80,  -2.35,   3.10,
            0.00,  -0.60,  -3.42,   1.90,  -0.12,  -0.10,   2.63,  -2.90,
            2.51,   2.80,  -0.64,   0.70,  -0.48,  -0.60,   2.19,  -2.90,
            0.24,  -0.30,   2.15,   2.90,   2.15,  -2.90,   0.52,   0.70,
            2.07,  -2.80,  -3.10,   0.00,   1.79,   0.00,   0.00,   0.07,

       /* 3156-3283 */
            0.00,  -0.04,   0.88,   0.00,  -3.46,   2.11,   2.80,  -0.36,
            0.50,   3.54,  -0.20,  -3.50,  -1.39,   1.50,  -1.91,  -2.10,
           -1.47,   2.00,   1.39,   1.90,   2.07,  -2.30,   0.91,   1.00,
            1.99,  -2.70,   3.30,   0.00,   0.60,  -0.44,  -0.70,  -1.95,
            2.60,   2.15,  -2.40,  -0.60,  -0.70,   3.30,   0.84,   0.00,
           -3.10,  -3.10,   0.00,  -0.72,  -0.32,   0.40,  -1.87,  -2.50,
            1.87,  -2.50,   0.32,   0.40,  -0.24,   0.30,  -1.87,  -2.50,
           -0.24,  -0.30,   1.87,  -2.50,  -2.70,   0.00,   1.55,   2.03,
            2.20,  -2.98,  -1.99,  -2.20,   0.12,  -0.10,  -0.40,   0.50,
            1.59,   2.10,   0.00,   0.00,  -1.79,   2.00,  -1.03,   1.40,
           -1.15,  -1.60,   0.32,   0.50,   1.39,  -1.90,   2.35,  -1.27,
            1.70,   0.60,   0.80,  -0.32,  -0.40,   1.35,  -1.80,   0.44,
            0.00,   2.23,  -0.84,   0.90,  -1.27,  -1.40,  -1.47,   1.60,
           -0.28,  -0.30,  -0.28,   0.40,  -1.27,  -1.70,   0.28,  -0.40,
           -1.43,  -1.50,   0.00,   0.00,  -1.27,  -1.70,   2.11,  -0.32,
           -0.40,  -1.23,   1.60,   1.19,  -1.30,  -0.72,  -0.80,   0.72,

       /* 3284-3411 */
           -0.80,  -1.15,  -1.30,  -1.35,  -1.50,  -1.19,  -1.60,  -0.12,
            0.20,   1.79,   0.00,  -0.88,  -0.28,   0.40,   1.11,   1.50,
           -1.83,   0.00,   0.56,  -0.12,   0.10,  -1.27,  -1.40,   0.00,
            0.00,   1.15,   1.50,  -0.12,   0.20,   1.11,   1.50,   0.36,
           -0.50,  -1.07,  -1.40,  -1.11,   1.50,   1.67,   0.00,   0.80,
           -1.11,   0.00,   1.43,   1.23,  -1.30,  -0.24,  -1.19,  -1.30,
           -0.24,   0.20,  -0.44,  -0.90,  -0.95,   1.10,   1.07,  -1.40,
            1.15,  -1.30,   1.03,  -1.10,  -0.56,  -0.60,  -0.68,   0.90,
           -0.76,  -1.00,  -0.24,  -0.30,   0.95,  -1.30,   0.56,   0.70,
            0.84,  -1.10,  -0.56,   0.00,  -1.55,   0.91,  -1.30,   0.28,
            0.30,   0.16,  -0.20,   0.95,   1.30,   0.40,  -0.50,  -0.88,
           -1.20,   0.95,  -1.10,  -0.48,  -0.50,   0.00,   0.00,  -1.07,
            1.20,   0.44,  -0.50,   0.95,   1.10,   0.00,   0.00,   0.92,
           -1.30,   0.95,   1.00,  -0.52,   0.60,   1.59,   0.24,  -0.40,
            0.91,   1.20,   0.84,  -1.10,  -0.44,  -0.60,   0.84,   1.10,
           -0.44,   0.60,  -0.44,   0.60,  -0.84,  -1.10,  -0.80,   0.00,

       /* 3412-3539 */
            1.35,   0.76,   0.20,  -0.91,  -1.00,   0.20,  -0.30,  -0.91,
           -1.20,  -0.95,   1.00,  -0.48,  -0.50,   0.88,   1.00,   0.48,
           -0.50,  -0.95,  -1.10,   0.20,  -0.20,  -0.99,   1.10,  -0.84,
            1.10,  -0.24,  -0.30,   0.20,  -0.30,   0.84,   1.10,  -1.39,
            0.00,  -0.28,  -0.16,   0.20,   0.84,   1.10,   0.00,   0.00,
            1.39,   0.00,   0.00,  -0.95,   1.00,   1.35,  -0.99,   0.00,
            0.88,  -0.52,   0.00,  -1.19,   0.20,   0.20,   0.76,  -1.00,
            0.00,   0.00,   0.76,   1.00,   0.00,   0.00,   0.76,   1.00,
           -0.76,   1.00,   0.00,   0.00,   1.23,   0.76,   0.80,  -0.32,
            0.40,  -0.72,   0.80,  -0.40,  -0.40,   0.00,   0.00,  -0.80,
           -0.90,  -0.68,   0.90,  -0.16,  -0.20,  -0.16,  -0.20,   0.68,
           -0.90,  -0.36,   0.50,  -0.56,  -0.80,   0.72,  -0.90,   0.44,
           -0.60,  -0.48,  -0.70,  -0.16,   0.00,  -1.11,   0.32,   0.00,
           -1.07,   0.60,  -0.80,  -0.28,  -0.40,  -0.64,   0.00,   0.91,
            1.11,   0.64,  -0.90,   0.76,  -0.80,   0.00,   0.00,  -0.76,
           -0.80,   1.03,   0.00,  -0.36,  -0.64,  -0.70,   0.36,  -0.40,

       /* 3540-3667 */
            1.07,   0.36,  -0.50,  -0.52,  -0.70,   0.60,   0.00,   0.88,
            0.95,   0.00,   0.48,   0.16,  -0.20,   0.60,   0.80,   0.16,
           -0.20,  -0.60,  -0.80,   0.00,  -1.00,   0.12,   0.20,   0.16,
           -0.20,   0.68,   0.70,   0.59,  -0.80,  -0.99,  -0.56,  -0.60,
            0.36,  -0.40,  -0.68,  -0.70,  -0.68,  -0.70,  -0.36,  -0.50,
           -0.44,   0.60,   0.64,   0.70,  -0.12,   0.10,  -0.52,   0.60,
            0.36,   0.40,   0.00,   0.00,   0.95,  -0.84,   0.00,   0.44,
            0.56,   0.60,   0.32,  -0.30,   0.00,   0.00,   0.60,   0.70,
            0.00,   0.00,   0.60,   0.70,  -0.12,  -0.20,   0.52,  -0.70,
            0.00,   0.00,   0.56,   0.70,  -0.12,   0.10,  -0.52,  -0.70,
            0.00,   0.00,   0.88,  -0.76,   0.00,  -0.44,   0.00,   0.00,
           -0.52,  -0.70,   0.52,  -0.70,   0.36,  -0.40,  -0.44,  -0.50,
            0.00,   0.00,   0.60,   0.60,   0.84,   0.00,   0.12,  -0.24,
            0.00,   0.80,  -0.56,   0.60,  -0.32,  -0.30,   0.48,  -0.50,
            0.28,  -0.30,  -0.48,  -0.50,   0.12,   0.20,   0.48,  -0.60,
            0.48,   0.60,  -0.12,   0.20,   0.24,   0.00,   0.76,  -0.52,

       /* 3668-3795 */
           -0.60,  -0.52,   0.60,   0.48,  -0.50,  -0.24,  -0.30,   0.12,
           -0.10,   0.48,   0.60,   0.52,  -0.20,   0.36,   0.40,  -0.44,
            0.50,  -0.24,  -0.30,  -0.48,  -0.60,  -0.44,  -0.60,  -0.12,
            0.10,   0.76,   0.76,   0.20,  -0.20,   0.48,   0.50,   0.40,
           -0.50,  -0.24,  -0.30,   0.44,  -0.60,   0.44,  -0.60,   0.36,
            0.00,  -0.64,   0.72,   0.00,  -0.12,   0.00,  -0.10,  -0.40,
           -0.60,  -0.20,  -0.20,  -0.44,   0.50,  -0.44,   0.50,   0.20,
            0.20,  -0.44,  -0.50,   0.20,  -0.20,  -0.20,   0.20,  -0.44,
           -0.50,   0.64,   0.00,   0.32,  -0.36,   0.50,  -0.20,  -0.30,
            0.12,  -0.10,   0.48,   0.50,  -0.12,   0.30,  -0.36,  -0.50,
            0.00,   0.00,   0.48,   0.50,  -0.48,   0.50,   0.68,   0.00,
           -0.12,   0.56,  -0.40,   0.44,  -0.50,  -0.12,  -0.10,   0.24,
            0.30,  -0.40,   0.40,   0.64,   0.00,  -0.24,   0.64,   0.00,
           -0.20,   0.00,   0.00,   0.44,  -0.50,   0.44,   0.50,  -0.12,
            0.20,  -0.36,  -0.50,   0.12,   0.00,   0.64,  -0.40,   0.50,
            0.00,   0.10,   0.00,   0.00,  -0.40,   0.50,   0.00,   0.00,

       /* 3796-3923 */
           -0.40,  -0.50,   0.56,   0.00,   0.28,   0.00,   0.10,   0.36,
            0.50,   0.00,  -0.10,   0.36,  -0.50,   0.36,   0.50,   0.00,
           -0.10,   0.24,  -0.20,  -0.36,  -0.40,   0.16,   0.20,   0.40,
           -0.40,   0.00,   0.00,  -0.36,  -0.50,  -0.36,  -0.50,  -0.32,
           -0.50,  -0.12,   0.10,   0.20,   0.20,  -0.36,   0.40,  -0.60,
            0.60,   0.28,   0.00,   0.52,   0.12,  -0.10,   0.40,   0.40,
            0.00,  -0.50,   0.20,  -0.20,  -0.32,   0.40,   0.16,   0.20,
           -0.16,   0.20,   0.32,   0.40,   0.56,   0.00,  -0.12,   0.32,
           -0.40,  -0.16,  -0.20,   0.00,   0.00,   0.40,   0.40,  -0.40,
           -0.40,  -0.40,   0.40,  -0.36,   0.40,   0.12,   0.10,   0.00,
            0.10,   0.36,   0.40,   0.00,  -0.10,   0.36,   0.40,  -0.36,
            0.40,   0.00,   0.10,   0.32,   0.00,   0.44,   0.12,   0.20,
            0.28,  -0.40,   0.00,   0.00,   0.36,   0.40,   0.32,  -0.40,
           -0.16,   0.12,   0.10,   0.32,  -0.40,   0.20,   0.30,  -0.24,
            0.30,   0.00,   0.10,   0.32,   0.40,   0.00,  -0.10,  -0.32,
           -0.40,  -0.32,   0.40,   0.00,   0.10,  -0.52,  -0.52,   0.52,

       /* 3924-4051 */
            0.32,  -0.40,   0.00,   0.00,   0.32,   0.40,   0.32,  -0.40,
            0.00,   0.00,  -0.32,  -0.40,  -0.32,   0.40,   0.32,   0.40,
            0.00,   0.00,   0.32,   0.40,   0.00,   0.00,  -0.32,  -0.40,
            0.00,   0.00,   0.32,   0.40,   0.16,   0.20,   0.32,  -0.30,
           -0.16,   0.00,  -0.48,  -0.20,   0.20,  -0.28,  -0.30,   0.28,
           -0.40,   0.00,   0.00,   0.28,  -0.40,   0.00,   0.00,   0.28,
           -0.40,   0.00,   0.00,  -0.28,  -0.40,   0.28,   0.40,  -0.28,
           -0.40,  -0.48,  -0.20,   0.20,   0.24,   0.30,   0.44,   0.00,
            0.16,   0.24,   0.30,   0.16,  -0.20,   0.24,   0.30,  -0.12,
            0.20,   0.20,   0.30,  -0.16,   0.20,   0.00,   0.00,   0.44,
           -0.32,   0.30,   0.24,   0.00,  -0.36,   0.36,   0.00,   0.24,
            0.12,  -0.20,   0.20,   0.30,  -0.12,   0.00,  -0.28,   0.30,
           -0.24,   0.30,   0.12,   0.10,  -0.28,  -0.30,  -0.28,   0.30,
            0.00,   0.00,  -0.28,  -0.30,   0.00,   0.00,  -0.28,  -0.30,
            0.00,   0.00,   0.28,   0.30,   0.00,   0.00,  -0.28,  -0.30,
           -0.28,   0.30,   0.00,   0.00,  -0.28,  -0.30,   0.00,   0.00,

       /* 4052-4179 */
            0.28,   0.30,   0.00,   0.00,  -0.28,   0.30,   0.28,  -0.30,
           -0.28,   0.30,   0.40,   0.40,  -0.24,   0.30,   0.00,  -0.10,
            0.16,   0.00,   0.36,  -0.20,   0.30,  -0.12,  -0.10,  -0.24,
           -0.30,   0.00,   0.00,  -0.24,   0.30,  -0.24,   0.30,   0.00,
            0.00,  -0.24,   0.30,  -0.24,   0.30,   0.24,  -0.30,   0.00,
            0.00,   0.24,  -0.30,   0.00,   0.00,   0.24,   0.30,   0.24,
           -0.30,   0.24,   0.30,  -0.24,   0.30,  -0.24,   0.30,  -0.20,
            0.20,  -0.16,  -0.20,   0.00,   0.00,  -0.32,   0.20,   0.00,
            0.10,   0.20,  -0.30,   0.20,  -0.20,   0.12,   0.20,  -0.16,
            0.20,   0.16,   0.20,   0.20,   0.30,   0.20,   0.30,   0.00,
            0.00,  -0.20,   0.30,   0.00,   0.00,   0.20,   0.30,  -0.20,
           -0.30,  -0.20,  -0.30,   0.20,  -0.30,   0.00,   0.00,   0.20,
            0.30,   0.00,   0.00,   0.20,   0.30,   0.00,   0.00,   0.20,
            0.30,   0.00,   0.00,   0.20,   0.30,   0.00,   0.00,   0.20,
           -0.30,   0.00,   0.00,  -0.20,  -0.30,   0.00,   0.00,  -0.20,
            0.30,   0.00,   0.00,  -0.20,   0.30,   0.00,   0.00,   0.36,

       /* 4180-4307 */
            0.00,   0.00,   0.36,   0.12,   0.10,  -0.24,   0.20,   0.12,
           -0.20,  -0.16,  -0.20,  -0.13,   0.10,   0.22,   0.21,   0.20,
            0.00,  -0.28,   0.32,   0.00,  -0.12,  -0.20,  -0.20,   0.12,
           -0.10,   0.12,   0.10,  -0.20,   0.20,   0.00,   0.00,  -0.32,
            0.32,   0.00,   0.00,   0.32,   0.32,   0.00,   0.00,  -0.24,
           -0.20,   0.24,   0.20,   0.20,   0.00,  -0.24,   0.00,   0.00,
           -0.24,  -0.20,   0.00,   0.00,   0.24,   0.20,  -0.24,  -0.20,
            0.00,   0.00,  -0.24,   0.20,   0.16,  -0.20,   0.12,   0.10,
            0.20,   0.20,   0.00,  -0.10,  -0.12,   0.10,  -0.16,  -0.20,
           -0.12,  -0.10,  -0.16,   0.20,   0.20,   0.20,   0.00,   0.00,
           -0.20,   0.20,  -0.20,   0.20,  -0.20,   0.20,  -0.20,   0.20,
            0.20,  -0.20,  -0.20,  -0.20,   0.00,   0.00,  -0.20,   0.20,
            0.20,   0.00,  -0.20,   0.00,   0.00,  -0.20,   0.20,  -0.20,
            0.20,  -0.20,  -0.20,  -0.20,  -0.20,   0.00,   0.00,   0.20,
            0.20,   0.20,   0.20,   0.12,  -0.20,  -0.12,  -0.10,   0.28,
           -0.28,   0.16,  -0.20,   0.00,  -0.10,   0.00,   0.10,  -0.16,

       /* 4308-4435 */
            0.20,   0.00,  -0.10,  -0.16,  -0.20,   0.00,  -0.10,   0.16,
           -0.20,   0.16,  -0.20,   0.00,   0.00,   0.16,   0.20,  -0.16,
            0.20,   0.00,   0.00,   0.16,   0.20,   0.16,  -0.20,   0.16,
           -0.20,  -0.16,   0.20,   0.16,  -0.20,   0.00,   0.00,   0.16,
            0.20,   0.00,   0.00,   0.16,   0.20,   0.00,   0.00,  -0.16,
           -0.20,   0.16,  -0.20,  -0.16,  -0.20,   0.00,   0.00,  -0.16,
           -0.20,   0.00,   0.00,  -0.16,   0.20,   0.00,   0.00,   0.16,
           -0.20,   0.16,   0.20,   0.16,   0.20,   0.00,   0.00,  -0.16,
           -0.20,   0.00,   0.00,  -0.16,  -0.20,   0.00,   0.00,   0.16,
            0.20,   0.16,   0.20,   0.00,   0.00,   0.16,   0.20,   0.16,
           -0.20,   0.16,   0.20,   0.00,   0.00,  -0.16,   0.20,   0.00,
            0.10,   0.12,  -0.20,   0.12,  -0.20,   0.00,  -0.10,   0.00,
           -0.10,   0.12,   0.20,   0.00,  -0.10,  -0.12,   0.20,  -0.15,
            0.20,  -0.24,   0.24,   0.00,   0.00,   0.24,   0.24,   0.12,
           -0.20,  -0.12,  -0.20,   0.00,   0.00,   0.12,   0.20,   0.12,
           -0.20,   0.12,   0.20,   0.12,   0.20,   0.12,   0.20,   0.12,

       /* 4436-4563 */
           -0.20,  -0.12,   0.20,   0.00,   0.00,   0.12,   0.20,   0.12,
            0.00,  -0.20,   0.00,   0.00,  -0.12,  -0.20,   0.12,  -0.20,
            0.00,   0.00,   0.12,   0.20,  -0.12,   0.20,  -0.12,   0.20,
            0.12,  -0.20,   0.00,   0.00,   0.12,   0.20,   0.20,   0.00,
            0.12,   0.00,   0.00,  -0.12,   0.20,   0.00,   0.00,  -0.12,
           -0.20,   0.00,   0.00,  -0.12,  -0.20,  -0.12,  -0.20,   0.00,
            0.00,   0.12,  -0.20,   0.12,  -0.20,   0.12,   0.20,  -0.12,
           -0.20,   0.00,   0.00,   0.12,  -0.20,   0.12,  -0.20,   0.12,
            0.20,   0.12,   0.00,   0.20,  -0.12,  -0.20,   0.00,   0.00,
            0.12,   0.20,  -0.16,   0.00,   0.16,  -0.20,   0.20,   0.00,
            0.00,  -0.20,   0.00,   0.00,  -0.20,   0.20,   0.00,   0.00,
            0.20,   0.20,  -0.20,   0.00,   0.00,  -0.20,   0.12,   0.00,
           -0.16,   0.20,   0.00,   0.00,   0.20,   0.12,  -0.10,   0.00,
            0.10,   0.16,  -0.16,  -0.16,  -0.16,  -0.16,  -0.16,   0.00,
            0.00,  -0.16,   0.00,   0.00,  -0.16,  -0.16,  -0.16,   0.00,
            0.00,  -0.16,   0.00,   0.00,   0.16,   0.00,   0.00,   0.16,

       /* 4564-4691 */
            0.00,   0.00,   0.16,   0.16,   0.00,   0.00,  -0.16,   0.00,
            0.00,  -0.16,  -0.16,   0.00,   0.00,   0.16,   0.00,   0.00,
           -0.16,  -0.16,   0.00,   0.00,  -0.16,  -0.16,   0.12,   0.10,
            0.12,  -0.10,   0.12,   0.10,   0.00,   0.00,   0.12,   0.10,
           -0.12,   0.10,   0.00,   0.00,   0.12,   0.10,   0.12,  -0.10,
            0.00,   0.00,  -0.12,  -0.10,   0.00,   0.00,   0.12,   0.10,
            0.12,   0.00,   0.00,   0.12,   0.00,   0.00,  -0.12,   0.00,
            0.00,   0.12,   0.12,   0.12,   0.12,   0.12,   0.00,   0.00,
            0.12,   0.00,   0.00,   0.12,   0.12,   0.00,   0.00,   0.12,
            0.00,   0.00,   0.12,  -0.12,  -0.12,   0.12,   0.12,  -0.12,
           -0.12,   0.00,   0.00,   0.12,  -0.12,   0.12,   0.12,  -0.12,
           -0.12,   0.00,   0.00,  -0.12,  -0.12,   0.00,   0.00,  -0.12,
            0.12,   0.00,   0.00,   0.12,   0.00,   0.00,   0.12,   0.00,
            0.00,   0.12,  -0.12,   0.00,   0.00,  -0.12,   0.12,  -0.12,
           -0.12,   0.12,   0.00,   0.00,   0.12,   0.12,   0.12,  -0.12,
            0.00,   0.00,  -0.12,  -0.12,  -0.12,   0.00,   0.00,  -0.12,

       /* 4692-NA */
           -0.12,   0.00,   0.00,   0.12,   0.12,   0.00,   0.00,  -0.12,
           -0.12,  -0.12,  -0.12,   0.12,   0.00,   0.00,   0.12,  -0.12,
            0.00,   0.00,  -0.12,  -0.12,   0.00,   0.00,   0.12,  -0.12,
           -0.12,  -0.12,  -0.12,   0.12,   0.12,  -0.12,  -0.12,   0.00,
            0.00,  -0.12,   0.00,   0.00,  -0.12,   0.12,   0.00,   0.00,
            0.12,   0.00,   0.00,  -0.12,  -0.12,   0.00,   0.00,  -0.12,
           -0.12,   0.12,   0.00,   0.00,   0.12,   0.12,   0.00,   0.00,
            0.12,   0.00,   0.00,   0.12,   0.12,   0.08,   0.00,   0.04
       };

        /// <summary>
        /// X,Y coordinates of celestial intermediate pole from series based
        /// on IAU 2006 precession and IAU 2000A nutation.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="x">CIP X coordinates</param>
        /// <param name="y">CIP Y coordinates</param>
        public static void iauXy06(double date1, double date2, out double x, out double y)
        {
            /* Maximum power of T in the polynomials for X and Y */
            //enum { MAXPT = 5 };
            int MAXPT = 5;


            /* Number of frequencies:  luni-solar */
            //static const int NFLS = (int)(sizeof mfals / sizeof(int) / 5);
            int NFLS = Xy06mfals.Length;

            /* Number of frequencies:  planetary */
            //static const int NFPL = (int)(sizeof mfapl / sizeof(int) / 14);
            int NFPL = Xy06mfapl.Length;

            /* Number of amplitude coefficients */
            //static const int NA = (int)(sizeof a / sizeof(double));
            int NA = Xy06a.Length;

            /* Amplitude usage: X or Y, sin or cos, power of T. */
            int[] jaxy = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };
            int[] jasc = { 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0 };
            int[] japt = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };

            /* Miscellaneous */
            double t, w, arg;
            double[] pt = new double[MAXPT + 1];
            double[] fa = new double[14];
            double[] xypr = new double[2];
            double[] xypl = new double[2];
            double[] xyls = new double[2];
            double[] sc = new double[2];
            int jpt, i, j, jxy, ialast, ifreq, m, ia, jsc;

            /* ------------------------------------------------------------------ */

            /* Interval between fundamental date J2000.0 and given date (JC). */
            t = ((date1 - DJ00) + date2) / DJC;

            /* Powers of T. */
            w = 1.0;
            for (jpt = 0; jpt <= MAXPT; jpt++)
            {
                pt[jpt] = w;
                w *= t;
            }

            /* Initialize totals in X and Y:  polynomial, luni-solar, planetary. */
            for (jxy = 0; jxy < 2; jxy++)
            {
                xypr[jxy] = 0.0;
                xyls[jxy] = 0.0;
                xypl[jxy] = 0.0;
            }

            /* --------------------------------- */
            /* Fundamental arguments (IERS 2003) */
            /* --------------------------------- */

            /* Mean anomaly of the Moon. */
            fa[0] = iauFal03(t);

            /* Mean anomaly of the Sun. */
            fa[1] = iauFalp03(t);

            /* Mean argument of the latitude of the Moon. */
            fa[2] = iauFaf03(t);

            /* Mean elongation of the Moon from the Sun. */
            fa[3] = iauFad03(t);

            /* Mean longitude of the ascending node of the Moon. */
            fa[4] = iauFaom03(t);

            /* Planetary longitudes, Mercury through Neptune. */
            fa[5] = iauFame03(t);
            fa[6] = iauFave03(t);
            fa[7] = iauFae03(t);
            fa[8] = iauFama03(t);
            fa[9] = iauFaju03(t);
            fa[10] = iauFasa03(t);
            fa[11] = iauFaur03(t);
            fa[12] = iauFane03(t);

            /* General accumulated precession in longitude. */
            fa[13] = iauFapa03(t);

            /* -------------------------------------- */
            /* Polynomial part of precession-nutation */
            /* -------------------------------------- */

            for (jxy = 0; jxy < 2; jxy++)
            {
                for (j = MAXPT; j >= 0; j--)
                {
                    xypr[jxy] += Xy06xyp[jxy][j] * pt[j];
                }
            }

            /* ---------------------------------- */
            /* Nutation periodic terms, planetary */
            /* ---------------------------------- */

            /* Work backwards through the coefficients per frequency list. */
            ialast = NA;
            for (ifreq = NFPL - 1; ifreq >= 0; ifreq--)
            {

                /* Obtain the argument functions. */
                arg = 0.0;
                for (i = 0; i < 14; i++)
                {
                    m = Xy06mfapl[ifreq][i];
                    if (m != 0) arg += (double)m * fa[i];
                }
                sc[0] = Math.Sin(arg);
                sc[1] = Math.Cos(arg);

                /* Work backwards through the amplitudes at this frequency. */
                ia = Xy06nc[ifreq + NFLS];
                for (i = ialast; i >= ia; i--)
                {

                    /* Coefficient number (0 = 1st). */
                    j = i - ia;

                    /* X or Y. */
                    jxy = jaxy[j];

                    /* Sin or cos. */
                    jsc = jasc[j];

                    /* Power of T. */
                    jpt = japt[j];

                    /* Accumulate the component. */
                    xypl[jxy] += Xy06a[i - 1] * sc[jsc] * pt[jpt];
                }
                ialast = ia - 1;
            }

            /* ----------------------------------- */
            /* Nutation periodic terms, luni-solar */
            /* ----------------------------------- */

            /* Continue working backwards through the number of coefficients list. */
            for (ifreq = NFLS - 1; ifreq >= 0; ifreq--)
            {

                /* Obtain the argument functions. */
                arg = 0.0;
                for (i = 0; i < 5; i++)
                {
                    m = Xy06mfals[ifreq][i];
                    if (m != 0) arg += (double)m * fa[i];
                }
                sc[0] = Math.Sin(arg);
                sc[1] = Math.Cos(arg);

                /* Work backwards through the amplitudes at this frequency. */
                ia = Xy06nc[ifreq];
                for (i = ialast; i >= ia; i--)
                {

                    /* Coefficient number (0 = 1st). */
                    j = i - ia;

                    /* X or Y. */
                    jxy = jaxy[j];

                    /* Sin or cos. */
                    jsc = jasc[j];

                    /* Power of T. */
                    jpt = japt[j];

                    /* Accumulate the component. */
                    xyls[jxy] += Xy06a[i - 1] * sc[jsc] * pt[jpt];
                }
                ialast = ia - 1;
            }

            /* ------------------------------------ */
            /* Results:  CIP unit vector components */
            /* ------------------------------------ */

            x = DAS2R * (xypr[0] + (xyls[0] + xypl[0]) / 1e6);
            y = DAS2R * (xypr[1] + (xyls[1] + xypl[1]) / 1e6);
        }


        /// <summary>
        /// For a given TT date, compute the X,Y coordinates of the Celestial
        /// Intermediate Pole and the CIO locator s, using the IAU 2000A
        /// precession-nutation model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="x">Celestial Intermediate Pole</param>
        /// <param name="y">Celestial Intermediate Pole</param>
        /// <param name="s">the CIO locator s</param>
        public static void iauXys00a(double date1, double date2,
               out double x, out double y, out double s)
        {
            double[][] rbpn;


            /* Form the bias-precession-nutation matrix, IAU 2000A. */
            iauPnm00a(date1, date2, out rbpn);

            /* Extract X,Y. */
            iauBpn2xy(rbpn, out x, out y);

            /* Obtain s. */
            s = iauS00(date1, date2, x, y);
        }


        /// <summary>
        /// For a given TT date, compute the X,Y coordinates of the Celestial
        /// Intermediate Pole and the CIO locator s, using the IAU 2000B
        /// precession-nutation model.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="x">Celestial Intermediate Pole</param>
        /// <param name="y">Celestial Intermediate Pole</param>
        /// <param name="s">the CIO locator s</param>
        public static void iauXys00b(double date1, double date2,
               out double x, out double y, out double s)
        {
            double[][] rbpn;


            /* Form the bias-precession-nutation matrix, IAU 2000A. */
            iauPnm00b(date1, date2, out rbpn);

            /* Extract X,Y. */
            iauBpn2xy(rbpn, out x, out y);

            /* Obtain s. */
            s = iauS00(date1, date2, x, y);
        }


        /// <summary>
        /// For a given TT date, compute the X,Y coordinates of the Celestial
        /// Intermediate Pole and the CIO locator s, using the IAU 2006
        /// precession and IAU 2000A nutation models.
        /// </summary>
        /// <param name="date1">TT as a 2-part Julian Date</param>
        /// <param name="date2">TT as a 2-part Julian Date</param>
        /// <param name="x">Celestial Intermediate Pole</param>
        /// <param name="y">Celestial Intermediate Pole</param>
        /// <param name="s">the CIO locator s</param>
        public static void iauXys06a(double date1, double date2,
               out double x, out double y, out double s)
        {
            double[][] rbpn;


            /* Form the bias-precession-nutation matrix, IAU 2006/2000A. */
            iauPnm06a(date1, date2, out rbpn);

            /* Extract X,Y. */
            iauBpn2xy(rbpn, out x, out y);

            /* Obtain s. */
            s = iauS06(date1, date2, x, y);
        }


    }
}
