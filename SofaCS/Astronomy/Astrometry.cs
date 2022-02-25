using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Apply aberration to transform natural direction into proper direction.
        /// </summary>
        /// <param name="pnat">double[3]   natural direction to the source (unit vector)</param>
        /// <param name="v">double[3]   observer barycentric velocity in units of c</param>
        /// <param name="s">distance between the Sun and the observer (au)</param>
        /// <param name="bm1">sqrt(1-|v|^2): reciprocal of Lorenz factor</param>
        /// <param name="ppr">double[3]   proper direction to source (unit vector)</param>
        public static void iauAb(in double[] pnat, in double[] v, double s, double bm1,
            out double[] ppr)
        {
            int i;
            double pdv, w1, w2, r2, w, r;
            double[] p = new double[3];

            ppr = new double[3];

            pdv = iauPdp(pnat, v);
            w1 = 1.0 + pdv / (1.0 + bm1);
            w2 = SRS / s;
            r2 = 0.0;
            for (i = 0; i < 3; i++)
            {
                w = pnat[i] * bm1 + w1 * v[i] + w2 * (v[i] - pdv * pnat[i]);
                p[i] = w;
                r2 += w * w;
            }
            r = Math.Sqrt(r2);
            for (i = 0; i < 3; i++)
            {
                ppr[i] = p[i] / r;
            }
        }


        /// <summary>
        /// For a geocentric observer, prepare star-independent astrometry
        /// parameters for transformations between ICRS and GCRS coordinates.
        /// The Earth ephemeris is supplied by the caller.
        /// <para></para>
        /// The parameters produced by this function are required in the
        /// parallax, light deflection and aberration parts of the astrometric
        /// transformation chain.
        /// </summary>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="ebpv">double[2][3] Earth barycentric pos/vel (au, au/day)</param>
        /// <param name="ehp">double[3]    Earth heliocentric position (au)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauApcg(double date1, double date2, in double[][] ebpv, in double[] ehp,
            ref iauASTROM astrom)
        {
            /* Geocentric observer */
            double[][] pv = new double[][] {
                new double[] { 0.0, 0.0, 0.0 },
                new double[] { 0.0, 0.0, 0.0 }
            };

            /* Compute the star-independent astrometry parameters. */
            iauApcs(date1, date2, pv, ebpv, ehp, ref astrom);
        }


        /// <summary>
        /// For a geocentric observer, prepare star-independent astrometry
        /// parameters for transformations between ICRS and GCRS coordinates.
        /// The caller supplies the date, and SOFA models are used to predict
        /// the Earth ephemeris.
        /// </summary>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauApcg13(double date1, double date2,
            ref iauASTROM astrom)
        {
            double[][] ehpv;
            double[][] ebpv;


            /* Earth barycentric & heliocentric position/velocity (au, au/d). */
            iauEpv00(date1, date2, out ehpv, out ebpv);

            /* Compute the star-independent astrometry parameters. */
            iauApcg(date1, date2, ebpv, ehpv[0], ref astrom);
        }


        /// <summary>
        /// For a terrestrial observer, prepare star-independent astrometry
        /// parameters for transformations between ICRS and geocentric CIRS
        /// coordinates.  The Earth ephemeris and CIP/CIO are supplied by the caller.
        /// </summary>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="ebpv">double[2][3] Earth barycentric position/velocity (au, au/day)</param>
        /// <param name="ehp">double[3]    Earth heliocentric position (au)</param>
        /// <param name="x">CIP X (components of unit vector)</param>
        /// <param name="y">CIP Y (components of unit vector)</param>
        /// <param name="s">the CIO locator s (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauApci(double date1, double date2,
            in double[][] ebpv, in double[] ehp, double x, double y, double s,
            ref iauASTROM astrom)
        {
            /* Star-independent astrometry parameters for geocenter. */
            iauApcg(date1, date2, ebpv, ehp, ref astrom);

            /* CIO based BPN matrix. */
            iauC2ixys(x, y, s, out astrom.bpn);
        }


        /// <summary>
        /// For a terrestrial observer, prepare star-independent astrometry
        /// parameters for transformations between ICRS and geocentric CIRS
        /// coordinates.  The caller supplies the date, and SOFA models are used
        /// to predict the Earth ephemeris and CIP/CIO.
        /// </summary>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="eo">equation of the origins (ERA-GST)</param>
        public static void iauApci13(double date1, double date2,
            ref iauASTROM astrom, out double eo)
        {
            double[][] ehpv;
            double[][] ebpv;
            double[][] r;
            double x, y, s;


            /* Earth barycentric & heliocentric position/velocity (au, au/d). */
            iauEpv00(date1, date2, out ehpv, out ebpv);

            /* Form the equinox based BPN matrix, IAU 2006/2000A. */
            iauPnm06a(date1, date2, out r);

            /* Extract CIP X,Y. */
            iauBpn2xy(r, out x, out y);

            /* Obtain CIO locator s. */
            s = iauS06(date1, date2, x, y);

            /* Compute the star-independent astrometry parameters. */
            iauApci(date1, date2, ebpv, ehpv[0], x, y, s, ref astrom);

            /* Equation of the origins. */
            eo = iauEors(r, s);
        }


        /// <summary>
        /// For a terrestrial observer, prepare star-independent astrometry
        /// parameters for transformations between ICRS and observed
        /// coordinates.  The caller supplies the Earth ephemeris, the Earth
        /// rotation information and the refraction constants as well as the
        /// site coordinates.
        /// </summary>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="ebpv">double[2][3] Earth barycentric PV (au, au/day)</param>
        /// <param name="ehp">double[3]    Earth heliocentric P (au)</param>
        /// <param name="x">CIP X (components of unit vector)</param>
        /// <param name="y">CIP Y (components of unit vector)</param>
        /// <param name="s">the CIO locator s (radians)</param>
        /// <param name="theta">Earth rotation angle (radians)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="sp">the TIO locator s' (radians)</param>
        /// <param name="refa">refraction constant A (radians)</param>
        /// <param name="refb">refraction constant B (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauApco(double date1, double date2,
             in double[][] ebpv, in double[] ehp,
             double x, double y, double s, double theta,
             double elong, double phi, double hm,
             double xp, double yp, double sp,
             double refa, double refb,
             ref iauASTROM astrom)
        {
            double[][] r;
            double a, b, eral, c;
            double[][] pvc;
            double[][] pv = { new double[3], new double[3] };


            /* Form the rotation matrix, CIRS to apparent [HA,Dec]. */
            iauIr(out r);
            iauRz(theta + sp, r);
            iauRy(-xp, r);
            iauRx(-yp, r);
            iauRz(elong, r);

            /* Solve for local Earth rotation angle. */
            a = r[0][0];
            b = r[0][1];
            eral = (a != 0.0 || b != 0.0) ? Math.Atan2(b, a) : 0.0;
            astrom.eral = eral;

            /* Solve for polar motion [X,Y] with respect to local meridian. */
            a = r[0][0];
            c = r[0][2];
            astrom.xpl = Math.Atan2(c, Math.Sqrt(a * a + b * b));
            a = r[1][2];
            b = r[2][2];
            astrom.ypl = (a != 0.0 || b != 0.0) ? -Math.Atan2(a, b) : 0.0;

            /* Adjusted longitude. */
            astrom.along = iauAnpm(eral - theta);

            /* Functions of latitude. */
            astrom.sphi = Math.Sin(phi);
            astrom.cphi = Math.Cos(phi);

            /* Refraction constants. */
            astrom.refa = refa;
            astrom.refb = refb;

            /* Disable the (redundant) diurnal aberration step. */
            astrom.diurab = 0.0;

            /* CIO based BPN matrix. */
            iauC2ixys(x, y, s, out r);

            /* Observer's geocentric position and velocity (m, m/s, CIRS). */
            iauPvtob(elong, phi, hm, xp, yp, sp, theta, out pvc);

            /* Rotate into GCRS. */
            iauTrxpv(r, pvc, ref pv);

            /* ICRS <-> GCRS parameters. */
            iauApcs(date1, date2, pv, ebpv, ehp, ref astrom);

            /* Store the CIO based BPN matrix. */
            iauCr(r, astrom.bpn);
        }


        /// <summary>
        /// For a terrestrial observer, prepare star-independent astrometry
        /// parameters for transformations between ICRS and observed
        /// coordinates.  The caller supplies UTC, site coordinates, ambient air
        /// conditions and observing wavelength, and SOFA models are used to
        /// obtain the Earth ephemeris, CIP/CIO and refraction constants.
        /// </summary>
        /// <param name="utc1">UTC as a 2-part...</param>
        /// <param name="utc2">...quasi Julian Date</param>
        /// <param name="dut1">UT1-UTC (seconds)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="phpa">pressure at the observer (hPa = mB)</param>
        /// <param name="tc">ambient temperature at the observer (deg C)</param>
        /// <param name="rh">relative humidity at the observer (range 0-1)</param>
        /// <param name="wl">wavelength (micrometers)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="eo">equation of the origins (ERA-GST)</param>
        /// <returns>status: +1 = dubious year
        /// <para>0 = OK</para>
        /// <para>-1 = unacceptable date</para>
        /// </returns>
        public static int iauApco13(double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              out iauASTROM astrom, out double eo)
        {
            int j;
            double tai1, tai2, tt1, tt2, ut11, ut12, x, y, s, theta, sp, refa, refb;
            double[][] ehpv;
            double[][] ebpv;
            double[][] r;

            astrom = new iauASTROM();
            eo = 0.0;

            /* UTC to other time scales. */
            j = iauUtctai(utc1, utc2, out tai1, out tai2);
            if (j < 0) return -1;
            j = iauTaitt(tai1, tai2, out tt1, out tt2);
            j = iauUtcut1(utc1, utc2, dut1, out ut11, out ut12);
            if (j < 0) return -1;

            /* Earth barycentric & heliocentric position/velocity (au, au/d). */
            iauEpv00(tt1, tt2, out ehpv, out ebpv);

            /* Form the equinox based BPN matrix, IAU 2006/2000A. */
            iauPnm06a(tt1, tt2, out r);

            /* Extract CIP X,Y. */
            iauBpn2xy(r, out x, out y);

            /* Obtain CIO locator s. */
            s = iauS06(tt1, tt2, x, y);

            /* Earth rotation angle. */
            theta = iauEra00(ut11, ut12);

            /* TIO locator s'. */
            sp = iauSp00(tt1, tt2);

            /* Refraction constants A and B. */
            iauRefco(phpa, tc, rh, wl, out refa, out refb);

            /* Compute the star-independent astrometry parameters. */
            iauApco(tt1, tt2, ebpv, ehpv[0], x, y, s, theta,
                    elong, phi, hm, xp, yp, sp, refa, refb, ref astrom);

            /* Equation of the origins. */
            eo = iauEors(r, s);

            /* Return any warning status. */
            return j;
        }


        /// <summary>
        /// For an observer whose geocentric position and velocity are known,
        /// prepare star-independent astrometry parameters for transformations
        /// between ICRS and GCRS.  The Earth ephemeris is supplied by the caller.
        /// </summary>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="pv">double[2][3] observer's geocentric pos/vel (m, m/s)</param>
        /// <param name="ebpv">double[2][3] Earth barycentric PV (au, au/day)</param>
        /// <param name="ehp">double[3]    Earth heliocentric P (au)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauApcs(double date1, double date2,
            in double[][] pv, in double[][] ebpv, in double[] ehp,
            ref iauASTROM astrom)
        {
            /* au/d to m/s */
            double AUDMS = DAU / DAYSEC;

            /* Light time for 1 au (day) */
            double CR = AULT / DAYSEC;

            int i;
            double dp, dv, v2, w;
            double[] pb = new double[3];
            double[] vb = new double[3];
            double[] ph = new double[3];


            /* Time since reference epoch, years (for proper motion calculation). */
            astrom.pmt = ((date1 - DJ00) + date2) / DJY;

            /* Adjust Earth ephemeris to observer. */
            for (i = 0; i < 3; i++)
            {
                dp = pv[0][i] / DAU;
                dv = pv[1][i] / AUDMS;
                pb[i] = ebpv[0][i] + dp;
                vb[i] = ebpv[1][i] + dv;
                ph[i] = ehp[i] + dp;
            }

            /* Barycentric position of observer (au). */
            iauCp(pb, astrom.eb);

            /* Heliocentric direction and distance (unit vector and au). */
            iauPn(ph, out astrom.em, ref astrom.eh);

            /* Barycentric vel. in units of c, and reciprocal of Lorenz factor. */
            v2 = 0.0;
            for (i = 0; i < 3; i++)
            {
                w = vb[i] * CR;
                astrom.v[i] = w;
                v2 += w * w;
            }
            astrom.bm1 = Math.Sqrt(1.0 - v2);

            /* Reset the NPB matrix. */
            iauIr(out astrom.bpn);
        }


        /// <summary>
        /// For an observer whose geocentric position and velocity are known,
        /// prepare star-independent astrometry parameters for transformations
        /// between ICRS and GCRS.  The Earth ephemeris is from SOFA models.
        /// </summary>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="pv">double[2][3] observer's geocentric pos/vel)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauApcs13(double date1, double date2, in double[][] pv,
            ref iauASTROM astrom)
        {
            double[][] ehpv;
            double[][] ebpv;


            /* Earth barycentric & heliocentric position/velocity (au, au/d). */
            iauEpv00(date1, date2, out ehpv, out ebpv);

            /* Compute the star-independent astrometry parameters. */
            iauApcs(date1, date2, pv, ebpv, ehpv[0], ref astrom);
        }


        /// <summary>
        /// In the star-independent astrometry parameters, update only the
        /// Earth rotation angle, supplied by the caller explicitly.
        /// </summary>
        /// <param name="theta">Earth rotation angle (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauAper(double theta, iauASTROM astrom)
        {
            astrom.eral = theta + astrom.along;
        }


        /// <summary>
        /// In the star-independent astrometry parameters, update only the
        /// Earth rotation angle.  The caller provides UT1, (n.b. not UTC).
        /// </summary>
        /// <param name="ut11">UT1 as a 2-part...</param>
        /// <param name="ut12">...Julian Date</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauAper13(double ut11, double ut12, iauASTROM astrom)
        {
            iauAper(iauEra00(ut11, ut12), astrom);
        }


        /// <summary>
        /// For a terrestrial observer, prepare star-independent astrometry
        /// parameters for transformations between CIRS and observed
        /// coordinates.  The caller supplies the Earth orientation information
        /// and the refraction constants as well as the site coordinates.
        /// </summary>
        /// <param name="sp">the TIO locator s' (radians)</param>
        /// <param name="theta">Earth rotation angle (radians)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">geodetic latitude (radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic2)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="refa">refraction constant A (radians)</param>
        /// <param name="refb">refraction constant B (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        public static void iauApio(double sp, double theta,
             double elong, double phi, double hm, double xp, double yp,
             double refa, double refb,
             out iauASTROM astrom)
        {
            double[][] r;
            double a, b, eral, c;
            double[][] pv;

            astrom = new iauASTROM();

            /* Form the rotation matrix, CIRS to apparent [HA,Dec]. */
            iauIr(out r);
            iauRz(theta + sp, r);
            iauRy(-xp, r);
            iauRx(-yp, r);
            iauRz(elong, r);

            /* Solve for local Earth rotation angle. */
            a = r[0][0];
            b = r[0][1];
            eral = (a != 0.0 || b != 0.0) ? Math.Atan2(b, a) : 0.0;
            astrom.eral = eral;

            /* Solve for polar motion [X,Y] with respect to local meridian. */
            a = r[0][0];
            c = r[0][2];
            astrom.xpl = Math.Atan2(c, Math.Sqrt(a * a + b * b));
            a = r[1][2];
            b = r[2][2];
            astrom.ypl = (a != 0.0 || b != 0.0) ? -Math.Atan2(a, b) : 0.0;

            /* Adjusted longitude. */
            astrom.along = iauAnpm(eral - theta);

            /* Functions of latitude. */
            astrom.sphi = Math.Sin(phi);
            astrom.cphi = Math.Cos(phi);

            /* Observer's geocentric position and velocity (m, m/s, CIRS). */
            iauPvtob(elong, phi, hm, xp, yp, sp, theta, out pv);

            /* Magnitude of diurnal aberration vector. */
            astrom.diurab = Math.Sqrt(pv[1][0] * pv[1][0] + pv[1][1] * pv[1][1]) / CMPS;

            /* Refraction constants. */
            astrom.refa = refa;
            astrom.refb = refb;
        }


        /// <summary>
        /// For a terrestrial observer, prepare star-independent astrometry
        /// parameters for transformations between CIRS and observed
        /// coordinates.  The caller supplies UTC, site coordinates, ambient air
        /// conditions and observing wavelength.
        /// </summary>
        /// <param name="utc1">UTC as a 2-part...</param>
        /// <param name="utc2">...quasi Julian Date</param>
        /// <param name="dut1">UT1-UTC (seconds)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">geodetic latitude (radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="phpa">pressure at the observer (hPa = mB)</param>
        /// <param name="tc">ambient temperature at the observer (deg C)</param>
        /// <param name="rh">relative humidity at the observer (range 0-1)</param>
        /// <param name="wl">wavelength (micrometers)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <returns>status: +1 = dubious year
        /// <para>0 = OK</para>
        /// <para>-1 = unacceptable date</para>
        /// </returns>
        public static int iauApio13(double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              out iauASTROM astrom)
        {
            int j;
            double tai1, tai2, tt1, tt2, ut11, ut12, sp, theta, refa, refb;

            astrom = new iauASTROM();

            /* UTC to other time scales. */
            j = iauUtctai(utc1, utc2, out tai1, out tai2);
            if (j < 0) return -1;
            j = iauTaitt(tai1, tai2, out tt1, out tt2);
            j = iauUtcut1(utc1, utc2, dut1, out ut11, out ut12);
            if (j < 0) return -1;

            /* TIO locator s'. */
            sp = iauSp00(tt1, tt2);

            /* Earth rotation angle. */
            theta = iauEra00(ut11, ut12);

            /* Refraction constants A and B. */
            iauRefco(phpa, tc, rh, wl, out refa, out refb);

            /* CIRS <-> observed astrometry parameters. */
            iauApio(sp, theta, elong, phi, hm, xp, yp, refa, refb, out astrom);

            /* Return any warning status. */
            return j;
        }


        /// <summary>
        /// Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS
        /// astrometric place.
        /// </summary>
        /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
        /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
        /// <param name="pr">RA proper motion (radians/year)</param>
        /// <param name="pd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="ra">ICRS astrometric RA (radians)</param>
        /// <param name="da">ICRS astrometric Dec (radians)</param>
        public static void iauAtcc13(double rc, double dc,
               double pr, double pd, double px, double rv,
               double date1, double date2,
               out double ra, out double da)
        {
            /* Star-independent astrometry parameters */
            iauASTROM astrom = new iauASTROM();

            double w;


            /* The transformation parameters. */
            iauApci13(date1, date2, ref astrom, out w);

            /* Catalog ICRS (epoch J2000.0) to astrometric. */
            iauAtccq(rc, dc, pr, pd, px, rv, astrom, out ra, out da);
        }


        /// <summary>
        /// Quick transformation of a star's ICRS catalog entry (epoch J2000.0)
        /// into ICRS astrometric place, given precomputed star-independent
        /// astrometry parameters.
        /// </summary>
        /// <param name="rc">ICRS RA at J2000.0 (radians)</param>
        /// <param name="dc">ICRS Dec at J2000.0 (radians)</param>
        /// <param name="pr">RA proper motion (radians/year)</param>
        /// <param name="pd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
        /// <param name="astrom">star-independent astrometry parameters:</param>
        /// <param name="ra">ICRS astrometric RA (radians)</param>
        /// <param name="da">ICRS astrometric Dec (radians)</param>
        public static void iauAtccq(double rc, double dc,
              double pr, double pd, double px, double rv,
              in iauASTROM astrom, out double ra, out double da)
        {
            double[] p;
            double w;


            /* Proper motion and parallax, giving BCRS coordinate direction. */
            iauPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb, out p);

            /* ICRS astrometric RA,Dec. */
            iauC2s(p, out w, out da);
            ra = iauAnp(w);
        }


        /// <summary>
        /// Transform ICRS star data, epoch J2000.0, to CIRS.
        /// </summary>
        /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
        /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
        /// <param name="pr">RA proper motion (radians/year)</param>
        /// <param name="pd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="ri">CIRS geocentric RA (radians)</param>
        /// <param name="di">CIRS geocentric Dec (radians)</param>
        /// <param name="eo">equation of the origins (ERA-GST)</param>
        public static void iauAtci13(double rc, double dc,
               double pr, double pd, double px, double rv,
               double date1, double date2,
               out double ri, out double di, out double eo)
        {
            /* Star-independent astrometry parameters */
            iauASTROM astrom = new iauASTROM();


            /* The transformation parameters. */
            iauApci13(date1, date2, ref astrom, out eo);

            /* ICRS (epoch J2000.0) to CIRS. */
            iauAtciq(rc, dc, pr, pd, px, rv, astrom, out ri, out di);
        }


        /// <summary>
        /// Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
        /// star-independent astrometry parameters.
        /// </summary>
        /// <param name="rc">ICRS RA at J2000.0 (radians)</param>
        /// <param name="dc">ICRS Dec at J2000.0 (radians)</param>
        /// <param name="pr">RA proper motion (radians/year)</param>
        /// <param name="pd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
        /// <param name="astrom">star-independent astrometry parameters:</param>
        /// <param name="ri">CIRS RA (radians)</param>
        /// <param name="di">CIRS Dec (radians)</param>
        public static void iauAtciq(double rc, double dc,
              double pr, double pd, double px, double rv,
              in iauASTROM astrom, out double ri, out double di)
        {
            double w;
            double[] pco;
            double[] pnat = new double[3];
            double[] ppr;
            double[] pi = new double[3];

            /* Proper motion and parallax, giving BCRS coordinate direction. */
            iauPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb, out pco);

            /* Light deflection by the Sun, giving BCRS natural direction. */
            iauLdsun(pco, astrom.eh, astrom.em, ref pnat);

            /* Aberration, giving GCRS proper direction. */
            iauAb(pnat, astrom.v, astrom.em, astrom.bm1, out ppr);

            /* Bias-precession-nutation, giving CIRS proper direction. */
            iauRxp(astrom.bpn, ppr, ref pi);

            /* CIRS RA,Dec. */
            iauC2s(pi, out w, out di);
            ri = iauAnp(w);
        }


        /// <summary>
        /// Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
        /// star-independent astrometry parameters plus a list of light-deflecting bodies.
        /// </summary>
        /// <param name="rc">ICRS RA at J2000.0 (radians)</param>
        /// <param name="dc">ICRS Dec at J2000.0 (radians)</param>
        /// <param name="pr">RA proper motion (radians/year)</param>
        /// <param name="pd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="n">number of bodies</param>
        /// <param name="b">data for each of the n bodies</param>
        /// <param name="ri">CIRS RA (radians)</param>
        /// <param name="di">CIRS Dec (radians)</param>
        public static void iauAtciqn(double rc, double dc, double pr, double pd,
               double px, double rv, in iauASTROM astrom,
               int n, in iauLDBODY[] b, out double ri, out double di)
        {
            double w;
            double[] pco;
            double[] pnat = new double[3];
            double[] ppr;
            double[] pi = new double[3];

            /* Proper motion and parallax, giving BCRS coordinate direction. */
            iauPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb, out pco);

            /* Light deflection, giving BCRS natural direction. */
            iauLdn(n, b, astrom.eb, pco, ref pnat);

            /* Aberration, giving GCRS proper direction. */
            iauAb(pnat, astrom.v, astrom.em, astrom.bm1, out ppr);

            /* Bias-precession-nutation, giving CIRS proper direction. */
            iauRxp(astrom.bpn, ppr, ref pi);

            /* CIRS RA,Dec. */
            iauC2s(pi, out w, out di);
            ri = iauAnp(w);
        }


        /// <summary>
        /// Quick ICRS to CIRS transformation, given precomputed 
        /// star-independent astrometry parameters, and assuming zero parallax and
        /// proper motion.
        /// </summary>
        /// <param name="rc">ICRS astrometric RA (radians)</param>
        /// <param name="dc">ICRS astrometric Dec (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="ri">CIRS RA (radians)</param>
        /// <param name="di">CIRS Dec (radians)</param>
        public static void iauAtciqz(double rc, double dc, in iauASTROM astrom,
               out double ri, out double di)
        {
            double w;
            double[] pco;
            double[] pnat = new double[3];
            double[] ppr;
            double[] pi = new double[3];


            /* BCRS coordinate direction (unit vector). */
            iauS2c(rc, dc, out pco);

            /* Light deflection by the Sun, giving BCRS natural direction. */
            iauLdsun(pco, astrom.eh, astrom.em, ref pnat);

            /* Aberration, giving GCRS proper direction. */
            iauAb(pnat, astrom.v, astrom.em, astrom.bm1, out ppr);

            /* Bias-precession-nutation, giving CIRS proper direction. */
            iauRxp(astrom.bpn, ppr, ref pi);

            /* CIRS RA,Dec. */
            iauC2s(pi, out w, out di);
            ri = iauAnp(w);
        }


        /// <summary>
        /// ICRS RA,Dec to observed place.  The caller supplies UTC, site
        /// coordinates, ambient air conditions and observing wavelength.
        /// </summary>
        /// <param name="rc">ICRS RA at J2000.0 (radians)</param>
        /// <param name="dc">ICRS Dec at J2000.0 (radians)</param>
        /// <param name="pr">RA proper motion (radians/year)</param>
        /// <param name="pd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
        /// <param name="utc1">UTC as a 2-part...</param>
        /// <param name="utc2">...quasi Julian Date</param>
        /// <param name="dut1">UT1-UTC (seconds)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="phpa">pressure at the observer (hPa = mB)</param>
        /// <param name="tc">ambient temperature at the observer (deg C)</param>
        /// <param name="rh">relative humidity at the observer (range 0-1)</param>
        /// <param name="wl">wavelength (micrometers)</param>
        /// <param name="aob">observed azimuth (radians: N=0,E=90)</param>
        /// <param name="zob">observed zenith distance (radians)</param>
        /// <param name="hob">observed hour angle (radians)</param>
        /// <param name="dob">observed declination (radians)</param>
        /// <param name="rob">observed right ascension (CIO-based, radians)</param>
        /// <param name="eo">equation of the origins (ERA-GST)</param>
        /// <returns>status: +1 = dubious year
        /// <para>0 = OK</para>
        /// <para>-1 = unacceptable date</para>
        /// </returns>
        public static int iauAtco13(double rc, double dc,
              double pr, double pd, double px, double rv,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              out double aob, out double zob, out double hob,
              out double dob, out double rob, out double eo)
        {
            int j;
            iauASTROM astrom;
            double ri, di;

            aob = 0;
            zob = 0;
            hob = 0;
            dob = 0;
            rob = 0;

            /* Star-independent astrometry parameters. */
            j = iauApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                          phpa, tc, rh, wl, out astrom, out eo);

            /* Abort if bad UTC. */
            if (j < 0) return j;

            /* Transform ICRS to CIRS. */
            iauAtciq(rc, dc, pr, pd, px, rv, astrom, out ri, out di);

            /* Transform CIRS to observed. */
            iauAtioq(ri, di, astrom, out aob, out zob, out hob, out dob, out rob);

            /* Return OK/warning status. */
            return j;
        }


        /// <summary>
        /// Transform star RA,Dec from geocentric CIRS to ICRS astrometric.
        /// </summary>
        /// <param name="ri">CIRS geocentric RA (radians)</param>
        /// <param name="di">CIRS geocentric Dec (radians)</param>
        /// <param name="date1">TDB as a 2-part...</param>
        /// <param name="date2">...Julian Date</param>
        /// <param name="rc">ICRS astrometric RA (radians)</param>
        /// <param name="dc">ICRS astrometric Dec (radians)</param>
        /// <param name="eo">equation of the origins (ERA-GST)</param>
        public static void iauAtic13(double ri, double di, double date1, double date2,
               out double rc, out double dc, out double eo)
        {
            /* Star-independent astrometry parameters */
            iauASTROM astrom = new iauASTROM();


            /* Star-independent astrometry parameters. */
            iauApci13(date1, date2, ref astrom, out eo);

            /* CIRS to ICRS astrometric. */
            iauAticq(ri, di, astrom, out rc, out dc);
        }


        /// <summary>
        /// Quick CIRS RA,Dec to ICRS astrometric place, given the star-independent
        /// astrometry parameters.
        /// </summary>
        /// <param name="ri">CIRS RA (radians)</param>
        /// <param name="di">CIRS Dec (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="rc">ICRS astrometric RA (radians)</param>
        /// <param name="dc">ICRS astrometric Dec (radians)</param>
        public static void iauAticq(double ri, double di, in iauASTROM astrom,
              out double rc, out double dc)
        {
            int j, i;
            double w, r2, r;
            double[] pi;
            double[] ppr = new double[3];
            double[] pnat = new double[3];
            double[] pco = new double[3];
            double[] d;
            double[] before = new double[3];
            double[] after = new double[3];

            /* CIRS RA,Dec to Cartesian. */
            iauS2c(ri, di, out pi);

            /* Bias-precession-nutation, giving GCRS proper direction. */
            iauTrxp(astrom.bpn, pi, ref ppr);

            /* Aberration, giving GCRS natural direction. */
            iauZp(out d);
            for (j = 0; j < 2; j++)
            {
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    w = ppr[i] - d[i];
                    before[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    before[i] /= r;
                }
                iauAb(before, astrom.v, astrom.em, astrom.bm1, out after);
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    d[i] = after[i] - before[i];
                    w = ppr[i] - d[i];
                    pnat[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    pnat[i] /= r;
                }
            }

            /* Light deflection by the Sun, giving BCRS coordinate direction. */
            iauZp(out d);
            for (j = 0; j < 5; j++)
            {
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    w = pnat[i] - d[i];
                    before[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    before[i] /= r;
                }
                iauLdsun(before, astrom.eh, astrom.em, ref after);
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    d[i] = after[i] - before[i];
                    w = pnat[i] - d[i];
                    pco[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    pco[i] /= r;
                }
            }

            /* ICRS astrometric RA,Dec. */
            iauC2s(pco, out w, out dc);
            rc = iauAnp(w);
        }


        /// <summary>
        /// Quick CIRS to ICRS astrometric place transformation, given the star-independent
        /// astrometry parameters plus a list of light-deflecting bodies.
        /// </summary>
        /// <param name="ri">CIRS RA (radians)</param>
        /// <param name="di">CIRS Dec (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="n">number of bodies</param>
        /// <param name="b">data for each of the n bodies</param>
        /// <param name="rc">ICRS astrometric RA (radians)</param>
        /// <param name="dc">ICRS astrometric Dec (radians)</param>
        public static void iauAticqn(double ri, double di, in iauASTROM astrom,
               int n, in iauLDBODY[] b, out double rc, out double dc)
        {
            int j, i;
            double w, r2, r;
            double[] pi;
            double[] ppr = new double[3];
            double[] pnat = new double[3];
            double[] pco = new double[3];
            double[] d;
            double[] before = new double[3];
            double[] after = new double[3];


            /* CIRS RA,Dec to Cartesian. */
            iauS2c(ri, di, out pi);

            /* Bias-precession-nutation, giving GCRS proper direction. */
            iauTrxp(astrom.bpn, pi, ref ppr);

            /* Aberration, giving GCRS natural direction. */
            iauZp(out d);
            for (j = 0; j < 2; j++)
            {
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    w = ppr[i] - d[i];
                    before[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    before[i] /= r;
                }
                iauAb(before, astrom.v, astrom.em, astrom.bm1, out after);
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    d[i] = after[i] - before[i];
                    w = ppr[i] - d[i];
                    pnat[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    pnat[i] /= r;
                }
            }

            /* Light deflection, giving BCRS coordinate direction. */
            iauZp(out d);
            for (j = 0; j < 5; j++)
            {
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    w = pnat[i] - d[i];
                    before[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    before[i] /= r;
                }
                iauLdn(n, b, astrom.eb, before, ref after);
                r2 = 0.0;
                for (i = 0; i < 3; i++)
                {
                    d[i] = after[i] - before[i];
                    w = pnat[i] - d[i];
                    pco[i] = w;
                    r2 += w * w;
                }
                r = Math.Sqrt(r2);
                for (i = 0; i < 3; i++)
                {
                    pco[i] /= r;
                }
            }

            /* ICRS astrometric RA,Dec. */
            iauC2s(pco, out w, out dc);
            rc = iauAnp(w);
        }


        /// <summary>
        /// CIRS RA,Dec to observed place.  The caller supplies UTC, site
        /// coordinates, ambient air conditions and observing wavelength.
        /// </summary>
        /// <param name="ri">CIRS right ascension (CIO-based, radians)</param>
        /// <param name="di">CIRS declination (radians)</param>
        /// <param name="utc1">UTC as a 2-part...</param>
        /// <param name="utc2">...quasi Julian Date</param>
        /// <param name="dut1">UT1-UTC (seconds)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">geodetic latitude (radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="phpa">pressure at the observer (hPa = mB)</param>
        /// <param name="tc">ambient temperature at the observer (deg C)</param>
        /// <param name="rh">relative humidity at the observer (range 0-1)</param>
        /// <param name="wl">wavelength (micrometers)</param>
        /// <param name="aob">observed azimuth (radians: N=0,E=90)</param>
        /// <param name="zob">observed zenith distance (radians)</param>
        /// <param name="hob">observed hour angle (radians)</param>
        /// <param name="dob">observed declination (radians)</param>
        /// <param name="rob">observed right ascension (CIO-based, radians)</param>
        /// <returns>status: +1 = dubious year
        /// <para>0 = OK</para>
        /// <para>-1 = unacceptable date</para>
        /// </returns>
        public static int iauAtio13(double ri, double di,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              out double aob, out double zob, out double hob,
              out double dob, out double rob)
        {
            int j;
            iauASTROM astrom;

            aob = 0; zob = 0; hob = 0; dob = 0; rob = 0;

            /* Star-independent astrometry parameters for CIRS->observed. */
            j = iauApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                          phpa, tc, rh, wl, out astrom);

            /* Abort if bad UTC. */
            if (j < 0) return j;

            /* Transform CIRS to observed. */
            iauAtioq(ri, di, astrom, out aob, out zob, out hob, out dob, out rob);

            /* Return OK/warning status. */
            return j;
        }


        /// <summary>
        /// Quick CIRS to observed place transformation.
        /// </summary>
        /// <param name="ri">CIRS right ascension</param>
        /// <param name="di">CIRS declination</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="aob">observed azimuth (radians: N=0,E=90)</param>
        /// <param name="zob">observed zenith distance (radians)</param>
        /// <param name="hob">observed hour angle (radians)</param>
        /// <param name="dob">observed declination (radians)</param>
        /// <param name="rob">observed right ascension (CIO-based, radians)</param>
        public static void iauAtioq(double ri, double di, in iauASTROM astrom,
              out double aob, out double zob,
              out double hob, out double dob, out double rob)
        {
            /* Minimum cos(alt) and sin(alt) for refraction purposes */
            double CELMIN = 1e-6;
            double SELMIN = 0.05;

            double x, y, z, sx, cx, sy, cy, xhd, yhd, zhd, f,
                   xhdt, yhdt, zhdt, xaet, yaet, zaet, azobs, r, tz, w, del,
                   cosdel, xaeo, yaeo, zaeo, zdobs, hmobs, dcobs, raobs;
            double[] v;


            /* CIRS RA,Dec to Cartesian -HA,Dec. */
            iauS2c(ri - astrom.eral, di, out v);
            x = v[0];
            y = v[1];
            z = v[2];

            /* Polar motion. */
            sx = Math.Sin(astrom.xpl);
            cx = Math.Cos(astrom.xpl);
            sy = Math.Sin(astrom.ypl);
            cy = Math.Cos(astrom.ypl);
            xhd = cx * x + sx * z;
            yhd = sx * sy * x + cy * y - cx * sy * z;
            zhd = -sx * cy * x + sy * y + cx * cy * z;

            /* Diurnal aberration. */
            f = (1.0 - astrom.diurab * yhd);
            xhdt = f * xhd;
            yhdt = f * (yhd + astrom.diurab);
            zhdt = f * zhd;

            /* Cartesian -HA,Dec to Cartesian Az,El (S=0,E=90). */
            xaet = astrom.sphi * xhdt - astrom.cphi * zhdt;
            yaet = yhdt;
            zaet = astrom.cphi * xhdt + astrom.sphi * zhdt;

            /* Azimuth (N=0,E=90). */
            azobs = (xaet != 0.0 || yaet != 0.0) ? Math.Atan2(yaet, -xaet) : 0.0;

            /* ---------- */
            /* Refraction */
            /* ---------- */

            /* Cosine and sine of altitude, with precautions. */
            r = Math.Sqrt(xaet * xaet + yaet * yaet);
            r = r > CELMIN ? r : CELMIN;
            z = zaet > SELMIN ? zaet : SELMIN;

            /* A*tan(z)+B*tan^3(z) model, with Newton-Raphson correction. */
            tz = r / z;
            w = astrom.refb * tz * tz;
            del = (astrom.refa + w) * tz /
                  (1.0 + (astrom.refa + 3.0 * w) / (z * z));

            /* Apply the change, giving observed vector. */
            cosdel = 1.0 - del * del / 2.0;
            f = cosdel - del * z / r;
            xaeo = xaet * f;
            yaeo = yaet * f;
            zaeo = cosdel * zaet + del * r;

            /* Observed ZD. */
            zdobs = Math.Atan2(Math.Sqrt(xaeo * xaeo + yaeo * yaeo), zaeo);

            /* Az/El vector to HA,Dec vector (both right-handed). */
            v[0] = astrom.sphi * xaeo + astrom.cphi * zaeo;
            v[1] = yaeo;
            v[2] = -astrom.cphi * xaeo + astrom.sphi * zaeo;

            /* To spherical -HA,Dec. */
            iauC2s(v, out hmobs, out dcobs);

            /* Right ascension (with respect to CIO). */
            raobs = astrom.eral + hmobs;

            /* Return the results. */
            aob = iauAnp(azobs);
            zob = zdobs;
            hob = -hmobs;
            dob = dcobs;
            rob = iauAnp(raobs);
        }


        /// <summary>
        /// Observed place at a groundbased site to to ICRS astrometric RA,Dec.
        /// The caller supplies UTC, site coordinates, ambient air conditions
        /// and observing wavelength.
        /// </summary>
        /// <param name="type">type of coordinates - "R", "H" or "A"</param>
        /// <param name="ob1">observed Az, HA or RA (radians; Az is N=0,E=90)</param>
        /// <param name="ob2">observed ZD or Dec (radians)</param>
        /// <param name="utc1">UTC as a 2-part...</param>
        /// <param name="utc2">...quasi Julian Date</param>
        /// <param name="dut1">UT1-UTC (seconds)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">geodetic latitude (radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="phpa">pressure at the observer (hPa = mB)</param>
        /// <param name="tc">ambient temperature at the observer (deg C)</param>
        /// <param name="rh">relative humidity at the observer (range 0-1)</param>
        /// <param name="wl">wavelength (micrometers)</param>
        /// <param name="rc">ICRS astrometric RA (radians)</param>
        /// <param name="dc">ICRS astrometric Dec (radians)</param>
        /// <returns>status: +1 = dubious year
        /// <para>0 = OK</para>
        /// <para>-1 = unacceptable date</para>
        /// </returns>
        public static int iauAtoc13(string type, double ob1, double ob2,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              out double rc, out double dc)
        {
            int j;
            iauASTROM astrom;
            double eo, ri, di;

            rc = 0;
            dc = 0;

            /* Star-independent astrometry parameters. */
            j = iauApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                          phpa, tc, rh, wl, out astrom, out eo);

            /* Abort if bad UTC. */
            if (j < 0) return j;

            /* Transform observed to CIRS. */
            iauAtoiq(type, ob1, ob2, astrom, out ri, out di);

            /* Transform CIRS to ICRS. */
            iauAticq(ri, di, astrom, out rc, out dc);

            /* Return OK/warning status. */
            return j;
        }


        /// <summary>
        /// Observed place to CIRS.  The caller supplies UTC, site coordinates,
        /// ambient air conditions and observing wavelength.
        /// </summary>
        /// <param name="type">type of coordinates - "R", "H" or "A"</param>
        /// <param name="ob1">observed Az, HA or RA (radians; Az is N=0,E=90)</param>
        /// <param name="ob2">observed ZD or Dec (radians)</param>
        /// <param name="utc1">UTC</param>
        /// <param name="utc2">UTC</param>
        /// <param name="dut1">UT1-UTC (seconds)</param>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">geodetic latitude (radians)</param>
        /// <param name="hm">height above ellipsoid (m, geodetic)</param>
        /// <param name="xp">polar motion coordinates (radians)</param>
        /// <param name="yp">polar motion coordinates (radians)</param>
        /// <param name="phpa">pressure at the observer (hPa = mB)</param>
        /// <param name="tc">ambient temperature at the observer (deg C)</param>
        /// <param name="rh">relative humidity at the observer (range 0-1)</param>
        /// <param name="wl">wavelength (micrometers)</param>
        /// <param name="ri">CIRS right ascension (CIO-based, radians)</param>
        /// <param name="di">CIRS declination (radians)</param>
        /// <returns>status: +1 = dubious year
        /// <para>0 = OK</para>
        /// <para>-1 = unacceptable date</para>
        /// </returns>
        public static int iauAtoi13(string type, double ob1, double ob2,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              out double ri, out double di)
        {
            int j;
            iauASTROM astrom;

            ri = 0;
            di = 0;

            /* Star-independent astrometry parameters for CIRS->observed. */
            j = iauApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                          phpa, tc, rh, wl, out astrom);

            /* Abort if bad UTC. */
            if (j < 0) return j;

            /* Transform observed to CIRS. */
            iauAtoiq(type, ob1, ob2, astrom, out ri, out di);

            /* Return OK/warning status. */
            return j;
        }


        /// <summary>
        /// Quick observed place to CIRS, given the star-independent astrometry parameters.
        /// </summary>
        /// <param name="type">type of coordinates - "R", "H" or "A"</param>
        /// <param name="ob1">observed Az, HA or RA (radians; Az is N=0,E=90)</param>
        /// <param name="ob2">observed ZD or Dec (radians)</param>
        /// <param name="astrom">star-independent astrometry parameters</param>
        /// <param name="ri">CIRS right ascension (CIO-based, radians)</param>
        /// <param name="di">CIRS declination (radians)</param>
        public static void iauAtoiq(string type,
              double ob1, double ob2, in iauASTROM astrom,
              out double ri, out double di)
        {
            /* Minimum sin(alt) for refraction purposes */
            double SELMIN = 0.05;

            char c;
            double c1, c2, sphi, cphi, ce, xaeo, yaeo, zaeo,
                   xmhdo, ymhdo, zmhdo, az, sz, zdo, refa, refb, tz, dref,
                   zdt, xaet, yaet, zaet, xmhda, ymhda, zmhda,
                   f, xhd, yhd, zhd, sx, cx, sy, cy, hma;
            double[] v = new double[3];


            /* Coordinate type. */
            c = type[0];

            /* Coordinates. */
            c1 = ob1;
            c2 = ob2;

            /* Sin, cos of latitude. */
            sphi = astrom.sphi;
            cphi = astrom.cphi;

            /* Standardize coordinate type. */
            if (c == 'r' || c == 'R')
            {
                c = 'R';
            }
            else if (c == 'h' || c == 'H')
            {
                c = 'H';
            }
            else
            {
                c = 'A';
            }

            /* If Az,ZD, convert to Cartesian (S=0,E=90). */
            if (c == 'A')
            {
                ce = Math.Sin(c2);
                xaeo = -Math.Cos(c1) * ce;
                yaeo = Math.Sin(c1) * ce;
                zaeo = Math.Cos(c2);

            }
            else
            {

                /* If RA,Dec, convert to HA,Dec. */
                if (c == 'R') c1 = astrom.eral - c1;

                /* To Cartesian -HA,Dec. */
                iauS2c(-c1, c2, out v);
                xmhdo = v[0];
                ymhdo = v[1];
                zmhdo = v[2];

                /* To Cartesian Az,El (S=0,E=90). */
                xaeo = sphi * xmhdo - cphi * zmhdo;
                yaeo = ymhdo;
                zaeo = cphi * xmhdo + sphi * zmhdo;
            }

            /* Azimuth (S=0,E=90). */
            az = (xaeo != 0.0 || yaeo != 0.0) ? Math.Atan2(yaeo, xaeo) : 0.0;

            /* Sine of observed ZD, and observed ZD. */
            sz = Math.Sqrt(xaeo * xaeo + yaeo * yaeo);
            zdo = Math.Atan2(sz, zaeo);

            /*
            ** Refraction
            ** ----------
            */

            /* Fast algorithm using two constant model. */
            refa = astrom.refa;
            refb = astrom.refb;
            tz = sz / (zaeo > SELMIN ? zaeo : SELMIN);
            dref = (refa + refb * tz * tz) * tz;
            zdt = zdo + dref;

            /* To Cartesian Az,ZD. */
            ce = Math.Sin(zdt);
            xaet = Math.Cos(az) * ce;
            yaet = Math.Sin(az) * ce;
            zaet = Math.Cos(zdt);

            /* Cartesian Az,ZD to Cartesian -HA,Dec. */
            xmhda = sphi * xaet + cphi * zaet;
            ymhda = yaet;
            zmhda = -cphi * xaet + sphi * zaet;

            /* Diurnal aberration. */
            f = (1.0 + astrom.diurab * ymhda);
            xhd = f * xmhda;
            yhd = f * (ymhda - astrom.diurab);
            zhd = f * zmhda;

            /* Polar motion. */
            sx = Math.Sin(astrom.xpl);
            cx = Math.Cos(astrom.xpl);
            sy = Math.Sin(astrom.ypl);
            cy = Math.Cos(astrom.ypl);
            v[0] = cx * xhd + sx * sy * yhd - sx * cy * zhd;
            v[1] = cy * yhd + sy * zhd;
            v[2] = sx * xhd - cx * sy * yhd + cx * cy * zhd;

            /* To spherical -HA,Dec. */
            iauC2s(v, out hma, out di);

            /* Right ascension. */
            ri = iauAnp(astrom.eral + hma);
        }


        /// <summary>
        /// Apply light deflection by a solar-system body, as part of
        /// transforming coordinate direction into natural direction.
        /// </summary>
        /// <param name="bm">mass of the gravitating body (solar masses)</param>
        /// <param name="p">double[3]  direction from observer to source (unit vector)</param>
        /// <param name="q">double[3]  direction from body to source (unit vector)</param>
        /// <param name="e">double[3]  direction from body to observer (unit vector)</param>
        /// <param name="em">distance from body to observer (au)</param>
        /// <param name="dlim">deflection limiter</param>
        /// <param name="p1">double[3]  observer to deflected source (unit vector)</param>
        public static void iauLd(double bm, in double[] p, in double[] q, in double[] e,
           double em, double dlim, ref double[] p1)
        {
            int i;
            double qdqpe, w;
            double[] qpe = new double[3];
            double[] eq = new double[3];
            double[] peq = new double[3];

            /* q . (q + e). */
            for (i = 0; i < 3; i++)
            {
                qpe[i] = q[i] + e[i];
            }
            qdqpe = iauPdp(q, qpe);

            /* 2 x G x bm / ( em x c^2 x ( q . (q + e) ) ). */
            w = bm * SRS / em / Math.Max(qdqpe, dlim);

            /* p x (e x q). */
            iauPxp(e, q, ref eq);
            iauPxp(p, eq, ref peq);

            /* Apply the deflection. */
            for (i = 0; i < 3; i++)
            {
                p1[i] = p[i] + w * peq[i];
            }
        }


        /// <summary>
        /// For a star, apply light deflection by multiple solar-system bodies,
        /// as part of transforming coordinate direction into natural direction.
        /// </summary>
        /// <param name="n">number of bodies</param>
        /// <param name="b">data for each of the n bodies</param>
        /// <param name="ob">double[3]     barycentric position of the observer (au)</param>
        /// <param name="sc">double[3]     observer to star coord direction (unit vector)</param>
        /// <param name="sn">double[3]      observer to deflected star (unit vector)</param>
        public static void iauLdn(int n, in iauLDBODY[] b, in double[] ob, in double[] sc,
            ref double[] sn)
        {
            /* Light time for 1 au (days) */
            double CR = AULT / DAYSEC;

            int i;
            double dt, em;
            double[] v = new double[3];
            double[] ev = new double[3];
            double[] e = new double[3];


            /* Star direction prior to deflection. */
            iauCp(sc, sn);

            /* Body by body. */
            for (i = 0; i < n; i++)
            {

                /* Body to observer vector at epoch of observation (au). */
                iauPmp(ob, b[i].pv[0], ref v);

                /* Minus the time since the light passed the body (days). */
                dt = iauPdp(sn, v) * CR;

                /* Neutralize if the star is "behind" the observer. */
                dt = Math.Min(dt, 0.0);

                /* Backtrack the body to the time the light was passing the body. */
                iauPpsp(v, -dt, b[i].pv[1], ref ev);

                /* Body to observer vector as magnitude and direction. */
                iauPn(ev, out em, ref e);

                /* Apply light deflection for this body. */
                iauLd(b[i].bm, sn, sn, e, em, b[i].dl, ref sn);

                /* Next body. */
            }
        }


        /// <summary>
        /// Deflection of starlight by the Sun.
        /// </summary>
        /// <param name="p">double[3]  direction from observer to star (unit vector)</param>
        /// <param name="e">double[3]  direction from Sun to observer (unit vector)</param>
        /// <param name="em">distance from Sun to observer (au)</param>
        /// <param name="p1">double[3]  observer to deflected star (unit vector)</param>
        public static void iauLdsun(in double[] p, in double[] e, double em, ref double[] p1)
        {
            double em2, dlim;


            /* Deflection limiter (smaller for distant observers). */
            em2 = em * em;
            if (em2 < 1.0) em2 = 1.0;
            dlim = 1e-6 / (em2 > 1.0 ? em2 : 1.0);

            /* Apply the deflection. */
            iauLd(1.0, p, p, e, em, dlim, ref p1);
        }


        /// <summary>
        /// Proper motion and parallax.
        /// </summary>
        /// <param name="rc">ICRS RA at catalog epoch (radians)</param>
        /// <param name="dc">ICRS Dec at catalog epoch (radians)</param>
        /// <param name="pr">RA proper motion (radians/year)</param>
        /// <param name="pd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
        /// <param name="pmt">proper motion time interval (SSB, Julian years)</param>
        /// <param name="pob">double[3]  SSB to observer vector (au)</param>
        /// <param name="pco">double[3]  coordinate direction (BCRS unit vector)</param>
        public static void iauPmpx(double rc, double dc, double pr, double pd,
             double px, double rv, double pmt, in double[] pob,
             out double[] pco)
        {
            pco = new double[3];
            /* Km/s to au/year */
            double VF = DAYSEC * DJM / DAU;

            /* Light time for 1 au, Julian years */
            double AULTY = AULT / DAYSEC / DJY;

            int i;
            double sr, cr, sd, cd, x, y, z, dt, pxr, w, pdz;
            double[] p = new double[3];
            double[] pm = new double[3];

            /* Spherical coordinates to unit vector (and useful functions). */
            sr = Math.Sin(rc);
            cr = Math.Cos(rc);
            sd = Math.Sin(dc);
            cd = Math.Cos(dc);
            p[0] = x = cr * cd;
            p[1] = y = sr * cd;
            p[2] = z = sd;

            /* Proper motion time interval (y) including Roemer effect. */
            dt = pmt + iauPdp(p, pob) * AULTY;

            /* Space motion (radians per year). */
            pxr = px * DAS2R;
            w = VF * rv * pxr;
            pdz = pd * z;
            pm[0] = -pr * y - pdz * cr + w * x;
            pm[1] = pr * x - pdz * sr + w * y;
            pm[2] = pd * cd + w * z;

            /* Coordinate direction of star (unit vector, BCRS). */
            for (i = 0; i < 3; i++)
            {
                p[i] += dt * pm[i] - pxr * pob[i];
            }
            iauPn(p, out w, ref pco);
        }


        /// <summary>
        /// Star proper motion:  update star catalog data for space motion, with
        /// special handling to handle the zero parallax case.
        /// </summary>
        /// <param name="ra1">right ascension (radians), before</param>
        /// <param name="dec1">declination (radians), before</param>
        /// <param name="pmr1">RA proper motion (radians/year), before</param>
        /// <param name="pmd1">Dec proper motion (radians/year), before</param>
        /// <param name="px1">parallax (arcseconds), before</param>
        /// <param name="rv1">radial velocity (km/s, +ve = receding), before</param>
        /// <param name="ep1a">"before" epoch, part A</param>
        /// <param name="ep1b">"before" epoch, part B</param>
        /// <param name="ep2a">"after" epoch, part A</param>
        /// <param name="ep2b">"after" epoch, part B</param>
        /// <param name="ra2">right ascension (radians), after</param>
        /// <param name="dec2">declination (radians), after</param>
        /// <param name="pmr2">RA proper motion (radians/year), after</param>
        /// <param name="pmd2">Dec proper motion (radians/year), after</param>
        /// <param name="px2">parallax (arcseconds), after</param>
        /// <param name="rv2">radial velocity (km/s, +ve = receding), after</param>
        /// <returns>status:
        /// <para>-1 = system error (should not occur)</para>
        /// <para> 0 = no warnings or errors</para>
        /// <para> 1 = distance overridden</para>
        /// <para> 2 = excessive velocity</para>
        /// <para> 4 = solution didn't converge</para>
        /// <para>else = binary logical OR of the above warnings</para>
        /// </returns>
        public static int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
              double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              out double ra2, out double dec2, out double pmr2, out double pmd2,
              out double px2, out double rv2)
        {
            /* Minimum allowed parallax (arcsec) */
            double PXMIN = 5e-7;

            /* Factor giving maximum allowed transverse speed of about 1% c */
            double F = 326.0;

            int jpx, j;
            double pm, px1a;


            /* Proper motion in one year (radians). */
            pm = iauSeps(ra1, dec1, ra1 + pmr1, dec1 + pmd1);

            /* Override the parallax to reduce the chances of a warning status. */
            jpx = 0;
            px1a = px1;
            pm *= F;
            if (px1a < pm) { jpx = 1; px1a = pm; }
            if (px1a < PXMIN) { jpx = 1; px1a = PXMIN; }

            /* Carry out the transformation using the modified parallax. */
            j = iauStarpm(ra1, dec1, pmr1, pmd1, px1a, rv1,
                          ep1a, ep1b, ep2a, ep2b,
                          out ra2, out dec2, out pmr2, out pmd2, out px2, out rv2);

            /* Revise and return the status. */
            if ((j % 2) == 0) j += jpx;
            return j;
        }


        /// <summary>
        /// Position and velocity of a terrestrial observing station.
        /// </summary>
        /// <param name="elong">longitude (radians, east +ve)</param>
        /// <param name="phi">latitude (geodetic, radians)</param>
        /// <param name="hm">height above ref. ellipsoid (geodetic, m)</param>
        /// <param name="xp">coordinates of the pole (radians)</param>
        /// <param name="yp">coordinates of the pole (radians)</param>
        /// <param name="sp">the TIO locator s' (radians)</param>
        /// <param name="theta">Earth rotation angle (radians)</param>
        /// <param name="pv">double[2][3] position/velocity vector (m, m/s, CIRS)</param>
        public static void iauPvtob(double elong, double phi, double hm,
              double xp, double yp, double sp, double theta,
              out double[][] pv)
        {
            pv = new double[][] { new double[3], new double[3] };

            /* Earth rotation rate in radians per UT1 second */
            double OM = 1.00273781191135448 * D2PI / DAYSEC;

            double x, y, z, s, c;
            double[] xyzm;
            double[] xyz = new double[3];
            double[][] rpm;


            /* Geodetic to geocentric transformation (WGS84). */
            iauGd2gc(ReferenceEllipsoids.WGS84, elong, phi, hm, out xyzm);

            /* Polar motion and TIO position. */
            iauPom00(xp, yp, sp, out rpm);
            iauTrxp(rpm, xyzm, ref xyz);
            x = xyz[0];
            y = xyz[1];
            z = xyz[2];

            /* Functions of ERA. */
            s = Math.Sin(theta);
            c = Math.Cos(theta);

            /* Position. */
            pv[0][0] = c * x - s * y;
            pv[0][1] = s * x + c * y;
            pv[0][2] = z;

            /* Velocity. */
            pv[1][0] = OM * (-s * x - c * y);
            pv[1][1] = OM * (c * x - s * y);
            pv[1][2] = 0.0;
        }


        /// <summary>
        /// Determine the constants A and B in the atmospheric refraction model
        /// dZ = A tan Z + B tan^3 Z.
        /// </summary>
        /// <param name="phpa">pressure at the observer (hPa = millibar)</param>
        /// <param name="tc">ambient temperature at the observer (deg C)</param>
        /// <param name="rh">relative humidity at the observer (range 0-1)</param>
        /// <param name="wl">wavelength (micrometers)</param>
        /// <param name="refa">tan Z coefficient (radians)</param>
        /// <param name="refb">tan^3 Z coefficient (radians)</param>
        public static void iauRefco(double phpa, double tc, double rh, double wl,
              out double refa, out double refb)
        {
            bool optic;
            double p, t, r, w, ps, pw, tk, wlsq, gamma, beta;


            /* Decide whether optical/IR or radio case:  switch at 100 microns. */
            optic = (wl <= 100.0);

            /* Restrict parameters to safe values. */
            t = Math.Max(tc, -150.0);
            t = Math.Min(t, 200.0);
            p = Math.Max(phpa, 0.0);
            p = Math.Min(p, 10000.0);
            r = Math.Max(rh, 0.0);
            r = Math.Min(r, 1.0);
            w = Math.Max(wl, 0.1);
            w = Math.Min(w, 1e6);

            /* Water vapour pressure at the observer. */
            if (p > 0.0)
            {
                ps = Math.Pow(10.0, (0.7859 + 0.03477 * t) /
                                    (1.0 + 0.00412 * t)) *
                           (1.0 + p * (4.5e-6 + 6e-10 * t * t));
                pw = r * ps / (1.0 - (1.0 - r) * ps / p);
            }
            else
            {
                pw = 0.0;
            }

            /* Refractive index minus 1 at the observer. */
            tk = t + 273.15;
            if (optic)
            {
                wlsq = w * w;
                gamma = ((77.53484e-6 +
                           (4.39108e-7 + 3.666e-9 / wlsq) / wlsq) * p
                              - 11.2684e-6 * pw) / tk;
            }
            else
            {
                gamma = (77.6890e-6 * p - (6.3938e-6 - 0.375463 / tk) * pw) / tk;
            }

            /* Formula for beta from Stone, with empirical adjustments. */
            beta = 4.4474e-6 * tk;
            if (!optic) beta -= 0.0074 * pw * beta;

            /* Refraction constants from Green. */
            refa = gamma * (1.0 - beta);
            refb = -gamma * (beta - gamma / 2.0);
        }


    }
}
