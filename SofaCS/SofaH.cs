using System;


namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Pi
        /// </summary>
        public static readonly double DPI = 3.141592653589793238462643;

        /// <summary>
        /// 2Pi
        /// </summary>
        public static readonly double D2PI = 6.283185307179586476925287;

        /// <summary>
        /// Radians to degrees
        /// </summary>
        public static readonly double DR2D = 57.29577951308232087679815;

        /// <summary>
        /// Degrees to radians
        /// </summary>
        public static readonly double DD2R = 1.745329251994329576923691e-2;

        /// <summary>
        /// Radians to arcseconds
        /// </summary>
        public static readonly double DR2AS = 206264.8062470963551564734;

        /// <summary>
        /// Arcseconds to radians
        /// </summary>
        public static readonly double DAS2R = 4.848136811095359935899141e-6;

        /// <summary>
        /// Seconds of time to radians
        /// </summary>
        public static readonly double DS2R = 7.272205216643039903848712e-5;

        /// <summary>
        /// Arcseconds in a full circle
        /// </summary>
        public static readonly double TURNAS = 1296000.0;

        /// <summary>
        /// Milliarcseconds to radians
        /// </summary>
        public static readonly double DMAS2R = DAS2R / 1e3;

        /// <summary>
        /// Length of tropical year B1900 (days)
        /// </summary>
        public static readonly double DTY = 365.242198781;

        /// <summary>
        /// Seconds per day.
        /// </summary>
        public static readonly double DAYSEC = 86400.0;

        /// <summary>
        /// Days per Julian year
        /// </summary>
        public static readonly double DJY = 365.25;

        /// <summary>
        /// Days per Julian century
        /// </summary>
        public static readonly double DJC = 36525.0;

        /// <summary>
        /// Days per Julian millennium
        /// </summary>
        public static readonly double DJM = 365250.0;

        /// <summary>
        /// Reference epoch (J2000.0), Julian Date
        /// </summary>
        public static readonly double DJ00 = 2451545.0;

        /// <summary>
        /// Julian Date of Modified Julian Date zero
        /// </summary>
        public static readonly double DJM0 = 2400000.5;

        /// <summary>
        /// Reference epoch (J2000.0), Modified Julian Date
        /// </summary>
        public static readonly double DJM00 = 51544.5;

        /// <summary>
        ///1977 Jan 1.0 as MJD 
        /// </summary>
        public static readonly double DJM77 = 43144.0;

        /// <summary>
        /// TT minus TAI (s)
        /// </summary>
        public static readonly double TTMTAI = 32.184;

        /// <summary>
        /// Astronomical unit (m, IAU 2012)
        /// </summary>
        public static readonly double DAU = 149597870.7e3;

        /// <summary>
        /// Speed of light (m/s)
        /// </summary>
        public static readonly double CMPS = 299792458.0;

        /// <summary>
        /// Light time for 1 au (s)
        /// </summary>
        public static readonly double AULT = DAU / CMPS;

        /// <summary>
        /// Speed of light (au per day)
        /// </summary>
        public static readonly double DC = DAYSEC / AULT;

        /// <summary>
        /// L_G = 1 - d(TT)/d(TCG)
        /// </summary>
        public static readonly double ELG = 6.969290134e-10;

        /// <summary>
        /// L_B = 1 - d(TDB)/d(TCB)
        /// </summary>
        public static readonly double ELB = (1.550519768e-8);

        /// <summary>
        /// TDB (s) at TAI 1977/1/1.0
        /// </summary>
        public static readonly double TDB0 = (-6.55e-5);

        /// <summary>
        /// Schwarzschild radius of the Sun (au)
        /// = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11
        /// </summary>
        public static readonly double SRS = 1.97412574336e-8;

        /// <summary>
        /// smallest such that 1.0+DBL_EPSILON != 1.0
        /// </summary>
        public static readonly double DBL_EPSILON = 2.2204460492503131e-016;

        /// <summary>
        /// truncate to nearest whole number towards zero (double)
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        public static double dint(double d)
        {
            return d < 0 ? Math.Ceiling(d) : Math.Floor(d);
        }

        /// <summary>
        /// round to nearest whole number (double)
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        public static double dnint(double d)
        {
            return Math.Abs(d) < 0.5 ? 0 : (d < 0 ? Math.Ceiling(d - 0.5) : Math.Floor(d + 0.5));
        }


        /// <summary>
        /// magnitude of A with sign of B (double)
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static double dsign(double A, double B)
        {
            return B < 0 ? -Math.Abs(A) : Math.Abs(A);
        }


        public enum ReferenceEllipsoids
        {
            WGS84=1,
            GRS80=2,
            WGS72=3
        }
    }

    /// <summary>
    /// Star-independent astrometry parameters
    /// </summary>
    public class iauASTROM
    {
        /// <summary>
        /// PM time interval (SSB, Julian years)
        /// </summary>
        public double pmt;

        /// <summary>
        /// SSB to observer (vector, au)
        /// </summary>
        public double[] eb = new double[3];

        /// <summary>
        /// Sun to observer (unit vector)
        /// </summary>
        public double[] eh = new double[3];

        /// <summary>
        /// distance from Sun to observer (au)
        /// </summary>
        public double em;

        /// <summary>
        /// barycentric observer velocity (vector, c)
        /// </summary>
        public double[] v = new double[3];

        /// <summary>
        /// sqrt(1-|v|^2): reciprocal of Lorenz factor
        /// </summary>
        public double bm1;

        /// <summary>
        /// bias-precession-nutation matrix
        /// </summary>
        public double[][] bpn = { new double[3], new double[3], new double[3] };

        /// <summary>
        /// longitude + s' + dERA(DUT) (radians)
        /// </summary>
        public double along;

        /// <summary>
        /// geodetic latitude (radians)
        /// </summary>
        public double phi;

        /// <summary>
        /// olar motion xp wrt local meridian (radians)
        /// </summary>
        public double xpl;

        /// <summary>
        /// polar motion yp wrt local meridian (radians)
        /// </summary>
        public double ypl;

        /// <summary>
        /// sine of geodetic latitude
        /// </summary>
        public double sphi;

        /// <summary>
        /// cosine of geodetic latitude
        /// </summary>
        public double cphi;

        /// <summary>
        /// magnitude of diurnal aberration vector
        /// </summary>
        public double diurab;

        /// <summary>
        /// "local" Earth rotation angle (radians)
        /// </summary>
        public double eral;

        /// <summary>
        /// refraction constant A (radians)
        /// </summary>
        public double refa;

        /// <summary>
        /// refraction constant B (radians)
        /// </summary>
        public double refb;
    }


    /// <summary>
    /// Body parameters for light deflection
    /// </summary>
    public class iauLDBODY
    {
        /// <summary>
        /// mass of the body (solar masses)
        /// </summary>
        public double bm;

        /// <summary>
        /// deflection limiter (radians^2/2)
        /// </summary>
        public double dl;

        /// <summary>
        /// barycentric PV of the body (au, au/day)
        /// </summary>
        public double[][] pv = { new double[3], new double[3] };
    }
}