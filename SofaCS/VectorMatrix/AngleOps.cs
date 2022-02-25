using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Decompose radians into degrees, arcminutes, arcseconds, fraction.
        /// </summary>
        /// <param name="ndp">resolution</param>
        /// <param name="angle">angle in radians</param>
        /// <param name="sign">'+' or '-'</param>
        /// <param name="idmsf">int[4]  degrees, arcminutes, arcseconds, fraction</param>
        public static void iauA2af(int ndp, double angle, out char sign, out int[] idmsf)
        {
            /* Hours to degrees * radians to turns */
            double F = 15.0 / D2PI;

            /* Scale then use days to h,m,s function. */
            iauD2tf(ndp, angle * F, out sign, out idmsf);
        }


        /// <summary>
        /// Decompose radians into hours, minutes, seconds, fraction.
        /// </summary>
        /// <param name="ndp">resolution</param>
        /// <param name="angle">angle in radians</param>
        /// <param name="sign">'+' or '-'</param>
        /// <param name="ihmsf">int[4]  hours, minutes, seconds, fraction</param>
        public static void iauA2tf(int ndp, double angle, out char sign, out int[] ihmsf)
        {
            /* Scale then use days to h,m,s function. */
            iauD2tf(ndp, angle / D2PI, out sign, out ihmsf);
        }


        /// <summary>
        /// Convert degrees, arcminutes, arcseconds to radians.
        /// </summary>
        /// <param name="s">sign:  '-' = negative, otherwise positive</param>
        /// <param name="ideg">degrees</param>
        /// <param name="iamin">arcminutes</param>
        /// <param name="asec">arcseconds</param>
        /// <param name="rad">angle in radians</param>
        /// <returns>status:  0 = OK
        /// <para>1 = ideg outside range 0-359</para>
        /// <para>2 = iamin outside range 0-59</para>
        /// <para>3 = asec outside range 0-59.999...</para>
        /// </returns>
        public static int iauAf2a(char s, int ideg, int iamin, double asec, out double rad)
        {
            /* Compute the interval. */
            rad = (s == '-' ? -1.0 : 1.0) *
                    (60.0 * (60.0 * ((double)Math.Abs(ideg)) +
                                      ((double)Math.Abs(iamin))) +
                                                 Math.Abs(asec)) * DAS2R;

            /* Validate arguments and return status. */
            if (ideg < 0 || ideg > 359) return 1;
            if (iamin < 0 || iamin > 59) return 2;
            if (asec < 0.0 || asec >= 60.0) return 3;
            return 0;
        }


        /// <summary>
        /// Normalize angle into the range 0 <= a < 2pi.
        /// </summary>
        /// <param name="a">angle (radians)</param>
        /// <returns>angle in range 0-2pi</returns>
        public static double iauAnp(double a)
        {
            double w;

            w = a % D2PI;
            if (w < 0) w += D2PI;

            return w;
        }


        /// <summary>
        /// Normalize angle into the range -pi <= a < +pi.
        /// </summary>
        /// <param name="a">angle (radians)</param>
        /// <returns>angle in range +/-pi</returns>
        public static double iauAnpm(double a)
        {
            double w;

            w = a % D2PI;
            if (Math.Abs(w) >= DPI) w -= dsign(D2PI, a);

            return w;
        }


        /// <summary>
        /// Decompose days to hours, minutes, seconds, fraction.
        /// </summary>
        /// <param name="ndp">resolution</param>
        /// <param name="days">interval in days</param>
        /// <param name="sign">'+' or '-'</param>
        /// <param name="ihmsf">int[4]  hours, minutes, seconds, fraction</param>
        public static void iauD2tf(int ndp, double days, out char sign, out int[] ihmsf)
        {
            int nrs, n;
            double rs, rm, rh, a, w, ah, am, As, af;

            ihmsf = new int[4];

            /* Handle sign. */
            sign = (char)((days >= 0.0) ? '+' : '-');

            /* Interval in seconds. */
            a = DAYSEC * Math.Abs(days);

            /* Pre-round if resolution coarser than 1s (then pretend ndp=1). */
            if (ndp < 0)
            {
                nrs = 1;
                for (n = 1; n <= -ndp; n++)
                {
                    nrs *= (n == 2 || n == 4) ? 6 : 10;
                }
                rs = (double)nrs;
                w = a / rs;
                a = rs * dnint(w);
            }

            /* Express the unit of each field in resolution units. */
            nrs = 1;
            for (n = 1; n <= ndp; n++)
            {
                nrs *= 10;
            }
            rs = (double)nrs;
            rm = rs * 60.0;
            rh = rm * 60.0;

            /* Round the interval and express in resolution units. */
            a = dnint(rs * a);

            /* Break into fields. */
            ah = a / rh;
            ah = dint(ah);
            a -= ah * rh;
            am = a / rm;
            am = dint(am);
            a -= am * rm;
            As = a / rs;
            As = dint(As);
            af = a - As * rs;

            /* Return results. */
            ihmsf[0] = (int)ah;
            ihmsf[1] = (int)am;
            ihmsf[2] = (int)As;
            ihmsf[3] = (int)af;
        }


        /// <summary>
        /// Convert hours, minutes, seconds to radians.
        /// </summary>
        /// <param name="s">sign:  '-' = negative, otherwise positive</param>
        /// <param name="ihour">hours</param>
        /// <param name="imin">minutes</param>
        /// <param name="sec">seconds</param>
        /// <param name="rad">angle in radians</param>
        /// <returns>status:
        /// <para>0 = OK</para>
        /// <para>1 = ihour outside range 0-23</para>
        /// <para>2 = imin outside range 0-59</para>
        /// <para>3 = sec outside range 0-59.999...</para>
        /// </returns>
        public static int iauTf2a(char s, int ihour, int imin, double sec, out double rad)
        {
            /* Compute the interval. */
            rad = (s == '-' ? -1.0 : 1.0) *
                    (60.0 * (60.0 * ((double)Math.Abs(ihour)) +
                                      ((double)Math.Abs(imin))) +
                                                 Math.Abs(sec)) * DS2R;

            /* Validate arguments and return status. */
            if (ihour < 0 || ihour > 23) return 1;
            if (imin < 0 || imin > 59) return 2;
            if (sec < 0.0 || sec >= 60.0) return 3;
            return 0;
        }


        /// <summary>
        /// Convert hours, minutes, seconds to days.
        /// </summary>
        /// <param name="s">sign:  '-' = negative, otherwise positive</param>
        /// <param name="ihour">hours</param>
        /// <param name="imin">minutes</param>
        /// <param name="sec">seconds</param>
        /// <param name="days">interval in days</param>
        /// <returns>status:
        /// <para>0 = OK</para>
        /// <para>1 = ihour outside range 0-23</para>
        /// <para>2 = imin outside range 0-59</para>
        /// <para>3 = sec outside range 0-59.999...</para>
        /// </returns>
        public static int iauTf2d(char s, int ihour, int imin, double sec, out double days)
        {
            /* Compute the interval. */
            days = (s == '-' ? -1.0 : 1.0) *
                     (60.0 * (60.0 * ((double)Math.Abs(ihour)) +
                                       ((double)Math.Abs(imin))) +
                                                  Math.Abs(sec)) / DAYSEC;

            /* Validate arguments and return status. */
            if (ihour < 0 || ihour > 23) return 1;
            if (imin < 0 || imin > 59) return 2;
            if (sec < 0.0 || sec >= 60.0) return 3;
            return 0;
        }


    }
}
