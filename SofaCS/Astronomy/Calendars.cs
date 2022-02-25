using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Gregorian Calendar to Julian Date.
        /// </summary>
        /// <param name="iy">year</param>
        /// <param name="im">month</param>
        /// <param name="id">day in Gregorian calendar</param>
        /// <param name="djm0">MJD zero-point: always 2400000.5</param>
        /// <param name="djm">Modified Julian Date for 0 hrs</param>
        /// <returns>status:
        /// <para>0 = OK</para>
        /// <para>-1 = bad year   (Note 3: JD not computed)</para>
        /// <para>-2 = bad month  (JD not computed)</para>
        /// <para>-3 = bad day    (JD computed)</para>
        /// </returns>
        public static int iauCal2jd(int iy, int im, int id, out double djm0, out double djm)
        {
            int j, ly, my;
            long iypmy;
            djm0 = 0;
            djm = 0;

            /* Earliest year allowed (4800BC) */
            int IYMIN = -4799;

            /* Month lengths in days */
            int[] mtab = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };


            /* Preset status. */
            j = 0;

            /* Validate year and month. */
            if (iy < IYMIN) return -1;
            if (im < 1 || im > 12) return -2;

            /* If February in a leap year, 1, otherwise 0. */
            ly = ((im == 2) && !((iy % 4) > 0) && (((iy % 100) > 0) || !((iy % 400) > 0))) ? 1 : 0;

            /* Validate day, taking into account leap years. */
            if ((id < 1) || (id > (mtab[im - 1] + ly))) j = -3;

            /* Return result. */
            my = (im - 14) / 12;
            iypmy = (long)(iy + my);
            djm0 = DJM0;
            djm = (double)((1461L * (iypmy + 4800L)) / 4L
                          + (367L * (long)(im - 2 - 12 * my)) / 12L
                          - (3L * ((iypmy + 4900L) / 100L)) / 4L
                          + (long)id - 2432076L);

            /* Return status. */
            return j;
        }


        /// <summary>
        /// Julian Date to Besselian Epoch.
        /// </summary>
        /// <param name="dj1">Julian Date</param>
        /// <param name="dj2">Julian Date</param>
        /// <returns>Besselian Epoch.</returns>
        public static double iauEpb(double dj1, double dj2)
        {
            /* J2000.0-B1900.0 (2415019.81352) in days */
            double D1900 = 36524.68648;

            return 1900.0 + ((dj1 - DJ00) + (dj2 + D1900)) / DTY;
        }


        /// <summary>
        /// Besselian Epoch to Julian Date.
        /// </summary>
        /// <param name="epb">Besselian Epoch (e.g. 1957.3)</param>
        /// <param name="djm0">MJD zero-point: always 2400000.5</param>
        /// <param name="djm">Modified Julian Date</param>
        public static void iauEpb2jd(double epb, out double djm0, out double djm)
        {
            djm0 = DJM0;
            djm = 15019.81352 + (epb - 1900.0) * DTY;
        }


        /// <summary>
        /// Julian Date to Julian Epoch.
        /// </summary>
        /// <param name="dj1">Julian Date</param>
        /// <param name="dj2">Julian Date</param>
        /// <returns>Julian Epoch</returns>
        public static double iauEpj(double dj1, double dj2)
        {
            double epj;


            epj = 2000.0 + ((dj1 - DJ00) + dj2) / DJY;

            return epj;
        }


        /// <summary>
        /// Julian Epoch to Julian Date.
        /// </summary>
        /// <param name="epj">Julian Epoch (e.g. 1996.8)</param>
        /// <param name="djm0">MJD zero-point: always 2400000.5</param>
        /// <param name="djm">Modified Julian Date</param>
        public static void iauEpj2jd(double epj, out double djm0, out double djm)
        {
            djm0 = DJM0;
            djm = DJM00 + (epj - 2000.0) * 365.25;
        }


        /// <summary>
        /// Julian Date to Gregorian year, month, day, and fraction of a day.
        /// </summary>
        /// <param name="dj1">Julian Date</param>
        /// <param name="dj2">Julian Date</param>
        /// <param name="iy">year</param>
        /// <param name="im">month</param>
        /// <param name="id">day</param>
        /// <param name="fd">fraction of day</param>
        /// <returns>status:
        /// <para>0 = OK</para>
        /// <para>-1 = unacceptable date</para>
        /// </returns>
        public static int iauJd2cal(double dj1, double dj2,
              out int iy, out int im, out int id, out double fd)
        {
            /* Minimum and maximum allowed JD */
            double DJMIN = -68569.5;
            double DJMAX = 1e9;

            long jd, i, l, n, k;
            double dj, f1, f2, d, s, cs, x, t, f;
            double[] v = new double[2];

            iy = 0;
            im = 0;
            id = 0;
            fd = 0;

            /* Verify date is acceptable. */
            dj = dj1 + dj2;
            if (dj < DJMIN || dj > DJMAX) return -1;

            /* Separate day and fraction (where -0.5 <= fraction < 0.5). */
            d = dnint(dj1);
            f1 = dj1 - d;
            jd = (long)d;
            d = dnint(dj2);
            f2 = dj2 - d;
            jd += (long)d;

            /* Compute f1+f2+0.5 using compensated summation (Klein 2006). */
            s = 0.5;
            cs = 0.0;
            v[0] = f1;
            v[1] = f2;
            for (i = 0; i < 2; i++)
            {
                x = v[i];
                t = s + x;
                cs += Math.Abs(s) >= Math.Abs(x) ? (s - t) + x : (x - t) + s;
                s = t;
                if (s >= 1.0)
                {
                    jd++;
                    s -= 1.0;
                }
            }
            f = s + cs;
            cs = f - s;

            /* Deal with negative f. */
            if (f < 0.0)
            {

                /* Compensated summation: assume that |s| <= 1.0. */
                f = s + 1.0;
                cs += (1.0 - f) + s;
                s = f;
                f = s + cs;
                cs = f - s;
                jd--;
            }

            /* Deal with f that is 1.0 or more (when rounded to double). */
            if ((f - 1.0) >= -DBL_EPSILON / 4.0)
            {

                /* Compensated summation: assume that |s| <= 1.0. */
                t = s - 1.0;
                cs += (s - t) - 1.0;
                s = t;
                f = s + cs;
                if (-DBL_EPSILON / 2.0 < f)
                {
                    jd++;
                    f = Math.Max(f, 0.0);
                }
            }

            /* Express day in Gregorian calendar. */
            l = jd + 68569L;
            n = (4L * l) / 146097L;
            l -= (146097L * n + 3L) / 4L;
            i = (4000L * (l + 1L)) / 1461001L;
            l -= (1461L * i) / 4L - 31L;
            k = (80L * l) / 2447L;
            id = (int)(l - (2447L * k) / 80L);
            l = k / 11L;
            im = (int)(k + 2L - 12L * l);
            iy = (int)(100L * (n - 49L) + i + l);
            fd = f;

            /* Success. */
            return 0;
        }


        /// <summary>
        /// Julian Date to Gregorian Calendar, expressed in a form convenient
        /// for formatting messages:  rounded to a specified precision.
        /// </summary>
        /// <param name="ndp">number of decimal places of days in fraction</param>
        /// <param name="dj1">dj1+dj2 = Julian Date</param>
        /// <param name="dj2">dj1+dj2 = Julian Date</param>
        /// <param name="iymdf">int[4]   year, month, day, fraction in Gregorian calendar</param>
        /// <returns>status:
        /// <para>-1 = date out of range</para>
        /// <para>0 = OK</para>
        /// <para>+1 = NDP not 0-9 (interpreted as 0)</para>
        /// </returns>
        public static int iauJdcalf(int ndp, double dj1, double dj2, out int[] iymdf)
        {
            int j, js;
            double denom, d1, d2, f1, f2, d, djd, f, rf;

            iymdf = new int[4];

            /* Denominator of fraction (e.g. 100 for 2 decimal places). */
            if ((ndp >= 0) && (ndp <= 9))
            {
                j = 0;
                denom = Math.Pow(10.0, ndp);
            }
            else
            {
                j = 1;
                denom = 1.0;
            }

            /* Copy the date, big then small. */
            if (Math.Abs(dj1) >= Math.Abs(dj2))
            {
                d1 = dj1;
                d2 = dj2;
            }
            else
            {
                d1 = dj2;
                d2 = dj1;
            }

            /* Realign to midnight (without rounding error). */
            d1 -= 0.5;

            /* Separate day and fraction (as precisely as possible). */
            d = dnint(d1);
            f1 = d1 - d;
            djd = d;
            d = dnint(d2);
            f2 = d2 - d;
            djd += d;
            d = dnint(f1 + f2);
            f = (f1 - d) + f2;
            if (f < 0.0)
            {
                f += 1.0;
                d -= 1.0;
            }
            djd += d;

            /* Round the total fraction to the specified number of places. */
            rf = dnint(f * denom) / denom;

            /* Re-align to noon. */
            djd += 0.5;

            /* Convert to Gregorian calendar. */
            js = iauJd2cal(djd, rf, out iymdf[0], out iymdf[1], out iymdf[2], out f);
            if (js == 0)
            {
                iymdf[3] = (int)dnint(f * denom);
            }
            else
            {
                j = js;
            }

            /* Return the status. */
            return j;
        }


    }
}
