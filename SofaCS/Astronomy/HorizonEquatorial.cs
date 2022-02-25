using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Horizon to equatorial coordinates:  transform azimuth and altitude to hour angle and declination.
        /// <para>All the arguments are angles in radians.</para>
        /// </summary>
        /// <param name="az">azimuth</param>
        /// <param name="el">altitude (informally, elevation)</param>
        /// <param name="phi">site latitude</param>
        /// <param name="ha">hour angle (local)</param>
        /// <param name="dec">declination</param>
        public static void iauAe2hd(double az, double el, double phi,
            out double ha, out double dec)
        {
            double sa, ca, se, ce, sp, cp, x, y, z, r;


            /* Useful trig functions. */
            sa = Math.Sin(az);
            ca = Math.Cos(az);
            se = Math.Sin(el);
            ce = Math.Cos(el);
            sp = Math.Sin(phi);
            cp = Math.Cos(phi);

            /* HA,Dec unit vector. */
            x = -ca * ce * sp + se * cp;
            y = -sa * ce;
            z = ca * ce * cp + se * sp;

            /* To spherical. */
            r = Math.Sqrt(x * x + y * y);
            ha = (r != 0.0) ? Math.Atan2(y, x) : 0.0;
            dec = Math.Atan2(z, r);
        }


        /// <summary>
        /// Equatorial to horizon coordinates:  transform hour angle and
        /// declination to azimuth and altitude.
        /// </summary>
        /// <param name="ha">hour angle (local)</param>
        /// <param name="dec">declination</param>
        /// <param name="phi">site latitude</param>
        /// <param name="az">azimuth</param>
        /// <param name="el">altitude (informally, elevation)</param>
        public static void iauHd2ae(double ha, double dec, double phi,
               out double az, out double el)
        {
            double sh, ch, sd, cd, sp, cp, x, y, z, r, a;


            /* Useful trig functions. */
            sh = Math.Sin(ha);
            ch = Math.Cos(ha);
            sd = Math.Sin(dec);
            cd = Math.Cos(dec);
            sp = Math.Sin(phi);
            cp = Math.Cos(phi);

            /* Az,Alt unit vector. */
            x = -ch * cd * sp + sd * cp;
            y = -sh * cd;
            z = ch * cd * cp + sd * sp;

            /* To spherical. */
            r = Math.Sqrt(x * x + y * y);
            a = (r != 0.0) ? Math.Atan2(y, x) : 0.0;
            az = (a < 0.0) ? a + D2PI : a;
            el = Math.Atan2(z, r);
        }


        /// <summary>
        /// Parallactic angle for a given hour angle and declination.
        /// </summary>
        /// <param name="ha">hour angle</param>
        /// <param name="dec">declination</param>
        /// <param name="phi">site latitude</param>
        /// <returns>parallactic angle</returns>
        public static double iauHd2pa(double ha, double dec, double phi)
        {
            double cp, cqsz, sqsz;


            cp = Math.Cos(phi);
            sqsz = cp * Math.Sin(ha);
            cqsz = Math.Sin(phi) * Math.Cos(dec) - cp * Math.Sin(dec) * Math.Cos(ha);
            return ((sqsz != 0.0 || cqsz != 0.0) ? Math.Atan2(sqsz, cqsz) : 0.0);
        }


    }
}
