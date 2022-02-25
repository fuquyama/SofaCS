using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Convert star position+velocity vector to catalog coordinates.
        /// </summary>
        /// <param name="pv">double[2][3]   pv-vector (au, au/day)</param>
        /// <param name="ra">right ascension (radians)</param>
        /// <param name="dec">declination (radians) </param>
        /// <param name="pmr">RA proper motion (radians/year)</param>
        /// <param name="pmd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcsec)</param>
        /// <param name="rv">radial velocity (km/s, positive = receding)</param>
        /// <returns>status:
        /// <para> 0 = OK</para>
        /// <para> 1 = distance overridden</para>
        /// <para>-1 = superluminal speed</para>
        /// <para>-2 = null position vector</para>
        /// </returns>
        public static int iauPvstar(in double[][] pv, out double ra, out double dec,
              out double pmr, out double pmd, out double px, out double rv)
        {
            double r, vr, vt, bett, betr, d, w, del,
                   a, rad, decd, rd;
            double[] x = new double[3];
            double[] ur = new double[3];
            double[] ut = new double[3];
            double[] usr = new double[3];
            double[] ust = new double[3];

            ra = 0; dec = 0; pmr = 0; pmd = 0; px = 0; rv = 0;


            /* Isolate the radial component of the velocity (au/day, inertial). */
            iauPn(pv[0], out r, ref x);
            vr = iauPdp(x, pv[1]);
            iauSxp(vr, x, ref ur);

            /* Isolate the transverse component of the velocity (au/day, inertial). */
            iauPmp(pv[1], ur, ref ut);
            vt = iauPm(ut);

            /* Special-relativity dimensionless parameters. */
            bett = vt / DC;
            betr = vr / DC;

            /* The inertial-to-observed correction terms. */
            d = 1.0 + betr;
            w = betr * betr + bett * bett;
            if (d == 0.0 || w > 1.0) return -1;
            del = -w / (Math.Sqrt(1.0 - w) + 1.0);

            /* Apply relativistic correction factor to radial velocity component. */
            w = (betr != 0) ? (betr - del) / (betr * d) : 1.0;
            iauSxp(w, ur, ref usr);

            /* Apply relativistic correction factor to tangential velocity */
            /* component.                                                  */
            iauSxp(1.0 / d, ut, ref ust);

            /* Combine the two to obtain the observed velocity vector (au/day). */
            iauPpp(usr, ust, ref pv[1]);

            /* Cartesian to spherical. */
            iauPv2s(pv, out a, out dec, out r, out rad, out decd, out rd);
            if (r == 0.0) return -2;

            /* Return RA in range 0 to 2pi. */
            ra = iauAnp(a);

            /* Return proper motions in radians per year. */
            pmr = rad * DJY;
            pmd = decd * DJY;

            /* Return parallax in arcsec. */
            px = DR2AS / r;

            /* Return radial velocity in km/s. */
            rv = 1e-3 * rd * DAU / DAYSEC;

            /* Success. */
            return 0;
        }


        /// <summary>
        /// Convert star catalog coordinates to position+velocity vector.
        /// </summary>
        /// <param name="ra">right ascension (radians)</param>
        /// <param name="dec">declination (radians)</param>
        /// <param name="pmr">RA proper motion (radians/year)</param>
        /// <param name="pmd">Dec proper motion (radians/year)</param>
        /// <param name="px">parallax (arcseconds)</param>
        /// <param name="rv">radial velocity (km/s, positive = receding)</param>
        /// <param name="pv">double[2][3]  pv-vector (au, au/day)</param>
        /// <returns>status:
        /// <para>0 = no warnings</para>
        /// <para>1 = distance overridden</para>
        /// <para>2 = excessive speed</para>
        /// <para>4 = solution didn't converge</para>
        /// <para>else = binary logical OR of the above</para>
        /// </returns>
        public static int iauStarpv(double ra, double dec,
              double pmr, double pmd, double px, double rv,
              out double[][] pv)
        {
            /* Smallest allowed parallax */
            double PXMIN = 1e-7;

            /* Largest allowed speed (fraction of c) */
            double VMAX = 0.5;

            /* Maximum number of iterations for relativistic solution */
            int IMAX = 100;

            int i, iwarn;
            double w, r, rd, rad, decd, v,
                   vsr, vst, betst, betsr, bett, betr,
                   dd, ddel,
                   d = 0.0, del = 0.0,       /* to prevent */
                   odd = 0.0, oddel = 0.0,   /* compiler   */
                   od = 0.0, odel = 0.0;     /* warnings   */
            double[] x = new double[3];
            double[] ur = new double[3];
            double[] ut = new double[3];
            double[] usr = new double[3];
            double[] ust = new double[3];

            /* Distance (au). */
            if (px >= PXMIN)
            {
                w = px;
                iwarn = 0;
            }
            else
            {
                w = PXMIN;
                iwarn = 1;
            }
            r = DR2AS / w;

            /* Radial velocity (au/day). */
            rd = DAYSEC * rv * 1e3 / DAU;

            /* Proper motion (radian/day). */
            rad = pmr / DJY;
            decd = pmd / DJY;

            /* To pv-vector (au,au/day). */
            iauS2pv(ra, dec, r, rad, decd, rd, out pv);

            /* If excessive velocity, arbitrarily set it to zero. */
            v = iauPm(pv[1]);
            if (v / DC > VMAX)
            {
                iauZp(out pv[1]);
                iwarn += 2;
            }

            /* Isolate the radial component of the velocity (au/day). */
            iauPn(pv[0], out w, ref x);
            vsr = iauPdp(x, pv[1]);
            iauSxp(vsr, x, ref usr);

            /* Isolate the transverse component of the velocity (au/day). */
            iauPmp(pv[1], usr, ref ust);
            vst = iauPm(ust);

            /* Special-relativity dimensionless parameters. */
            betsr = vsr / DC;
            betst = vst / DC;

            /* Determine the inertial-to-observed relativistic correction terms. */
            bett = betst;
            betr = betsr;
            for (i = 0; i < IMAX; i++)
            {
                d = 1.0 + betr;
                w = betr * betr + bett * bett;
                del = -w / (Math.Sqrt(1.0 - w) + 1.0);
                betr = d * betsr + del;
                bett = d * betst;
                if (i > 0)
                {
                    dd = Math.Abs(d - od);
                    ddel = Math.Abs(del - odel);
                    if ((i > 1) && (dd >= odd) && (ddel >= oddel)) break;
                    odd = dd;
                    oddel = ddel;
                }
                od = d;
                odel = del;
            }
            if (i >= IMAX) iwarn += 4;

            /* Replace observed radial velocity with inertial value. */
            w = (betsr != 0.0) ? d + del / betsr : 1.0;
            iauSxp(w, usr, ref ur);

            /* Replace observed tangential velocity with inertial value. */
            iauSxp(d, ust, ref ut);

            /* Combine the two to obtain the inertial space velocity. */
            iauPpp(ur, ut, ref pv[1]);

            /* Return the status. */
            return iwarn;
        }


    }
}
