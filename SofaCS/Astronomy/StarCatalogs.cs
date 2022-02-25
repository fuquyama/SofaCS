using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Convert B1950.0 FK4 star catalog data to J2000.0 FK5.
        /// </summary>
        /// <param name="r1950">B1950.0 RA (rad)</param>
        /// <param name="d1950">B1950.0 Dec (rad)</param>
        /// <param name="dr1950">B1950.0 RA proper motions (rad/trop.yr)</param>
        /// <param name="dd1950">B1950.0 Dec proper motions (rad/trop.yr)</param>
        /// <param name="p1950">parallax (arcsec)</param>
        /// <param name="v1950">radial velocity (km/s, +ve = moving away)</param>
        /// <param name="r2000">J2000.0 RA (rad)</param>
        /// <param name="d2000">J2000.0 Dec (rad)</param>
        /// <param name="dr2000">J2000.0 RA proper motions (rad/Jul.yr)</param>
        /// <param name="dd2000">J2000.0 Dec proper motions (rad/Jul.yr)</param>
        /// <param name="p2000">parallax (arcsec)</param>
        /// <param name="v2000">radial velocity (km/s, +ve = moving away)</param>
        public static void iauFk425(double r1950, double d1950,
              double dr1950, double dd1950,
              double p1950, double v1950,
              out double r2000, out double d2000,
              out double dr2000, out double dd2000,
              out double p2000, out double v2000)
        {
            /* Radians per year to arcsec per century */
            double PMF = 100.0 * DR2AS;

            /* Small number to avoid arithmetic problems */
            double TINY = 1e-30;

            /* Miscellaneous */
            double r, d, ur, ud, px, rv, pxvf, w, rd;
            int i, j, k, l;

            /* Pv-vectors */
            double[][] r0;
            double[][] pv1 = { new double[3], new double[3] };
            double[][] pv2 = { new double[3], new double[3] };

            /*
            ** CANONICAL CONSTANTS (Seidelmann 1992)
            */

            /* Km per sec to AU per tropical century */
            /* = 86400 * 36524.2198782 / 149597870.7 */
            double VF = 21.095;

            /* Constant pv-vector (cf. Seidelmann 3.591-2, vectors A and Adot) */
            double[][] a =  {
               new double[] { -1.62557e-6, -0.31919e-6, -0.13843e-6 },
                     new double[] { +1.245e-3,   -1.580e-3,   -0.659e-3   }
            };

            /* 3x2 matrix of pv-vectors (cf. Seidelmann 3.591-4, matrix M) */
            double[][][][] em = {

                new double[][][]{
                    new double[][]{
                        new double[] { +0.9999256782,     -0.0111820611,     -0.0048579477     },
                        new double[]{ +0.00000242395018, -0.00000002710663, -0.00000001177656 }
                    },

                    new double[][] {
                        new double[] { +0.0111820610,     +0.9999374784,     -0.0000271765     },
                        new double[] { +0.00000002710663, +0.00000242397878, -0.00000000006587 }
                    },

                    new double[][]{
                        new double[] { +0.0048579479,     -0.0000271474,     +0.9999881997,    },
                        new double[] { +0.00000001177656, -0.00000000006582, +0.00000242410173 }
                    }
                },

                new double[][][]{
                    new double[][] {
                        new double[] { -0.000551,         -0.238565,         +0.435739        },
                        new double[] { +0.99994704,       -0.01118251,       -0.00485767       }
                    },

                    new double[][]{
                        new double[] { +0.238514,         -0.002667,         -0.008541        },
                        new double[]{ +0.01118251,       +0.99995883,       -0.00002718       }
                    },

                    new double[][]{
                        new double[]  { -0.435623,         +0.012254,         +0.002117         },
                        new double[] { +0.00485767,       -0.00002714,       +1.00000956       }
                    }
                }

            };

            /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* The FK4 data (units radians and arcsec per tropical century). */
            r = r1950;
            d = d1950;
            ur = dr1950 * PMF;
            ud = dd1950 * PMF;
            px = p1950;
            rv = v1950;

            /* Express as a pv-vector. */
            pxvf = px * VF;
            w = rv * pxvf;
            iauS2pv(r, d, 1.0, ur, ud, w, out r0);

            /* Allow for E-terms (cf. Seidelmann 3.591-2). */
            iauPvmpv(r0, a, ref pv1);
            iauSxp(iauPdp(r0[0], a[0]), r0[0], ref pv2[0]);
            iauSxp(iauPdp(r0[0], a[1]), r0[0], ref pv2[1]);
            iauPvppv(pv1, pv2, ref pv1);

            /* Convert pv-vector to Fricke system (cf. Seidelmann 3.591-3). */
            for (i = 0; i < 2; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    w = 0.0;
                    for (k = 0; k < 2; k++)
                    {
                        for (l = 0; l < 3; l++)
                        {
                            w += em[i][j][k][l] * pv1[k][l];
                        }
                    }
                    pv2[i][j] = w;
                }
            }

            /* Revert to catalog form. */
            iauPv2s(pv2, out r, out d, out w, out ur, out ud, out rd);
            if (px > TINY)
            {
                rv = rd / pxvf;
                px = px / w;
            }

            /* Return the results. */
            r2000 = iauAnp(r);
            d2000 = d;
            dr2000 = ur / PMF;
            dd2000 = ud / PMF;
            v2000 = rv;
            p2000 = px;
        }


        /// <summary>
        /// Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
        /// proper motion in the FK5 system.
        /// </summary>
        /// <param name="r1950">B1950.0 FK4 RA at epoch (rad)</param>
        /// <param name="d1950">B1950.0 FK4 Dec at epoch (rad)</param>
        /// <param name="bepoch">Besselian epoch (e.g. 1979.3)</param>
        /// <param name="r2000">J2000.0 FK5 RA (rad)</param>
        /// <param name="d2000">J2000.0 FK5 Dec (rad)</param>
        public static void iauFk45z(double r1950, double d1950, double bepoch,
              out double r2000, out double d2000)
        {
            /* Radians per year to arcsec per century */
            double PMF = 100.0 * DR2AS;

            /* Position and position+velocity vectors */
            double[] r0;
            double[] p = new double[3];
            double[][] pv = { new double[3], new double[3] };

            /* Miscellaneous */
            double w, djm0, djm;
            int i, j, k;

            /*
            ** CANONICAL CONSTANTS (Seidelmann 1992)
            */

            /* Vectors A and Adot (Seidelmann 3.591-2) */
            double[] a = { -1.62557e-6, -0.31919e-6, -0.13843e-6 };
            double[] ad = { +1.245e-3, -1.580e-3, -0.659e-3 };

            /* 3x2 matrix of p-vectors (cf. Seidelmann 3.591-4, matrix M) */
            double[][][] em = {
            new double[][] {
                new double[] { +0.9999256782, -0.0111820611, -0.0048579477 },
                new double[] { +0.0111820610, +0.9999374784, -0.0000271765 },
                new double[] { +0.0048579479, -0.0000271474, +0.9999881997 }
            },
            new double[][] {
                new double[] { -0.000551,     -0.238565,     +0.435739     },
                new double[] { +0.238514,     -0.002667,     -0.008541     },
                new double[] { -0.435623,     +0.012254,     +0.002117     }
                }
            };

            /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* Spherical coordinates to p-vector. */
            iauS2c(r1950, d1950, out r0);

            /* Adjust p-vector A to give zero proper motion in FK5. */
            w = (bepoch - 1950) / PMF;
            iauPpsp(a, w, ad, ref p);

            /* Remove E-terms. */
            iauPpsp(p, -iauPdp(r0, p), r0, ref p);
            iauPmp(r0, p, ref p);

            /* Convert to Fricke system pv-vector (cf. Seidelmann 3.591-3). */
            for (i = 0; i < 2; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    w = 0.0;
                    for (k = 0; k < 3; k++)
                    {
                        w += em[i][j][k] * p[k];
                    }
                    pv[i][j] = w;
                }
            }

            /* Allow for fictitious proper motion. */
            iauEpb2jd(bepoch, out djm0, out djm);
            w = (iauEpj(djm0, djm) - 2000.0) / PMF;
            iauPvu(w, pv, ref pv);

            /* Revert to spherical coordinates. */
            iauC2s(pv[0], out w, out d2000);
            r2000 = iauAnp(w);
        }


        /// <summary>
        /// Convert J2000.0 FK5 star catalog data to B1950.0 FK4.
        /// </summary>
        /// <param name="r2000">J2000.0 RA (rad)</param>
        /// <param name="d2000">J2000.0 Dec (rad)</param>
        /// <param name="dr2000">J2000.0 RA proper motions (rad/Jul.yr)</param>
        /// <param name="dd2000">J2000.0 Dec proper motions (rad/Jul.yr)</param>
        /// <param name="p2000">parallax (arcsec)</param>
        /// <param name="v2000">radial velocity (km/s, +ve = moving away)</param>
        /// <param name="r1950">B1950.0 RA (rad)</param>
        /// <param name="d1950">B1950.0 Dec (rad)</param>
        /// <param name="dr1950">B1950.0 RA proper motions (rad/trop.yr)</param>
        /// <param name="dd1950">B1950.0 Dec proper motions (rad/trop.yr)</param>
        /// <param name="p1950">parallax (arcsec)</param>
        /// <param name="v1950">radial velocity (km/s, +ve = moving away)</param>
        public static void iauFk524(double r2000, double d2000,
              double dr2000, double dd2000,
              double p2000, double v2000,
              out double r1950, out double d1950,
              out double dr1950, out double dd1950,
              out double p1950, out double v1950)
        {
            /* Radians per year to arcsec per century */
            double PMF = 100.0 * DR2AS;

            /* Small number to avoid arithmetic problems */
            double TINY = 1e-30;

            /* Miscellaneous */
            double r, d, ur, ud, px, rv, pxvf, w, rd;
            int i, j, k, l;

            /* Vectors, p and pv */
            double[][] r0;
            double[][] r1 = { new double[3], new double[3] };
            double[] p1 = new double[3];
            double[] p2 = new double[3];
            double[][] pv = { new double[3], new double[3] };

            /*
            ** CANONICAL CONSTANTS (Seidelmann 1992)
            */

            /* Km per sec to AU per tropical century */
            /* = 86400 * 36524.2198782 / 149597870.7 */
            double VF = 21.095;

            /* Constant pv-vector (cf. Seidelmann 3.591-2, vectors A and Adot) */
            double[][] a = {
                new double[] { -1.62557e-6, -0.31919e-6, -0.13843e-6 },
                new double[] { +1.245e-3,   -1.580e-3,   -0.659e-3   }
            };

            /* 3x2 matrix of pv-vectors (cf. Seidelmann 3.592-1, matrix M^-1) */
            double[][][][] em = {

                new double[][][]{
                    new double[][] {
                        new double[] { +0.9999256795,     +0.0111814828,     +0.0048590039,    },
                        new double[]{ -0.00000242389840, -0.00000002710544, -0.00000001177742 }
                    },

                    new double[][] {
                        new double[] { -0.0111814828,     +0.9999374849,     -0.0000271771,    },
                        new double[]{ +0.00000002710544, -0.00000242392702, +0.00000000006585 }
                    },

                    new double[][]{
                        new double[]  { -0.0048590040,     -0.0000271557,     +0.9999881946,    },
                        new double[]{ +0.00000001177742, +0.00000000006585, -0.00000242404995 }
                    }
                    },

                new double[][][]{
                    new double[][] {
                        new double[] { -0.000551,         +0.238509,         -0.435614,        },
                        new double[]{ +0.99990432,       +0.01118145,       +0.00485852       }
                    },

                    new double[][]{
                        new double[]  { -0.238560,         -0.002667,         +0.012254,        },
                        new double[]{ -0.01118145,       +0.99991613,       -0.00002717       }
                    },

                    new double[][]{
                        new double[]    { +0.435730,         -0.008541,         +0.002117,        },
                        new double[]{ -0.00485852,       -0.00002716,       +0.99996684       }
                    }
                }

            };

            /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* The FK5 data (units radians and arcsec per Julian century). */
            r = r2000;
            d = d2000;
            ur = dr2000 * PMF;
            ud = dd2000 * PMF;
            px = p2000;
            rv = v2000;

            /* Express as a pv-vector. */
            pxvf = px * VF;
            w = rv * pxvf;
            iauS2pv(r, d, 1.0, ur, ud, w, out r0);

            /* Convert pv-vector to Bessel-Newcomb system (cf. Seidelmann 3.592-1). */
            for (i = 0; i < 2; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    w = 0.0;
                    for (k = 0; k < 2; k++)
                    {
                        for (l = 0; l < 3; l++)
                        {
                            w += em[i][j][k][l] * r0[k][l];
                        }
                    }
                    r1[i][j] = w;
                }
            }

            /* Apply E-terms (equivalent to Seidelmann 3.592-3, one iteration). */

            /* Direction. */
            w = iauPm(r1[0]);
            iauSxp(iauPdp(r1[0], a[0]), r1[0], ref p1);
            iauSxp(w, a[0], ref p2);
            iauPmp(p2, p1, ref p1);
            iauPpp(r1[0], p1, ref p1);

            /* Recompute length. */
            w = iauPm(p1);

            /* Direction. */
            iauSxp(iauPdp(r1[0], a[0]), r1[0], ref p1);
            iauSxp(w, a[0], ref p2);
            iauPmp(p2, p1, ref p1);
            iauPpp(r1[0], p1, ref pv[0]);

            /* Derivative. */
            iauSxp(iauPdp(r1[0], a[1]), pv[0], ref p1);
            iauSxp(w, a[1], ref p2);
            iauPmp(p2, p1, ref p1);
            iauPpp(r1[1], p1, ref pv[1]);

            /* Revert to catalog form. */
            iauPv2s(pv, out r, out d, out w, out ur, out ud, out rd);
            if (px > TINY)
            {
                rv = rd / pxvf;
                px = px / w;
            }

            /* Return the results. */
            r1950 = iauAnp(r);
            d1950 = d;
            dr1950 = ur / PMF;
            dd1950 = ud / PMF;
            p1950 = px;
            v1950 = rv;
        }


        /// <summary>
        /// Transform FK5 (J2000.0) star data into the Hipparcos system.
        /// </summary>
        /// <param name="r5">RA (radians)</param>
        /// <param name="d5">Dec (radians)</param>
        /// <param name="dr5">proper motion in RA (dRA/dt, rad/Jyear)</param>
        /// <param name="dd5">proper motion in Dec (dDec/dt, rad/Jyear)</param>
        /// <param name="px5">parallax (arcsec)</param>
        /// <param name="rv5">radial velocity (km/s, positive = receding)</param>
        /// <param name="rh">RA (radians)</param>
        /// <param name="dh">Dec (radians)</param>
        /// <param name="drh">proper motion in RA (dRA/dt, rad/Jyear)</param>
        /// <param name="ddh">proper motion in Dec (dDec/dt, rad/Jyear)</param>
        /// <param name="pxh">parallax (arcsec)</param>
        /// <param name="rvh">radial velocity (km/s, positive = receding)</param>
        public static void iauFk52h(double r5, double d5,
              double dr5, double dd5, double px5, double rv5,
              out double rh, out double dh,
              out double drh, out double ddh, out double pxh, out double rvh)
        {
            int i;
            double[][] pv5;
            double[][] r5h;
            double[] s5h;
            double[] wxp = new double[3];
            double[] vv = new double[3];
            double[][] pvh = { new double[3], new double[3] };


            /* FK5 barycentric position/velocity pv-vector (normalized). */
            iauStarpv(r5, d5, dr5, dd5, px5, rv5, out pv5);

            /* FK5 to Hipparcos orientation matrix and spin vector. */
            iauFk5hip(out r5h, out s5h);

            /* Make spin units per day instead of per year. */
            for (i = 0; i < 3; s5h[i++] /= 365.25) ;

            /* Orient the FK5 position into the Hipparcos system. */
            iauRxp(r5h, pv5[0], ref pvh[0]);

            /* Apply spin to the position giving an extra space motion component. */
            iauPxp(pv5[0], s5h, ref wxp);

            /* Add this component to the FK5 space motion. */
            iauPpp(wxp, pv5[1], ref vv);

            /* Orient the FK5 space motion into the Hipparcos system. */
            iauRxp(r5h, vv, ref pvh[1]);

            /* Hipparcos pv-vector to spherical. */
            iauPvstar(pvh, out rh, out dh, out drh, out ddh, out pxh, out rvh);
        }


        /// <summary>
        /// Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
        /// proper motion in FK5 and parallax.
        /// </summary>
        /// <param name="r2000">J2000.0 FK5 RA (rad)</param>
        /// <param name="d2000">J2000.0 FK5 Dec (rad)</param>
        /// <param name="bepoch">Besselian epoch (e.g. 1950.0)</param>
        /// <param name="r1950">B1950.0 FK4 RA (rad) at epoch BEPOCH</param>
        /// <param name="d1950">B1950.0 FK4 Dec (rad) at epoch BEPOCH</param>
        /// <param name="dr1950">B1950.0 FK4 RA proper motions (rad/trop.yr)</param>
        /// <param name="dd1950">B1950.0 FK4 Dec proper motions (rad/trop.yr)</param>
        public static void iauFk54z(double r2000, double d2000, double bepoch,
              out double r1950, out double d1950,
              out double dr1950, out double dd1950)
        {
            double r, d, pr, pd, px, rv, w;
            double[] p;
            double[] v = new double[3];
            int i;


            /* FK5 equinox J2000.0 to FK4 equinox B1950.0. */
            iauFk524(r2000, d2000, 0.0, 0.0, 0.0, 0.0,
                     out r, out d, out pr, out pd, out px, out rv);

            /* Spherical to Cartesian. */
            iauS2c(r, d, out p);

            /* Fictitious proper motion (radians per year). */
            v[0] = -pr * p[1] - pd * Math.Cos(r) * Math.Sin(d);
            v[1] = pr * p[0] - pd * Math.Sin(r) * Math.Sin(d);
            v[2] = pd * Math.Cos(d);

            /* Apply the motion. */
            w = bepoch - 1950.0;
            for (i = 0; i < 3; i++)
            {
                p[i] += w * v[i];
            }

            /* Cartesian to spherical. */
            iauC2s(p, out w, out d1950);
            r1950 = iauAnp(w);

            /* Fictitious proper motion. */
            dr1950 = pr;
            dd1950 = pd;
        }


        /// <summary>
        /// FK5 to Hipparcos rotation and spin.
        /// </summary>
        /// <param name="r5h">double[3][3]  r-matrix: FK5 rotation wrt Hipparcos</param>
        /// <param name="s5h">double[3]     r-vector: FK5 spin wrt Hipparcos</param>
        public static void iauFk5hip(out double[][] r5h, out double[] s5h)
        {
            double[] v = new double[3];
            s5h = new double[3];

            /* FK5 wrt Hipparcos orientation and spin (radians, radians/year) */
            double epx, epy, epz;
            double omx, omy, omz;


            epx = -19.9e-3 * DAS2R;
            epy = -9.1e-3 * DAS2R;
            epz = 22.9e-3 * DAS2R;

            omx = -0.30e-3 * DAS2R;
            omy = 0.60e-3 * DAS2R;
            omz = 0.70e-3 * DAS2R;

            /* FK5 to Hipparcos orientation expressed as an r-vector. */
            v[0] = epx;
            v[1] = epy;
            v[2] = epz;

            /* Re-express as an r-matrix. */
            iauRv2m(v, out r5h);

            /* Hipparcos wrt FK5 spin expressed as an r-vector. */
            s5h[0] = omx;
            s5h[1] = omy;
            s5h[2] = omz;
        }


        /// <summary>
        /// Transform an FK5 (J2000.0) star position into the system of the
        /// Hipparcos catalogue, assuming zero Hipparcos proper motion.
        /// </summary>
        /// <param name="r5">FK5 RA (radians), equinox J2000.0, at date</param>
        /// <param name="d5">FK5 Dec (radians), equinox J2000.0, at date</param>
        /// <param name="date1">TDB date</param>
        /// <param name="date2">TDB date</param>
        /// <param name="rh">Hipparcos RA (radians)</param>
        /// <param name="dh">Hipparcos Dec (radians)</param>
        public static void iauFk5hz(double r5, double d5, double date1, double date2,
              out double rh, out double dh)
        {
            double t, w;
            double[] p5e;
            double[][] r5h;
            double[] s5h;
            double[] vst = new double[3];
            double[][] rst;
            double[] p5 = new double[3];
            double[] ph = new double[3];

            /* Interval from given date to fundamental epoch J2000.0 (JY). */
            t = -((date1 - DJ00) + date2) / DJY;

            /* FK5 barycentric position vector. */
            iauS2c(r5, d5, out p5e);

            /* FK5 to Hipparcos orientation matrix and spin vector. */
            iauFk5hip(out r5h, out s5h);

            /* Accumulated Hipparcos wrt FK5 spin over that interval. */
            iauSxp(t, s5h, ref vst);

            /* Express the accumulated spin as a rotation matrix. */
            iauRv2m(vst, out rst);

            /* Derotate the vector's FK5 axes back to date. */
            iauTrxp(rst, p5e, ref p5);

            /* Rotate the vector into the Hipparcos system. */
            iauRxp(r5h, p5, ref ph);

            /* Hipparcos vector to spherical. */
            iauC2s(ph, out w, out dh);
            rh = iauAnp(w);
        }


        /// <summary>
        /// Transform Hipparcos star data into the FK5 (J2000.0) system.
        /// </summary>
        /// <param name="rh">RA (radians)</param>
        /// <param name="dh">Dec (radians)</param>
        /// <param name="drh">proper motion in RA (dRA/dt, rad/Jyear)</param>
        /// <param name="ddh">proper motion in Dec (dDec/dt, rad/Jyear)</param>
        /// <param name="pxh">parallax (arcsec)</param>
        /// <param name="rvh">radial velocity (km/s, positive = receding)</param>
        /// <param name="r5">RA (radians)</param>
        /// <param name="d5">Dec (radians)</param>
        /// <param name="dr5">proper motion in RA (dRA/dt, rad/Jyear)</param>
        /// <param name="dd5">proper motion in Dec (dDec/dt, rad/Jyear)</param>
        /// <param name="px5">parallax (arcsec)</param>
        /// <param name="rv5">radial velocity (km/s, positive = receding)</param>
        public static void iauH2fk5(double rh, double dh,
              double drh, double ddh, double pxh, double rvh,
              out double r5, out double d5,
              out double dr5, out double dd5, out double px5, out double rv5)
        {
            int i;
            double[][] pvh;
            double[][] r5h;
            double[] s5h;
            double[] sh = new double[3];
            double[] wxp = new double[3];
            double[] vv = new double[3];
            double[][] pv5 = { new double[3], new double[3] };


            /* Hipparcos barycentric position/velocity pv-vector (normalized). */
            iauStarpv(rh, dh, drh, ddh, pxh, rvh, out pvh);

            /* FK5 to Hipparcos orientation matrix and spin vector. */
            iauFk5hip(out r5h, out s5h);

            /* Make spin units per day instead of per year. */
            for (i = 0; i < 3; s5h[i++] /= 365.25) ;

            /* Orient the spin into the Hipparcos system. */
            iauRxp(r5h, s5h, ref sh);

            /* De-orient the Hipparcos position into the FK5 system. */
            iauTrxp(r5h, pvh[0], ref pv5[0]);

            /* Apply spin to the position giving an extra space motion component. */
            iauPxp(pvh[0], sh, ref wxp);

            /* Subtract this component from the Hipparcos space motion. */
            iauPmp(pvh[1], wxp, ref vv);

            /* De-orient the Hipparcos space motion into the FK5 system. */
            iauTrxp(r5h, vv, ref pv5[1]);

            /* FK5 pv-vector to spherical. */
            iauPvstar(pv5, out r5, out d5, out dr5, out dd5, out px5, out rv5);
        }


        /// <summary>
        /// Transform a Hipparcos star position into FK5 J2000.0, assuming
        /// zero Hipparcos proper motion.
        /// </summary>
        /// <param name="rh">Hipparcos RA (radians)</param>
        /// <param name="dh">Hipparcos Dec (radians)</param>
        /// <param name="date1">TDB date</param>
        /// <param name="date2">TDB date</param>
        /// <param name="r5">RA (radians)</param>
        /// <param name="d5">Dec (radians)</param>
        /// <param name="dr5">FK5 RA proper motion (rad/year)</param>
        /// <param name="dd5">Dec proper motion (rad/year)</param>
        public static void iauHfk5z(double rh, double dh, double date1, double date2,
              out double r5, out double d5, out double dr5, out double dd5)
        {
            double t, w, r, v;
            double[] ph;
            double[][] r5h;
            double[] s5h;
            double[] sh = new double[3];
            double[] vst = new double[3];
            double[][] rst;
            double[][] r5ht = { new double[3], new double[3], new double[3] };
            double[][] pv5e = { new double[3], new double[3] };
            double[] vv = new double[3];


            /* Time interval from fundamental epoch J2000.0 to given date (JY). */
            t = ((date1 - DJ00) + date2) / DJY;

            /* Hipparcos barycentric position vector (normalized). */
            iauS2c(rh, dh, out ph);

            /* FK5 to Hipparcos orientation matrix and spin vector. */
            iauFk5hip(out r5h, out s5h);

            /* Rotate the spin into the Hipparcos system. */
            iauRxp(r5h, s5h, ref sh);

            /* Accumulated Hipparcos wrt FK5 spin over that interval. */
            iauSxp(t, s5h, ref vst);

            /* Express the accumulated spin as a rotation matrix. */
            iauRv2m(vst, out rst);

            /* Rotation matrix:  accumulated spin, then FK5 to Hipparcos. */
            iauRxr(r5h, rst, ref r5ht);

            /* De-orient & de-spin the Hipparcos position into FK5 J2000.0. */
            iauTrxp(r5ht, ph, ref pv5e[0]);

            /* Apply spin to the position giving a space motion. */
            iauPxp(sh, ph, ref vv);

            /* De-orient & de-spin the Hipparcos space motion into FK5 J2000.0. */
            iauTrxp(r5ht, vv, ref pv5e[1]);

            /* FK5 position/velocity pv-vector to spherical. */
            iauPv2s(pv5e, out w, out d5, out r, out dr5, out dd5, out v);
            r5 = iauAnp(w);
        }


        /// <summary>
        /// Star proper motion:  update star catalog data for space motion.
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
        public static int iauStarpm(double ra1, double dec1,
              double pmr1, double pmd1, double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              out double ra2, out double dec2,
              out double pmr2, out double pmd2, out double px2, out double rv2)
        {
            double[][] pv1;
            double tl1, dt;
            double[][] pv = { new double[3], new double[3] };
            double r2, rdv, v2, c2mv2, tl2;
            double[][] pv2 = { new double[3], new double[3] };
            int j1, j2, j;

            ra2 = 0; dec2 = 0; pmr2 = 0; pmd2 = 0; px2 = 0; rv2 = 0;


            /* RA,Dec etc. at the "before" epoch to space motion pv-vector. */
            j1 = iauStarpv(ra1, dec1, pmr1, pmd1, px1, rv1, out pv1);

            /* Light time when observed (days). */
            tl1 = iauPm(pv1[0]) / DC;

            /* Time interval, "before" to "after" (days). */
            dt = (ep2a - ep1a) + (ep2b - ep1b);

            /* Move star along track from the "before" observed position to the */
            /* "after" geometric position. */
            iauPvu(dt + tl1, pv1, ref pv);

            /* From this geometric position, deduce the observed light time (days) */
            /* at the "after" epoch (with theoretically unneccessary error check). */
            r2 = iauPdp(pv[0], pv[0]);
            rdv = iauPdp(pv[0], pv[1]);
            v2 = iauPdp(pv[1], pv[1]);
            c2mv2 = DC * DC - v2;
            if (c2mv2 <= 0) return -1;
            tl2 = (-rdv + Math.Sqrt(rdv * rdv + c2mv2 * r2)) / c2mv2;

            /* Move the position along track from the observed place at the */
            /* "before" epoch to the observed place at the "after" epoch. */
            iauPvu(dt + (tl1 - tl2), pv1, ref pv2);

            /* Space motion pv-vector to RA,Dec etc. at the "after" epoch. */
            j2 = iauPvstar(pv2, out ra2, out dec2, out pmr2, out pmd2, out px2, out rv2);

            /* Final status. */
            j = (j2 == 0) ? j1 : -1;

            return j;
        }


    }
}
