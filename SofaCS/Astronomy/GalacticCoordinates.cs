using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Transformation from Galactic Coordinates to ICRS.
        /// </summary>
        /// <param name="dl">galactic longitude (radians)</param>
        /// <param name="db">galactic latitude (radians)</param>
        /// <param name="dr">ICRS right ascension (radians)</param>
        /// <param name="dd">ICRS declination (radians)</param>
        public static void iauG2icrs(double dl, double db, out double dr, out double dd)
        {
            double[] v1;
            double[] v2 = new double[3];

            /*
            **  L2,B2 system of galactic coordinates in the form presented in the
            **  Hipparcos Catalogue.  In degrees:
            **
            **  P = 192.85948    right ascension of the Galactic north pole in ICRS
            **  Q =  27.12825    declination of the Galactic north pole in ICRS
            **  R =  32.93192    Galactic longitude of the ascending node of
            **                   the Galactic equator on the ICRS equator
            **
            **  ICRS to galactic rotation matrix, obtained by computing
            **  R_3(-R) R_1(pi/2-Q) R_3(pi/2+P) to the full precision shown:
            */
            double[][] r = {
                new double[] {
                    -0.054875560416215368492398900454,
                    -0.873437090234885048760383168409,
                    -0.483835015548713226831774175116
                },
                new double[] {
                    +0.494109427875583673525222371358,
                    -0.444829629960011178146614061616,
                    +0.746982244497218890527388004556
                },
                new double[] {
                    -0.867666149019004701181616534570,
                    -0.198076373431201528180486091412,
                    +0.455983776175066922272100478348
                }
            };


            /* Spherical to Cartesian. */
            iauS2c(dl, db, out v1);

            /* Galactic to ICRS. */
            iauTrxp(r, v1, ref v2);

            /* Cartesian to spherical. */
            iauC2s(v2, out dr, out dd);

            /* Express in conventional ranges. */
            dr = iauAnp(dr);
            dd = iauAnpm(dd);
        }


        /// <summary>
        /// Transformation from ICRS to Galactic Coordinates.
        /// </summary>
        /// <param name="dr">ICRS right ascension (radians)</param>
        /// <param name="dd">ICRS declination (radians)</param>
        /// <param name="dl">galactic longitude (radians)</param>
        /// <param name="db">galactic latitude (radians)</param>
        public static void iauIcrs2g(double dr, double dd, out double dl, out double db)
        {
            double[] v1;
            double[] v2 = new double[3];

            /*
            **  L2,B2 system of galactic coordinates in the form presented in the
            **  Hipparcos Catalogue.  In degrees:
            **
            **  P = 192.85948    right ascension of the Galactic north pole in ICRS
            **  Q =  27.12825    declination of the Galactic north pole in ICRS
            **  R =  32.93192    Galactic longitude of the ascending node of
            **                   the Galactic equator on the ICRS equator
            **
            **  ICRS to galactic rotation matrix, obtained by computing
            **  R_3(-R) R_1(pi/2-Q) R_3(pi/2+P) to the full precision shown:
            */
            double[][] r = {
                new double[]{
                    -0.054875560416215368492398900454,
                    -0.873437090234885048760383168409,
                    -0.483835015548713226831774175116
                },
                new double[]{
                    +0.494109427875583673525222371358,
                    -0.444829629960011178146614061616,
                    +0.746982244497218890527388004556
                },
                new double[]{
                    -0.867666149019004701181616534570,
                    -0.198076373431201528180486091412,
                    +0.455983776175066922272100478348
                }
            };


            /* Spherical to Cartesian. */
            iauS2c(dr, dd, out v1);

            /* ICRS to Galactic. */
            iauRxp(r, v1, ref v2);

            /* Cartesian to spherical. */
            iauC2s(v2, out dl, out db);

            /* Express in conventional ranges. */
            dl = iauAnp(dl);
            db = iauAnpm(db);
        }


    }
}
