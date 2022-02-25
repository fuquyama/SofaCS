using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// P-vector to spherical coordinates.
        /// </summary>
        /// <param name="p">double[3]    p-vector</param>
        /// <param name="theta">longitude angle (radians)</param>
        /// <param name="phi">latitude angle (radians)</param>
        public static void iauC2s(in double[] p, out double theta, out double phi)
        {
            double x, y, z, d2;

            x = p[0];
            y = p[1];
            z = p[2];
            d2 = x * x + y * y;

            theta = (d2 == 0.0) ? 0.0 : Math.Atan2(y, x);
            phi = (z == 0.0) ? 0.0 : Math.Atan2(z, Math.Sqrt(d2));
        }


        /// <summary>
        /// P-vector to spherical polar coordinates.
        /// </summary>
        /// <param name="p">double[3]    p-vector</param>
        /// <param name="theta">longitude angle (radians)</param>
        /// <param name="phi">latitude angle (radians)</param>
        /// <param name="r">radial distance</param>
        public static void iauP2s(in double[] p, out double theta, out double phi, out double r)
        {
            iauC2s(p, out theta, out phi);
            r = iauPm(p);
        }


        /// <summary>
        /// Convert position/velocity from Cartesian to spherical coordinates.
        /// </summary>
        /// <param name="pv">double[2][3]  pv-vector</param>
        /// <param name="theta">longitude angle (radians)</param>
        /// <param name="phi">latitude angle (radians)</param>
        /// <param name="r">radial distance</param>
        /// <param name="td">rate of change of theta</param>
        /// <param name="pd">rate of change of phi</param>
        /// <param name="rd">rate of change of r</param>
        public static void iauPv2s(in double[][] pv,
             out double theta, out double phi, out double r,
             out double td, out double pd, out double rd)
        {
            double x, y, z, xd, yd, zd, rxy2, rxy, r2, rtrue, rw, xyp;


            /* Components of position/velocity vector. */
            x = pv[0][0];
            y = pv[0][1];
            z = pv[0][2];
            xd = pv[1][0];
            yd = pv[1][1];
            zd = pv[1][2];

            /* Component of r in XY plane squared. */
            rxy2 = x * x + y * y;

            /* Modulus squared. */
            r2 = rxy2 + z * z;

            /* Modulus. */
            rtrue = Math.Sqrt(r2);

            /* If null vector, move the origin along the direction of movement. */
            rw = rtrue;
            if (rtrue == 0.0)
            {
                x = xd;
                y = yd;
                z = zd;
                rxy2 = x * x + y * y;
                r2 = rxy2 + z * z;
                rw = Math.Sqrt(r2);
            }

            /* Position and velocity in spherical coordinates. */
            rxy = Math.Sqrt(rxy2);
            xyp = x * xd + y * yd;
            if (rxy2 != 0.0)
            {
                theta = Math.Atan2(y, x);
                phi = Math.Atan2(z, rxy);
                td = (x * yd - y * xd) / rxy2;
                pd = (zd * rxy2 - z * xyp) / (r2 * rxy);
            }
            else
            {
                theta = 0.0;
                phi = (z != 0.0) ? Math.Atan2(z, rxy) : 0.0;
                td = 0.0;
                pd = 0.0;
            }
            r = rtrue;
            rd = (rw != 0.0) ? (xyp + z * zd) / rw : 0.0;
        }


        /// <summary>
        /// Convert spherical coordinates to Cartesian.
        /// </summary>
        /// <param name="theta">longitude angle (radians)</param>
        /// <param name="phi">latitude angle (radians)</param>
        /// <param name="c">double[3]    direction cosines</param>
        public static void iauS2c(double theta, double phi, out double[] c)
        {
            c = new double[3];

            double cp;

            cp = Math.Cos(phi);
            c[0] = Math.Cos(theta) * cp;
            c[1] = Math.Sin(theta) * cp;
            c[2] = Math.Sin(phi);
        }


        /// <summary>
        /// Convert spherical polar coordinates to p-vector.
        /// </summary>
        /// <param name="theta">longitude angle (radians)</param>
        /// <param name="phi">latitude angle (radians)</param>
        /// <param name="r">radial distance</param>
        /// <param name="p">double[3]    Cartesian coordinates</param>
        public static void iauS2p(double theta, double phi, double r, out double[] p)
        {
            p = new double[3];
            double[] u;


            iauS2c(theta, phi, out u);
            iauSxp(r, u, ref p);
        }


        /// <summary>
        /// Convert position/velocity from spherical to Cartesian coordinates.
        /// </summary>
        /// <param name="theta">longitude angle (radians)</param>
        /// <param name="phi">latitude angle (radians)</param>
        /// <param name="r">radial distance</param>
        /// <param name="td">rate of change of theta</param>
        /// <param name="pd">rate of change of phi</param>
        /// <param name="rd">rate of change of r</param>
        /// <param name="pv">double[2][3]    pv-vector</param>
        public static void iauS2pv(double theta, double phi, double r,
             double td, double pd, double rd,
             out double[][] pv)
        {
            pv = new double[][] { new double[3], new double[3] };

            double st, ct, sp, cp, rcp, x, y, rpd, w;


            st = Math.Sin(theta);
            ct = Math.Cos(theta);
            sp = Math.Sin(phi);
            cp = Math.Cos(phi);
            rcp = r * cp;
            x = rcp * ct;
            y = rcp * st;
            rpd = r * pd;
            w = rpd * sp - cp * rd;

            pv[0][0] = x;
            pv[0][1] = y;
            pv[0][2] = r * sp;
            pv[1][0] = -y * td - w * ct;
            pv[1][1] = x * td - w * st;
            pv[1][2] = rpd * cp + sp * rd;
        }


    }
}
