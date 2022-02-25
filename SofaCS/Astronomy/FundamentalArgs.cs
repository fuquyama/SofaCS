using System;

namespace SofaCS
{
    public static partial class Sofa
    {
        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean elongation of the Moon from the Sun.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>D, radians</returns>
        public static double iauFad03(double t)
        {
            double a;


            /* Mean elongation of the Moon from the Sun (IERS Conventions 2003). */
            a = ((1072260.703692 +
                t * (1602961601.2090 +
                t * (-6.3706 +
                t * (0.006593 +
                t * (-0.00003169))))) % TURNAS) * DAS2R;

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Earth.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Earth, radians</returns>
        public static double iauFae03(double t)
        {
            double a;


            /* Mean longitude of Earth (IERS Conventions 2003). */
            a = (1.753470314 + 628.3075849991 * t) % D2PI;

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of the Moon minus mean longitude of the ascending node.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>F, radians</returns>
        public static double iauFaf03(double t)
        {
            double a;


            /* Mean longitude of the Moon minus that of the ascending node */
            /* (IERS Conventions 2003).                                    */
            a = ((335779.526232 +
                      t * (1739527262.8478 +
                      t * (-12.7512 +
                      t * (-0.001037 +
                      t * (0.00000417))))) % TURNAS) * DAS2R;

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Jupiter.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Jupiter, radians</returns>
        public static double iauFaju03(double t)
        {
            double a;


            /* Mean longitude of Jupiter (IERS Conventions 2003). */
            a = ((0.599546497 + 52.9690962641 * t) % D2PI);

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean anomaly of the Moon.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>l, radians</returns>
        public static double iauFal03(double t)
        {
            double a;


            /* Mean anomaly of the Moon (IERS Conventions 2003). */
            a = ((485868.249036 +
                      t * (1717915923.2178 +
                      t * (31.8792 +
                      t * (0.051635 +
                      t * (-0.00024470))))) % TURNAS) * DAS2R;

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean anomaly of the Sun.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>l', radians</returns>
        public static double iauFalp03(double t)
        {
            double a;


            /* Mean anomaly of the Sun (IERS Conventions 2003). */
            a = ((1287104.793048 +
                      t * (129596581.0481 +
                      t * (-0.5532 +
                      t * (0.000136 +
                      t * (-0.00001149))))) % TURNAS) * DAS2R;

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Mars.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Mars, radians</returns>
        public static double iauFama03(double t)
        {
            double a;


            /* Mean longitude of Mars (IERS Conventions 2003). */
            a = ((6.203480913 + 334.0612426700 * t) % D2PI);

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Mercury.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Mercury, radians</returns>
        public static double iauFame03(double t)
        {
            double a;


            /* Mean longitude of Mercury (IERS Conventions 2003). */
            a = ((4.402608842 + 2608.7903141574 * t) % D2PI);

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Neptune.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Neptune, radians</returns>
        public static double iauFane03(double t)
        {
            double a;


            /* Mean longitude of Neptune (IERS Conventions 2003). */
            a = ((5.311886287 + 3.8133035638 * t) % D2PI);

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of the Moon's ascending node.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>Omega, radians</returns>
        public static double iauFaom03(double t)
        {
            double a;


            /* Mean longitude of the Moon's ascending node */
            /* (IERS Conventions 2003).                    */
            a = ((450160.398036 +
                      t * (-6962890.5431 +
                      t * (7.4722 +
                      t * (0.007702 +
                      t * (-0.00005939))))) % TURNAS) * DAS2R;

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// general accumulated precession in longitude.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>general precession in longitude, radians</returns>
        public static double iauFapa03(double t)
        {
            double a;


            /* General accumulated precession in longitude. */
            a = (0.024381750 + 0.00000538691 * t) * t;

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Saturn.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Saturn, radians</returns>
        public static double iauFasa03(double t)
        {
            double a;


            /* Mean longitude of Saturn (IERS Conventions 2003). */
            a = ((0.874016757 + 21.3299104960 * t) % D2PI);

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Uranus.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Uranus, radians</returns>
        public static double iauFaur03(double t)
        {
            double a;


            /* Mean longitude of Uranus (IERS Conventions 2003). */
            a = ((5.481293872 + 7.4781598567 * t) % D2PI);

            return a;
        }


        /// <summary>
        /// Fundamental argument, IERS Conventions (2003):
        /// mean longitude of Venus.
        /// </summary>
        /// <param name="t">TDB, Julian centuries since J2000.0</param>
        /// <returns>mean longitude of Venus, radians</returns>
        public static double iauFave03(double t)
        {
            double a;


            /* Mean longitude of Venus (IERS Conventions 2003). */
            a = ((3.176146697 + 1021.3285546211 * t) % D2PI);

            return a;
        }


    }
}
