using System;
using SofaCS;

namespace Validate_the_SOFA_C_functions
{
    internal class Program
    {
        static int verbose = 0;

        static void viv(int ival, int ivalok,
                string func, string test, ref int status)
        /*
        **  - - - -
        **   v i v
        **  - - - -
        **
        **  Validate an integer result.
        **
        **  Internal function used by t_sofa_c program.
        **
        **  Given:
        **     ival     int          value computed by function under test
        **     ivalok   int          correct value
        **     func     char[]       name of function under test
        **     test     char[]       name of individual test
        **
        **  Given and returned:
        **     status   int          set to TRUE if test fails
        **
        **  This revision:  2013 August 7
        */
        {

            if (ival != ivalok)
            {
                status = 1;
                Console.Write("{0} failed: {1} want {2} got {3}\n",
                       func, test, ivalok, ival);
            }
            else if (verbose > 0)
            {
                Console.Write("{0} passed: {1} want {2} got {3}\n",
                              func, test, ivalok, ival);
            }

        }


        static void vvd(double val, double valok, double dval,
                string func, string test, ref int status)
        /*
        **  - - - -
        **   v v d
        **  - - - -
        **
        **  Validate a double result.
        **
        **  Internal function used by t_sofa_c program.
        **
        **  Given:
        **     val      double       value computed by function under test
        **     valok    double       expected value
        **     dval     double       maximum allowable error
        **     func     char[]       name of function under test
        **     test     char[]       name of individual test
        **
        **  Given and returned:
        **     status   int          set to TRUE if test fails
        **
        **  This revision:  2016 April 21
        */
        {
            double a, f;   /* absolute and fractional error */


            a = val - valok;
            if (a != 0.0 &&  Math.Abs(a) > Math.Abs(dval))
            {
                f = Math.Abs(valok / a);
                status = 1;
                Console.Write("{0} failed: {1} want {2:G20} got {3:G20} (1/{4:G3})\n",
             func, test, valok, val, f);
            }
            else if (verbose > 0)
            {
                Console.Write("{0} passed: {1} want {2:G20} got {3:G20}\n",
             func, test, valok, val);
            }

        }



        static void t_a2af(ref int status)
        /*
        **  - - - - - - -
        **   t _ a 2 a f
        **  - - - - - - -
        **
        **  Test iauA2af function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauA2af, viv
        **
        **  This revision:  2013 August 7
        */
        {
            int[] idmsf = new int[4];
            char s;


            Sofa.iauA2af(4, 2.345, out s, out idmsf);

            viv(s, '+', "iauA2af", "s", ref status);

            viv(idmsf[0], 134, "iauA2af", "0", ref status);
            viv(idmsf[1], 21, "iauA2af", "1", ref status);
            viv(idmsf[2], 30, "iauA2af", "2", ref status);
            viv(idmsf[3], 9706, "iauA2af", "3", ref status);

        }

        static void t_a2tf(ref int status)
        /*
        **  - - - - - - -
        **   t _ a 2 t f
        **  - - - - - - -
        **
        **  Test iauA2tf function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauA2tf, viv
        **
        **  This revision:  2013 August 7
        */
        {
            int[] ihmsf = new int[4];
            char s;


            Sofa.iauA2tf(4, -3.01234, out s, out ihmsf);

            viv((int)s, '-', "iauA2tf", "s", ref status);

            viv(ihmsf[0], 11, "iauA2tf", "0", ref status);
            viv(ihmsf[1], 30, "iauA2tf", "1", ref status);
            viv(ihmsf[2], 22, "iauA2tf", "2", ref status);
            viv(ihmsf[3], 6484, "iauA2tf", "3", ref status);

        }

        static void t_ab(ref int status)
        /*
        **  - - - - -
        **   t _ a b
        **  - - - - -
        **
        **  Test iauAb function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAb, vvd
        **
        **  This revision:  2013 October 1
        */
        {
            double s, bm1;
            double[] pnat = new double[3];
            double[] v = new double[3];
            double[] ppr;


            pnat[0] = -0.76321968546737951;
            pnat[1] = -0.60869453983060384;
            pnat[2] = -0.21676408580639883;
            v[0] = 2.1044018893653786e-5;
            v[1] = -8.9108923304429319e-5;
            v[2] = -3.8633714797716569e-5;
            s = 0.99980921395708788;
            bm1 = 0.99999999506209258;

            Sofa.iauAb(pnat, v, s, bm1, out ppr);

            vvd(ppr[0], -0.7631631094219556269, 1e-12, "iauAb", "1", ref status);
            vvd(ppr[1], -0.6087553082505590832, 1e-12, "iauAb", "2", ref status);
            vvd(ppr[2], -0.2167926269368471279, 1e-12, "iauAb", "3", ref status);

        }

        static void t_ae2hd(ref int status)
        /*
        **  - - - - - - - -
        **   t _ a e 2 h d
        **  - - - - - - - -
        **
        **  Test iauAe2hd function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAe2hd and vvd
        **
        **  This revision:  2017 October 21
        */
        {
            double a, e, p, h, d;


            a = 5.5;
            e = 1.1;
            p = 0.7;

            Sofa.iauAe2hd(a, e, p, out h, out d);

            vvd(h, 0.5933291115507309663, 1e-14, "iauAe2hd", "h", ref status);
            vvd(d, 0.9613934761647817620, 1e-14, "iauAe2hd", "d", ref status);

        }

        static void t_af2a(ref int status)
        /*
        **  - - - - - - -
        **   t _ a f 2 a
        **  - - - - - - -
        **
        **  Test iauAf2a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAf2a, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double a;
            int j;


            j = Sofa.iauAf2a('-', 45, 13, 27.2, out a);

            vvd(a, -0.7893115794313644842, 1e-12, "iauAf2a", "a", ref status);
            viv(j, 0, "iauAf2a", "j", ref status);

        }

        static void t_anp(ref int status)
        /*
        **  - - - - - -
        **   t _ a n p
        **  - - - - - -
        **
        **  Test iauAnp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAnp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauAnp(-0.1), 6.183185307179586477, 1e-12, "iauAnp", "", ref status);
        }

        static void t_anpm(ref int status)
        /*
        **  - - - - - - -
        **   t _ a n p m
        **  - - - - - - -
        **
        **  Test iauAnpm function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAnpm, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauAnpm(-4.0), 2.283185307179586477, 1e-12, "iauAnpm", "", ref status);
        }

        static void t_apcg(ref int status)
        /*
        **  - - - - - - -
        **   t _ a p c g
        **  - - - - - - -
        **
        **  Test iauApcg function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApcg, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2;
            double[][] ebpv = { new double[3], new double[3] };
            double[] ehp = new double[3];
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;
            ebpv[0][0] = 0.901310875;
            ebpv[0][1] = -0.417402664;
            ebpv[0][2] = -0.180982288;
            ebpv[1][0] = 0.00742727954;
            ebpv[1][1] = 0.0140507459;
            ebpv[1][2] = 0.00609045792;
            ehp[0] = 0.903358544;
            ehp[1] = -0.415395237;
            ehp[2] = -0.180084014;

            Sofa.iauApcg(date1, date2, ebpv, ehp, ref astrom);

            vvd(astrom.pmt, 12.65133794027378508, 1e-11,
                            "iauApcg", "pmt", ref status);
            vvd(astrom.eb[0], 0.901310875, 1e-12,
                              "iauApcg", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.417402664, 1e-12,
                              "iauApcg", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.180982288, 1e-12,
                              "iauApcg", "eb(3)", ref status);
            vvd(astrom.eh[0], 0.8940025429324143045, 1e-12,
                              "iauApcg", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.4110930268679817955, 1e-12,
                              "iauApcg", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.1782189004872870264, 1e-12,
                              "iauApcg", "eh(3)", ref status);
            vvd(astrom.em, 1.010465295811013146, 1e-12,
                           "iauApcg", "em", ref status);
            vvd(astrom.v[0], 0.4289638913597693554e-4, 1e-16,
                             "iauApcg", "v(1)", ref status);
            vvd(astrom.v[1], 0.8115034051581320575e-4, 1e-16,
                             "iauApcg", "v(2)", ref status);
            vvd(astrom.v[2], 0.3517555136380563427e-4, 1e-16,
                             "iauApcg", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999951686012981, 1e-12,
                            "iauApcg", "bm1", ref status);
            vvd(astrom.bpn[0][0], 1.0, 0.0,
                                  "iauApcg", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0.0, 0.0,
                                  "iauApcg", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0.0, 0.0,
                                  "iauApcg", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], 0.0, 0.0,
                                  "iauApcg", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 1.0, 0.0,
                                  "iauApcg", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], 0.0, 0.0,
                                  "iauApcg", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], 0.0, 0.0,
                                  "iauApcg", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0.0, 0.0,
                                  "iauApcg", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 1.0, 0.0,
                                  "iauApcg", "bpn(3,3)", ref status);

        }

        static void t_apcg13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a p c g 1 3
        **  - - - - - - - - -
        **
        **  Test iauApcg13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApcg13, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2;
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;

            Sofa.iauApcg13(date1, date2, ref astrom);

            vvd(astrom.pmt, 12.65133794027378508, 1e-11,
                            "iauApcg13", "pmt", ref status);
            vvd(astrom.eb[0], 0.9013108747340644755, 1e-12,
                            "iauApcg13", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.4174026640406119957, 1e-12,
                            "iauApcg13", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.1809822877867817771, 1e-12,
                            "iauApcg13", "eb(3)", ref status);
            vvd(astrom.eh[0], 0.8940025429255499549, 1e-12,
                            "iauApcg13", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.4110930268331896318, 1e-12,
                            "iauApcg13", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.1782189006019749850, 1e-12,
                            "iauApcg13", "eh(3)", ref status);
            vvd(astrom.em, 1.010465295964664178, 1e-12,
                            "iauApcg13", "em", ref status);
            vvd(astrom.v[0], 0.4289638912941341125e-4, 1e-16,
                            "iauApcg13", "v(1)", ref status);
            vvd(astrom.v[1], 0.8115034032405042132e-4, 1e-16,
                            "iauApcg13", "v(2)", ref status);
            vvd(astrom.v[2], 0.3517555135536470279e-4, 1e-16,
                            "iauApcg13", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999951686013142, 1e-12,
                            "iauApcg13", "bm1", ref status);
            vvd(astrom.bpn[0][0], 1.0, 0.0,
                                  "iauApcg13", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0.0, 0.0,
                                  "iauApcg13", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0.0, 0.0,
                                  "iauApcg13", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], 0.0, 0.0,
                                  "iauApcg13", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 1.0, 0.0,
                                  "iauApcg13", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], 0.0, 0.0,
                                  "iauApcg13", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], 0.0, 0.0,
                                  "iauApcg13", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0.0, 0.0,
                                  "iauApcg13", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 1.0, 0.0,
                                  "iauApcg13", "bpn(3,3)", ref status);

        }

        static void t_apci(ref int status)
        /*
        **  - - - - - - -
        **   t _ a p c i
        **  - - - - - - -
        **
        **  Test iauApci function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApci, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2, x, y, s;
            double[][] ebpv = { new double[3], new double[3] };
            double[] ehp = new double[3];
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;
            ebpv[0][0] = 0.901310875;
            ebpv[0][1] = -0.417402664;
            ebpv[0][2] = -0.180982288;
            ebpv[1][0] = 0.00742727954;
            ebpv[1][1] = 0.0140507459;
            ebpv[1][2] = 0.00609045792;
            ehp[0] = 0.903358544;
            ehp[1] = -0.415395237;
            ehp[2] = -0.180084014;
            x = 0.0013122272;
            y = -2.92808623e-5;
            s = 3.05749468e-8;

            Sofa.iauApci(date1, date2, ebpv, ehp, x, y, s, ref astrom);

            vvd(astrom.pmt, 12.65133794027378508, 1e-11,
                            "iauApci", "pmt", ref status);
            vvd(astrom.eb[0], 0.901310875, 1e-12,
                              "iauApci", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.417402664, 1e-12,
                              "iauApci", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.180982288, 1e-12,
                              "iauApci", "eb(3)", ref status);
            vvd(astrom.eh[0], 0.8940025429324143045, 1e-12,
                              "iauApci", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.4110930268679817955, 1e-12,
                              "iauApci", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.1782189004872870264, 1e-12,
                              "iauApci", "eh(3)", ref status);
            vvd(astrom.em, 1.010465295811013146, 1e-12,
                           "iauApci", "em", ref status);
            vvd(astrom.v[0], 0.4289638913597693554e-4, 1e-16,
                             "iauApci", "v(1)", ref status);
            vvd(astrom.v[1], 0.8115034051581320575e-4, 1e-16,
                             "iauApci", "v(2)", ref status);
            vvd(astrom.v[2], 0.3517555136380563427e-4, 1e-16,
                             "iauApci", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999951686012981, 1e-12,
                            "iauApci", "bm1", ref status);
            vvd(astrom.bpn[0][0], 0.9999991390295159156, 1e-12,
                                  "iauApci", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0.4978650072505016932e-7, 1e-12,
                                  "iauApci", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0.1312227200000000000e-2, 1e-12,
                                  "iauApci", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], -0.1136336653771609630e-7, 1e-12,
                                  "iauApci", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 0.9999999995713154868, 1e-12,
                                  "iauApci", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], -0.2928086230000000000e-4, 1e-12,
                                  "iauApci", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], -0.1312227200895260194e-2, 1e-12,
                                  "iauApci", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0.2928082217872315680e-4, 1e-12,
                                  "iauApci", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 0.9999991386008323373, 1e-12,
                                  "iauApci", "bpn(3,3)", ref status);

        }

        static void t_apci13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a p c i 1 3
        **  - - - - - - - - -
        **
        **  Test iauApci13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApci13, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2, eo;
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;

            Sofa.iauApci13(date1, date2, ref astrom, out eo);

            vvd(astrom.pmt, 12.65133794027378508, 1e-11,
                            "iauApci13", "pmt", ref status);
            vvd(astrom.eb[0], 0.9013108747340644755, 1e-12,
                              "iauApci13", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.4174026640406119957, 1e-12,
                              "iauApci13", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.1809822877867817771, 1e-12,
                              "iauApci13", "eb(3)", ref status);
            vvd(astrom.eh[0], 0.8940025429255499549, 1e-12,
                              "iauApci13", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.4110930268331896318, 1e-12,
                              "iauApci13", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.1782189006019749850, 1e-12,
                              "iauApci13", "eh(3)", ref status);
            vvd(astrom.em, 1.010465295964664178, 1e-12,
                           "iauApci13", "em", ref status);
            vvd(astrom.v[0], 0.4289638912941341125e-4, 1e-16,
                             "iauApci13", "v(1)", ref status);
            vvd(astrom.v[1], 0.8115034032405042132e-4, 1e-16,
                             "iauApci13", "v(2)", ref status);
            vvd(astrom.v[2], 0.3517555135536470279e-4, 1e-16,
                             "iauApci13", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999951686013142, 1e-12,
                            "iauApci13", "bm1", ref status);
            vvd(astrom.bpn[0][0], 0.9999992060376761710, 1e-12,
                                  "iauApci13", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0.4124244860106037157e-7, 1e-12,
                                  "iauApci13", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0.1260128571051709670e-2, 1e-12,
                                  "iauApci13", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], -0.1282291987222130690e-7, 1e-12,
                                  "iauApci13", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 0.9999999997456835325, 1e-12,
                                  "iauApci13", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], -0.2255288829420524935e-4, 1e-12,
                                  "iauApci13", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], -0.1260128571661374559e-2, 1e-12,
                                  "iauApci13", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0.2255285422953395494e-4, 1e-12,
                                  "iauApci13", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 0.9999992057833604343, 1e-12,
                                  "iauApci13", "bpn(3,3)", ref status);
            vvd(eo, -0.2900618712657375647e-2, 1e-12,
                    "iauApci13", "eo", ref status);

        }

        static void t_apco(ref int status)
        /*
        **  - - - - - - -
        **   t _ a p c o
        **  - - - - - - -
        **
        **  Test iauApco function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApco, vvd
        **
        **  This revision:  2021 January 5
        */
        {
            double date1, date2, x, y, s,
                theta, elong, phi, hm, xp, yp, sp, refa, refb;
            double[][] ebpv = { new double[3], new double[3] };
            double[] ehp = new double[3];
            iauASTROM astrom = new iauASTROM();


            date1 = 2456384.5;
            date2 = 0.970031644;
            ebpv[0][0] = -0.974170438;
            ebpv[0][1] = -0.211520082;
            ebpv[0][2] = -0.0917583024;
            ebpv[1][0] = 0.00364365824;
            ebpv[1][1] = -0.0154287319;
            ebpv[1][2] = -0.00668922024;
            ehp[0] = -0.973458265;
            ehp[1] = -0.209215307;
            ehp[2] = -0.0906996477;
            x = 0.0013122272;
            y = -2.92808623e-5;
            s = 3.05749468e-8;
            theta = 3.14540971;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            sp = -3.01974337e-11;
            refa = 0.000201418779;
            refb = -2.36140831e-7;

            Sofa.iauApco(date1, date2, ebpv, ehp, x, y, s,
                    theta, elong, phi, hm, xp, yp, sp,
                    refa, refb, ref astrom);

            vvd(astrom.pmt, 13.25248468622587269, 1e-11,
                            "iauApco", "pmt", ref status);
            vvd(astrom.eb[0], -0.9741827110630322720, 1e-12,
                              "iauApco", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.2115130190135344832, 1e-12,
                              "iauApco", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.09179840186949532298, 1e-12,
                              "iauApco", "eb(3)", ref status);
            vvd(astrom.eh[0], -0.9736425571689739035, 1e-12,
                              "iauApco", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.2092452125849330936, 1e-12,
                              "iauApco", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.09075578152243272599, 1e-12,
                              "iauApco", "eh(3)", ref status);
            vvd(astrom.em, 0.9998233241709957653, 1e-12,
                           "iauApco", "em", ref status);
            vvd(astrom.v[0], 0.2078704992916728762e-4, 1e-16,
                             "iauApco", "v(1)", ref status);
            vvd(astrom.v[1], -0.8955360107151952319e-4, 1e-16,
                             "iauApco", "v(2)", ref status);
            vvd(astrom.v[2], -0.3863338994288951082e-4, 1e-16,
                             "iauApco", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999950277561236, 1e-12,
                            "iauApco", "bm1", ref status);
            vvd(astrom.bpn[0][0], 0.9999991390295159156, 1e-12,
                                  "iauApco", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0.4978650072505016932e-7, 1e-12,
                                  "iauApco", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0.1312227200000000000e-2, 1e-12,
                                  "iauApco", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], -0.1136336653771609630e-7, 1e-12,
                                  "iauApco", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 0.9999999995713154868, 1e-12,
                                  "iauApco", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], -0.2928086230000000000e-4, 1e-12,
                                  "iauApco", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], -0.1312227200895260194e-2, 1e-12,
                                  "iauApco", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0.2928082217872315680e-4, 1e-12,
                                  "iauApco", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 0.9999991386008323373, 1e-12,
                                  "iauApco", "bpn(3,3)", ref status);
            vvd(astrom.along, -0.5278008060295995734, 1e-12,
                              "iauApco", "along", ref status);
            vvd(astrom.xpl, 0.1133427418130752958e-5, 1e-17,
                            "iauApco", "xpl", ref status);
            vvd(astrom.ypl, 0.1453347595780646207e-5, 1e-17,
                            "iauApco", "ypl", ref status);
            vvd(astrom.sphi, -0.9440115679003211329, 1e-12,
                             "iauApco", "sphi", ref status);
            vvd(astrom.cphi, 0.3299123514971474711, 1e-12,
                             "iauApco", "cphi", ref status);
            vvd(astrom.diurab, 0, 0,
                               "iauApco", "diurab", ref status);
            vvd(astrom.eral, 2.617608903970400427, 1e-12,
                             "iauApco", "eral", ref status);
            vvd(astrom.refa, 0.2014187790000000000e-3, 1e-15,
                             "iauApco", "refa", ref status);
            vvd(astrom.refb, -0.2361408310000000000e-6, 1e-18,
                             "iauApco", "refb", ref status);

        }

        static void t_apco13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a p c o 1 3
        **  - - - - - - - - -
        **
        **  Test iauApco13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApco13, vvd, viv
        **
        **  This revision:  2021 January 5
        */
        {
            double utc1, utc2, dut1, elong, phi, hm, xp, yp,
                   phpa, tc, rh, wl, eo;
            iauASTROM astrom;
            int j;


            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;

            j = Sofa.iauApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                          phpa, tc, rh, wl, out astrom, out eo);

            vvd(astrom.pmt, 13.25248468622475727, 1e-11,
                            "iauApco13", "pmt", ref status);
            vvd(astrom.eb[0], -0.9741827107320875162, 1e-12,
                            "iauApco13", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.2115130190489716682, 1e-12,
                              "iauApco13", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.09179840189496755339, 1e-12,
                              "iauApco13", "eb(3)", ref status);
            vvd(astrom.eh[0], -0.9736425572586935247, 1e-12,
                              "iauApco13", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.2092452121603336166, 1e-12,
                              "iauApco13", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.09075578153885665295, 1e-12,
                              "iauApco13", "eh(3)", ref status);
            vvd(astrom.em, 0.9998233240913898141, 1e-12,
                           "iauApco13", "em", ref status);
            vvd(astrom.v[0], 0.2078704994520489246e-4, 1e-16,
                             "iauApco13", "v(1)", ref status);
            vvd(astrom.v[1], -0.8955360133238868938e-4, 1e-16,
                             "iauApco13", "v(2)", ref status);
            vvd(astrom.v[2], -0.3863338993055887398e-4, 1e-16,
                             "iauApco13", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999950277561004, 1e-12,
                            "iauApco13", "bm1", ref status);
            vvd(astrom.bpn[0][0], 0.9999991390295147999, 1e-12,
                                  "iauApco13", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0.4978650075315529277e-7, 1e-12,
                                  "iauApco13", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0.001312227200850293372, 1e-12,
                                  "iauApco13", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], -0.1136336652812486604e-7, 1e-12,
                                  "iauApco13", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 0.9999999995713154865, 1e-12,
                                  "iauApco13", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], -0.2928086230975367296e-4, 1e-12,
                                  "iauApco13", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], -0.001312227201745553566, 1e-12,
                                  "iauApco13", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0.2928082218847679162e-4, 1e-12,
                                  "iauApco13", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 0.9999991386008312212, 1e-12,
                                  "iauApco13", "bpn(3,3)", ref status);
            vvd(astrom.along, -0.5278008060295995733, 1e-12,
                              "iauApco13", "along", ref status);
            vvd(astrom.xpl, 0.1133427418130752958e-5, 1e-17,
                            "iauApco13", "xpl", ref status);
            vvd(astrom.ypl, 0.1453347595780646207e-5, 1e-17,
                            "iauApco13", "ypl", ref status);
            vvd(astrom.sphi, -0.9440115679003211329, 1e-12,
                             "iauApco13", "sphi", ref status);
            vvd(astrom.cphi, 0.3299123514971474711, 1e-12,
                             "iauApco13", "cphi", ref status);
            vvd(astrom.diurab, 0, 0,
                               "iauApco13", "diurab", ref status);
            vvd(astrom.eral, 2.617608909189664000, 1e-12,
                             "iauApco13", "eral", ref status);
            vvd(astrom.refa, 0.2014187785940396921e-3, 1e-15,
                             "iauApco13", "refa", ref status);
            vvd(astrom.refb, -0.2361408314943696227e-6, 1e-18,
                             "iauApco13", "refb", ref status);
            vvd(eo, -0.003020548354802412839, 1e-14,
                    "iauApco13", "eo", ref status);
            viv(j, 0, "iauApco13", "j", ref status);

        }

        static void t_apcs(ref int status)
        /*
        **  - - - - - - -
        **   t _ a p c s
        **  - - - - - - -
        **
        **  Test iauApcs function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApcs, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2;
            double[][] pv = { new double[3], new double[3] };
            double[][] ebpv = { new double[3], new double[3] };
            double[] ehp = new double[3];
            iauASTROM astrom = new iauASTROM();


            date1 = 2456384.5;
            date2 = 0.970031644;
            pv[0][0] = -1836024.09;
            pv[0][1] = 1056607.72;
            pv[0][2] = -5998795.26;
            pv[1][0] = -77.0361767;
            pv[1][1] = -133.310856;
            pv[1][2] = 0.0971855934;
            ebpv[0][0] = -0.974170438;
            ebpv[0][1] = -0.211520082;
            ebpv[0][2] = -0.0917583024;
            ebpv[1][0] = 0.00364365824;
            ebpv[1][1] = -0.0154287319;
            ebpv[1][2] = -0.00668922024;
            ehp[0] = -0.973458265;
            ehp[1] = -0.209215307;
            ehp[2] = -0.0906996477;

            Sofa.iauApcs(date1, date2, pv, ebpv, ehp, ref astrom);

            vvd(astrom.pmt, 13.25248468622587269, 1e-11,
                            "iauApcs", "pmt", ref status);
            vvd(astrom.eb[0], -0.9741827110629881886, 1e-12,
                              "iauApcs", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.2115130190136415986, 1e-12,
                              "iauApcs", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.09179840186954412099, 1e-12,
                              "iauApcs", "eb(3)", ref status);
            vvd(astrom.eh[0], -0.9736425571689454706, 1e-12,
                              "iauApcs", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.2092452125850435930, 1e-12,
                              "iauApcs", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.09075578152248299218, 1e-12,
                              "iauApcs", "eh(3)", ref status);
            vvd(astrom.em, 0.9998233241709796859, 1e-12,
                           "iauApcs", "em", ref status);
            vvd(astrom.v[0], 0.2078704993282685510e-4, 1e-16,
                             "iauApcs", "v(1)", ref status);
            vvd(astrom.v[1], -0.8955360106989405683e-4, 1e-16,
                             "iauApcs", "v(2)", ref status);
            vvd(astrom.v[2], -0.3863338994289409097e-4, 1e-16,
                             "iauApcs", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999950277561237, 1e-12,
                            "iauApcs", "bm1", ref status);
            vvd(astrom.bpn[0][0], 1, 0,
                                  "iauApcs", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0, 0,
                                  "iauApcs", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0, 0,
                                  "iauApcs", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], 0, 0,
                                  "iauApcs", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 1, 0,
                                  "iauApcs", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], 0, 0,
                                  "iauApcs", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], 0, 0,
                                  "iauApcs", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0, 0,
                                  "iauApcs", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 1, 0,
                                  "iauApcs", "bpn(3,3)", ref status);

        }

        static void t_apcs13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a p c s 1 3
        **  - - - - - - - - -
        **
        **  Test iauApcs13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApcs13, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2;
            double[][] pv = { new double[3], new double[3] };
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;
            pv[0][0] = -6241497.16;
            pv[0][1] = 401346.896;
            pv[0][2] = -1251136.04;
            pv[1][0] = -29.264597;
            pv[1][1] = -455.021831;
            pv[1][2] = 0.0266151194;

            Sofa.iauApcs13(date1, date2, pv, ref astrom);

            vvd(astrom.pmt, 12.65133794027378508, 1e-11,
                            "iauApcs13", "pmt", ref status);
            vvd(astrom.eb[0], 0.9012691529025250644, 1e-12,
                              "iauApcs13", "eb(1)", ref status);
            vvd(astrom.eb[1], -0.4173999812023194317, 1e-12,
                              "iauApcs13", "eb(2)", ref status);
            vvd(astrom.eb[2], -0.1809906511146429670, 1e-12,
                              "iauApcs13", "eb(3)", ref status);
            vvd(astrom.eh[0], 0.8939939101760130792, 1e-12,
                              "iauApcs13", "eh(1)", ref status);
            vvd(astrom.eh[1], -0.4111053891734021478, 1e-12,
                              "iauApcs13", "eh(2)", ref status);
            vvd(astrom.eh[2], -0.1782336880636997374, 1e-12,
                              "iauApcs13", "eh(3)", ref status);
            vvd(astrom.em, 1.010428384373491095, 1e-12,
                           "iauApcs13", "em", ref status);
            vvd(astrom.v[0], 0.4279877294121697570e-4, 1e-16,
                             "iauApcs13", "v(1)", ref status);
            vvd(astrom.v[1], 0.7963255087052120678e-4, 1e-16,
                             "iauApcs13", "v(2)", ref status);
            vvd(astrom.v[2], 0.3517564013384691531e-4, 1e-16,
                             "iauApcs13", "v(3)", ref status);
            vvd(astrom.bm1, 0.9999999952947980978, 1e-12,
                            "iauApcs13", "bm1", ref status);
            vvd(astrom.bpn[0][0], 1, 0,
                                  "iauApcs13", "bpn(1,1)", ref status);
            vvd(astrom.bpn[1][0], 0, 0,
                                  "iauApcs13", "bpn(2,1)", ref status);
            vvd(astrom.bpn[2][0], 0, 0,
                                  "iauApcs13", "bpn(3,1)", ref status);
            vvd(astrom.bpn[0][1], 0, 0,
                                  "iauApcs13", "bpn(1,2)", ref status);
            vvd(astrom.bpn[1][1], 1, 0,
                                  "iauApcs13", "bpn(2,2)", ref status);
            vvd(astrom.bpn[2][1], 0, 0,
                                  "iauApcs13", "bpn(3,2)", ref status);
            vvd(astrom.bpn[0][2], 0, 0,
                                  "iauApcs13", "bpn(1,3)", ref status);
            vvd(astrom.bpn[1][2], 0, 0,
                                  "iauApcs13", "bpn(2,3)", ref status);
            vvd(astrom.bpn[2][2], 1, 0,
                                  "iauApcs13", "bpn(3,3)", ref status);

        }

        static void t_aper(ref int status)
        /*
        **  - - - - - - -
        **   t _ a p e r
        **  - - - - - - -
        *
        **  Test iauAper function.
        *
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        *
        **  Called:  iauAper, vvd
        *
        **  This revision:  2013 October 3
        */
        {
            double theta;
            iauASTROM astrom = new iauASTROM();


            astrom.along = 1.234;
            theta = 5.678;

            Sofa.iauAper(theta, astrom);

            vvd(astrom.eral, 6.912000000000000000, 1e-12,
                             "iauAper", "pmt", ref status);

        }

        static void t_aper13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a p e r 1 3
        **  - - - - - - - - -
        **
        **  Test iauAper13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAper13, vvd
        **
        **  This revision:  2013 October 3
        */
        {
            double ut11, ut12;
            iauASTROM astrom = new iauASTROM();


            astrom.along = 1.234;
            ut11 = 2456165.5;
            ut12 = 0.401182685;

            Sofa.iauAper13(ut11, ut12, astrom);

            vvd(astrom.eral, 3.316236661789694933, 1e-12,
                             "iauAper13", "pmt", ref status);

        }

        static void t_apio(ref int status)
        /*
        **  - - - - - - -
        **   t _ a p i o
        **  - - - - - - -
        **
        **  Test iauApio function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApio, vvd
        **
        **  This revision:  2021 January 5
        */
        {
            double sp, theta, elong, phi, hm, xp, yp, refa, refb;
            iauASTROM astrom;


            sp = -3.01974337e-11;
            theta = 3.14540971;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            refa = 0.000201418779;
            refb = -2.36140831e-7;

            Sofa.iauApio(sp, theta, elong, phi, hm, xp, yp, refa, refb, out astrom);

            vvd(astrom.along, -0.5278008060295995734, 1e-12,
                              "iauApio", "along", ref status);
            vvd(astrom.xpl, 0.1133427418130752958e-5, 1e-17,
                            "iauApio", "xpl", ref status);
            vvd(astrom.ypl, 0.1453347595780646207e-5, 1e-17,
                            "iauApio", "ypl", ref status);
            vvd(astrom.sphi, -0.9440115679003211329, 1e-12,
                             "iauApio", "sphi", ref status);
            vvd(astrom.cphi, 0.3299123514971474711, 1e-12,
                             "iauApio", "cphi", ref status);
            vvd(astrom.diurab, 0.5135843661699913529e-6, 1e-12,
                               "iauApio", "diurab", ref status);
            vvd(astrom.eral, 2.617608903970400427, 1e-12,
                             "iauApio", "eral", ref status);
            vvd(astrom.refa, 0.2014187790000000000e-3, 1e-15,
                             "iauApio", "refa", ref status);
            vvd(astrom.refb, -0.2361408310000000000e-6, 1e-18,
                             "iauApio", "refb", ref status);

        }

        static void t_apio13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a p i o 1 3
        **  - - - - - - - - -
        **
        **  Test iauApio13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApio13, vvd, viv
        **
        **  This revision:  2021 January 5
        */
        {
            double utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl;
            int j;
            iauASTROM astrom;


            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;

            j = Sofa.iauApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                          phpa, tc, rh, wl, out astrom);

            vvd(astrom.along, -0.5278008060295995733, 1e-12,
                              "iauApio13", "along", ref status);
            vvd(astrom.xpl, 0.1133427418130752958e-5, 1e-17,
                            "iauApio13", "xpl", ref status);
            vvd(astrom.ypl, 0.1453347595780646207e-5, 1e-17,
                            "iauApio13", "ypl", ref status);
            vvd(astrom.sphi, -0.9440115679003211329, 1e-12,
                             "iauApio13", "sphi", ref status);
            vvd(astrom.cphi, 0.3299123514971474711, 1e-12,
                             "iauApio13", "cphi", ref status);
            vvd(astrom.diurab, 0.5135843661699913529e-6, 1e-12,
                               "iauApio13", "diurab", ref status);
            vvd(astrom.eral, 2.617608909189664000, 1e-12,
                             "iauApio13", "eral", ref status);
            vvd(astrom.refa, 0.2014187785940396921e-3, 1e-15,
                             "iauApio13", "refa", ref status);
            vvd(astrom.refb, -0.2361408314943696227e-6, 1e-18,
                             "iauApio13", "refb", ref status);
            viv(j, 0, "iauApio13", "j", ref status);

        }

        static void t_atcc13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t c c 1 3
        **  - - - - - - - - -
        **
        **  Test iauAtcc13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAtcc13, vvd
        **
        **  This revision:  2021 April 18
        */
        {
            double rc, dc, pr, pd, px, rv, date1, date2, ra, da;


            rc = 2.71;
            dc = 0.174;
            pr = 1e-5;
            pd = 5e-6;
            px = 0.1;
            rv = 55.0;
            date1 = 2456165.5;
            date2 = 0.401182685;

            Sofa.iauAtcc13(rc, dc, pr, pd, px, rv, date1, date2, out ra, out da);

            vvd(ra, 2.710126504531372384, 1e-12,
                    "iauAtcc13", "ra", ref status);
            vvd(da, 0.1740632537628350152, 1e-12,
                    "iauAtcc13", "da", ref status);

        }

        static void t_atccq(ref int status)
        /*
        **  - - - - - - - -
        **   t _ a t c c q
        **  - - - - - - - -
        **
        **  Test iauAtccq function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApcc13, iauAtccq, vvd
        **
        **  This revision:  2021 April 18
        */
        {
            double date1, date2, eo, rc, dc, pr, pd, px, rv, ra, da;
            iauASTROM astrom = new iauASTROM();

            date1 = 2456165.5;
            date2 = 0.401182685;
            Sofa.iauApci13(date1, date2, ref astrom, out eo);
            rc = 2.71;
            dc = 0.174;
            pr = 1e-5;
            pd = 5e-6;
            px = 0.1;
            rv = 55.0;

            Sofa.iauAtccq(rc, dc, pr, pd, px, rv, astrom, out ra, out da);

            vvd(ra, 2.710126504531372384, 1e-12, "iauAtccq", "ra", ref status);
            vvd(da, 0.1740632537628350152, 1e-12, "iauAtccq", "da", ref status);

        }

        static void t_atci13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t c i 1 3
        **  - - - - - - - - -
        **
        **  Test iauAtci13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAtci13, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double rc, dc, pr, pd, px, rv, date1, date2, ri, di, eo;


            rc = 2.71;
            dc = 0.174;
            pr = 1e-5;
            pd = 5e-6;
            px = 0.1;
            rv = 55.0;
            date1 = 2456165.5;
            date2 = 0.401182685;

            Sofa.iauAtci13(rc, dc, pr, pd, px, rv, date1, date2, out ri, out di, out eo);

            vvd(ri, 2.710121572968696744, 1e-12,
                    "iauAtci13", "ri", ref status);
            vvd(di, 0.1729371367219539137, 1e-12,
                    "iauAtci13", "di", ref status);
            vvd(eo, -0.002900618712657375647, 1e-14,
                    "iauAtci13", "eo", ref status);

        }

        static void t_atciq(ref int status)
        /*
        **  - - - - - - - -
        **   t _ a t c i q
        **  - - - - - - - -
        **
        **  Test iauAtciq function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApci13, iauAtciq, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2, eo, rc, dc, pr, pd, px, rv, ri, di;
            iauASTROM astrom = new iauASTROM();

            date1 = 2456165.5;
            date2 = 0.401182685;
            Sofa.iauApci13(date1, date2, ref astrom, out eo);
            rc = 2.71;
            dc = 0.174;
            pr = 1e-5;
            pd = 5e-6;
            px = 0.1;
            rv = 55.0;

            Sofa.iauAtciq(rc, dc, pr, pd, px, rv, astrom, out ri, out di);

            vvd(ri, 2.710121572968696744, 1e-12, "iauAtciq", "ri", ref status);
            vvd(di, 0.1729371367219539137, 1e-12, "iauAtciq", "di", ref status);

        }

        static void t_atciqn(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t c i q n
        **  - - - - - - - - -
        **
        **  Test iauAtciqn function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApci13, iauAtciqn, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            iauLDBODY[] b = new iauLDBODY[] { new iauLDBODY(), new iauLDBODY(), new iauLDBODY() };
            double date1, date2, eo, rc, dc, pr, pd, px, rv, ri, di;
            iauASTROM astrom = new iauASTROM();

            date1 = 2456165.5;
            date2 = 0.401182685;
            Sofa.iauApci13(date1, date2, ref astrom, out eo);
            rc = 2.71;
            dc = 0.174;
            pr = 1e-5;
            pd = 5e-6;
            px = 0.1;
            rv = 55.0;
            b[0].bm = 0.00028574;
            b[0].dl = 3e-10;
            b[0].pv[0][0] = -7.81014427;
            b[0].pv[0][1] = -5.60956681;
            b[0].pv[0][2] = -1.98079819;
            b[0].pv[1][0] = 0.0030723249;
            b[0].pv[1][1] = -0.00406995477;
            b[0].pv[1][2] = -0.00181335842;
            b[1].bm = 0.00095435;
            b[1].dl = 3e-9;
            b[1].pv[0][0] = 0.738098796;
            b[1].pv[0][1] = 4.63658692;
            b[1].pv[0][2] = 1.9693136;
            b[1].pv[1][0] = -0.00755816922;
            b[1].pv[1][1] = 0.00126913722;
            b[1].pv[1][2] = 0.000727999001;
            b[2].bm = 1.0;
            b[2].dl = 6e-6;
            b[2].pv[0][0] = -0.000712174377;
            b[2].pv[0][1] = -0.00230478303;
            b[2].pv[0][2] = -0.00105865966;
            b[2].pv[1][0] = 6.29235213e-6;
            b[2].pv[1][1] = -3.30888387e-7;
            b[2].pv[1][2] = -2.96486623e-7;

            Sofa.iauAtciqn(rc, dc, pr, pd, px, rv, astrom, 3, b, out ri, out di);

            vvd(ri, 2.710122008104983335, 1e-12, "iauAtciqn", "ri", ref status);
            vvd(di, 0.1729371916492767821, 1e-12, "iauAtciqn", "di", ref status);

        }

        static void t_atciqz(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t c i q z
        **  - - - - - - - - -
        **
        **  Test iauAtciqz function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApci13, iauAtciqz, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2, eo, rc, dc, ri, di;
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;
            Sofa.iauApci13(date1, date2, ref astrom, out eo);
            rc = 2.71;
            dc = 0.174;

            Sofa.iauAtciqz(rc, dc, astrom, out ri, out di);

            vvd(ri, 2.709994899247256984, 1e-12, "iauAtciqz", "ri", ref status);
            vvd(di, 0.1728740720984931891, 1e-12, "iauAtciqz", "di", ref status);

        }

        static void t_atco13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t c o 1 3
        **  - - - - - - - - -
        **
        **  Test iauAtco13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAtco13, vvd, viv
        **
        **  This revision:  2021 January 5
        */
        {
            double rc, dc, pr, pd, px, rv, utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                   aob, zob, hob, dob, rob, eo;
            int j;


            rc = 2.71;
            dc = 0.174;
            pr = 1e-5;
            pd = 5e-6;
            px = 0.1;
            rv = 55.0;
            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;

            j = Sofa.iauAtco13(rc, dc, pr, pd, px, rv,
                          utc1, utc2, dut1, elong, phi, hm, xp, yp,
                          phpa, tc, rh, wl,
                          out aob, out zob, out hob, out dob, out rob, out eo);

            vvd(aob, 0.9251774485485515207e-1, 1e-12, "iauAtco13", "aob", ref status);
            vvd(zob, 1.407661405256499357, 1e-12, "iauAtco13", "zob", ref status);
            vvd(hob, -0.9265154431529724692e-1, 1e-12, "iauAtco13", "hob", ref status);
            vvd(dob, 0.1716626560072526200, 1e-12, "iauAtco13", "dob", ref status);
            vvd(rob, 2.710260453504961012, 1e-12, "iauAtco13", "rob", ref status);
            vvd(eo, -0.003020548354802412839, 1e-14, "iauAtco13", "eo", ref status);
            viv(j, 0, "iauAtco13", "j", ref status);

        }

        static void t_atic13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t i c 1 3
        **  - - - - - - - - -
        **
        **  Test iauAtic13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAtic13, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double ri, di, date1, date2, rc, dc, eo;


            ri = 2.710121572969038991;
            di = 0.1729371367218230438;
            date1 = 2456165.5;
            date2 = 0.401182685;

            Sofa.iauAtic13(ri, di, date1, date2, out rc, out dc, out eo);

            vvd(rc, 2.710126504531716819, 1e-12, "iauAtic13", "rc", ref status);
            vvd(dc, 0.1740632537627034482, 1e-12, "iauAtic13", "dc", ref status);
            vvd(eo, -0.002900618712657375647, 1e-14, "iauAtic13", "eo", ref status);

        }

        static void t_aticq(ref int status)
        /*
        **  - - - - - - - -
        **   t _ a t i c q
        **  - - - - - - - -
        **
        **  Test iauAticq function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApci13, iauAticq, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2, eo, ri, di, rc, dc;
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;
            Sofa.iauApci13(date1, date2, ref astrom, out eo);
            ri = 2.710121572969038991;
            di = 0.1729371367218230438;

            Sofa.iauAticq(ri, di, astrom, out rc, out dc);

            vvd(rc, 2.710126504531716819, 1e-12, "iauAticq", "rc", ref status);
            vvd(dc, 0.1740632537627034482, 1e-12, "iauAticq", "dc", ref status);

        }

        static void t_aticqn(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t i c q n
        **  - - - - - - - - -
        **
        **  Test iauAticqn function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApci13, iauAticqn, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double date1, date2, eo, ri, di, rc, dc;
            iauLDBODY[] b = new iauLDBODY[] { new iauLDBODY(), new iauLDBODY(), new iauLDBODY() };
            iauASTROM astrom = new iauASTROM();


            date1 = 2456165.5;
            date2 = 0.401182685;
            Sofa.iauApci13(date1, date2, ref astrom, out eo);
            ri = 2.709994899247599271;
            di = 0.1728740720983623469;
            b[0].bm = 0.00028574;
            b[0].dl = 3e-10;
            b[0].pv[0][0] = -7.81014427;
            b[0].pv[0][1] = -5.60956681;
            b[0].pv[0][2] = -1.98079819;
            b[0].pv[1][0] = 0.0030723249;
            b[0].pv[1][1] = -0.00406995477;
            b[0].pv[1][2] = -0.00181335842;
            b[1].bm = 0.00095435;
            b[1].dl = 3e-9;
            b[1].pv[0][0] = 0.738098796;
            b[1].pv[0][1] = 4.63658692;
            b[1].pv[0][2] = 1.9693136;
            b[1].pv[1][0] = -0.00755816922;
            b[1].pv[1][1] = 0.00126913722;
            b[1].pv[1][2] = 0.000727999001;
            b[2].bm = 1.0;
            b[2].dl = 6e-6;
            b[2].pv[0][0] = -0.000712174377;
            b[2].pv[0][1] = -0.00230478303;
            b[2].pv[0][2] = -0.00105865966;
            b[2].pv[1][0] = 6.29235213e-6;
            b[2].pv[1][1] = -3.30888387e-7;
            b[2].pv[1][2] = -2.96486623e-7;

            Sofa.iauAticqn(ri, di, astrom, 3, b, out rc, out dc);

            vvd(rc, 2.709999575033027333, 1e-12, "iauAtciqn", "rc", ref status);
            vvd(dc, 0.1739999656316469990, 1e-12, "iauAtciqn", "dc", ref status);

        }

        static void t_atio13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t i o 1 3
        **  - - - - - - - - -
        **
        **  Test iauAtio13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAtio13, vvd, viv
        **
        **  This revision:  2021 January 5
        */
        {
            double ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp,
                   phpa, tc, rh, wl, aob, zob, hob, dob, rob;
            int j;


            ri = 2.710121572969038991;
            di = 0.1729371367218230438;
            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;

            j = Sofa.iauAtio13(ri, di, utc1, utc2, dut1, elong, phi, hm,
                          xp, yp, phpa, tc, rh, wl,
                          out aob, out zob, out hob, out dob, out rob);

            vvd(aob, 0.9233952224895122499e-1, 1e-12, "iauAtio13", "aob", ref status);
            vvd(zob, 1.407758704513549991, 1e-12, "iauAtio13", "zob", ref status);
            vvd(hob, -0.9247619879881698140e-1, 1e-12, "iauAtio13", "hob", ref status);
            vvd(dob, 0.1717653435756234676, 1e-12, "iauAtio13", "dob", ref status);
            vvd(rob, 2.710085107988480746, 1e-12, "iauAtio13", "rob", ref status);
            viv(j, 0, "iauAtio13", "j", ref status);

        }

        static void t_atioq(ref int status)
        /*
        **  - - - - - - - -
        **   t _ a t i o q
        **  - - - - - - - -
        **
        **  Test iauAtioq function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauApio13, iauAtioq, vvd, viv
        **
        **  This revision:  2021 January 5
        */
        {
            double utc1, utc2, dut1, elong, phi, hm, xp, yp,
                   phpa, tc, rh, wl, ri, di, aob, zob, hob, dob, rob;
            iauASTROM astrom;


            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;
            Sofa.iauApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                             phpa, tc, rh, wl, out astrom);
            ri = 2.710121572969038991;
            di = 0.1729371367218230438;

            Sofa.iauAtioq(ri, di, astrom, out aob, out zob, out hob, out dob, out rob);

            vvd(aob, 0.9233952224895122499e-1, 1e-12, "iauAtioq", "aob", ref status);
            vvd(zob, 1.407758704513549991, 1e-12, "iauAtioq", "zob", ref status);
            vvd(hob, -0.9247619879881698140e-1, 1e-12, "iauAtioq", "hob", ref status);
            vvd(dob, 0.1717653435756234676, 1e-12, "iauAtioq", "dob", ref status);
            vvd(rob, 2.710085107988480746, 1e-12, "iauAtioq", "rob", ref status);

        }

        static void t_atoc13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t o c 1 3
        **  - - - - - - - - -
        **
        **  Test iauAtoc13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAtoc13, vvd, viv
        **
        **  This revision:  2021 January 5
        */
        {
            double utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                   ob1, ob2, rc, dc;
            int j;


            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;

            ob1 = 2.710085107986886201;
            ob2 = 0.1717653435758265198;
            j = Sofa.iauAtoc13("R", ob1, ob2, utc1, utc2, dut1,
                            elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                            out rc, out dc);
            vvd(rc, 2.709956744659136129, 1e-12, "iauAtoc13", "R/rc", ref status);
            vvd(dc, 0.1741696500898471362, 1e-12, "iauAtoc13", "R/dc", ref status);
            viv(j, 0, "iauAtoc13", "R/j", ref status);

            ob1 = -0.09247619879782006106;
            ob2 = 0.1717653435758265198;
            j = Sofa.iauAtoc13("H", ob1, ob2, utc1, utc2, dut1,
                            elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                            out rc, out dc);
            vvd(rc, 2.709956744659734086, 1e-12, "iauAtoc13", "H/rc", ref status);
            vvd(dc, 0.1741696500898471362, 1e-12, "iauAtoc13", "H/dc", ref status);
            viv(j, 0, "iauAtoc13", "H/j", ref status);

            ob1 = 0.09233952224794989993;
            ob2 = 1.407758704513722461;
            j = Sofa.iauAtoc13("A", ob1, ob2, utc1, utc2, dut1,
                            elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                            out rc, out dc);
            vvd(rc, 2.709956744659734086, 1e-12, "iauAtoc13", "A/rc", ref status);
            vvd(dc, 0.1741696500898471366, 1e-12, "iauAtoc13", "A/dc", ref status);
            viv(j, 0, "iauAtoc13", "A/j", ref status);

        }

        static void t_atoi13(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ a t o i 1 3
        **  - - - - - - - - -
        **
        **  Test iauAtoi13 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauAtoi13, vvd, viv
        **
        **  This revision:  2021 January 5
        */
        {
            double utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                   ob1, ob2, ri, di;
            int j;


            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;

            ob1 = 2.710085107986886201;
            ob2 = 0.1717653435758265198;
            j = Sofa.iauAtoi13("R", ob1, ob2, utc1, utc2, dut1,
                            elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                            out ri, out di);
            vvd(ri, 2.710121574447540810, 1e-12, "iauAtoi13", "R/ri", ref status);
            vvd(di, 0.1729371839116608778, 1e-12, "iauAtoi13", "R/di", ref status);
            viv(j, 0, "iauAtoi13", "R/J", ref status);

            ob1 = -0.09247619879782006106;
            ob2 = 0.1717653435758265198;
            j = Sofa.iauAtoi13("H", ob1, ob2, utc1, utc2, dut1,
                            elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                            out ri, out di);
            vvd(ri, 2.710121574448138676, 1e-12, "iauAtoi13", "H/ri", ref status);
            vvd(di, 0.1729371839116608778, 1e-12, "iauAtoi13", "H/di", ref status);
            viv(j, 0, "iauAtoi13", "H/J", ref status);

            ob1 = 0.09233952224794989993;
            ob2 = 1.407758704513722461;
            j = Sofa.iauAtoi13("A", ob1, ob2, utc1, utc2, dut1,
                            elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                            out ri, out di);
            vvd(ri, 2.710121574448138676, 1e-12, "iauAtoi13", "A/ri", ref status);
            vvd(di, 0.1729371839116608781, 1e-12, "iauAtoi13", "A/di", ref status);
            viv(j, 0, "iauAtoi13", "A/J", ref status);

        }

        static void t_atoiq(ref int status)
        /*
        **  - - - - - - - -
        **   t _ a t o i q
        **  - - - - - - - -
        *
        **  Test iauAtoiq function.
        *
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        *
        **  Called:  iauApio13, iauAtoiq, vvd
        *
        **  This revision:  2021 January 5
        */
        {
            double utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
                   ob1, ob2, ri, di;
            iauASTROM astrom;


            utc1 = 2456384.5;
            utc2 = 0.969254051;
            dut1 = 0.1550675;
            elong = -0.527800806;
            phi = -1.2345856;
            hm = 2738.0;
            xp = 2.47230737e-7;
            yp = 1.82640464e-6;
            phpa = 731.0;
            tc = 12.8;
            rh = 0.59;
            wl = 0.55;
            Sofa.iauApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                             phpa, tc, rh, wl, out astrom);

            ob1 = 2.710085107986886201;
            ob2 = 0.1717653435758265198;
            Sofa.iauAtoiq("R", ob1, ob2, astrom, out ri, out di);
            vvd(ri, 2.710121574447540810, 1e-12,
                    "iauAtoiq", "R/ri", ref status);
            vvd(di, 0.17293718391166087785, 1e-12,
                    "iauAtoiq", "R/di", ref status);

            ob1 = -0.09247619879782006106;
            ob2 = 0.1717653435758265198;
            Sofa.iauAtoiq("H", ob1, ob2, astrom, out ri, out di);
            vvd(ri, 2.710121574448138676, 1e-12,
                    "iauAtoiq", "H/ri", ref status);
            vvd(di, 0.1729371839116608778, 1e-12,
                    "iauAtoiq", "H/di", ref status);

            ob1 = 0.09233952224794989993;
            ob2 = 1.407758704513722461;
            Sofa.iauAtoiq("A", ob1, ob2, astrom, out ri, out di);
            vvd(ri, 2.710121574448138676, 1e-12,
                    "iauAtoiq", "A/ri", ref status);
            vvd(di, 0.1729371839116608781, 1e-12,
                    "iauAtoiq", "A/di", ref status);

        }

        static void t_bi00(ref int status)
        /*
        **  - - - - - - -
        **   t _ b i 0 0
        **  - - - - - - -
        **
        **  Test iauBi00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauBi00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsibi, depsbi, dra;

            Sofa.iauBi00(out dpsibi, out depsbi, out dra);

            vvd(dpsibi, -0.2025309152835086613e-6, 1e-12,
               "iauBi00", "dpsibi", ref status);
            vvd(depsbi, -0.3306041454222147847e-7, 1e-12,
               "iauBi00", "depsbi", ref status);
            vvd(dra, -0.7078279744199225506e-7, 1e-12,
               "iauBi00", "dra", ref status);
        }

        static void t_bp00(ref int status)
        /*
        **  - - - - - - -
        **   t _ b p 0 0
        **  - - - - - - -
        **
        **  Test iauBp00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauBp00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rb, rp, rbp;


            Sofa.iauBp00(2400000.5, 50123.9999, out rb, out rp, out rbp);

            vvd(rb[0][0], 0.9999999999999942498, 1e-12,
                "iauBp00", "rb11", ref status);
            vvd(rb[0][1], -0.7078279744199196626e-7, 1e-16,
                "iauBp00", "rb12", ref status);
            vvd(rb[0][2], 0.8056217146976134152e-7, 1e-16,
                "iauBp00", "rb13", ref status);
            vvd(rb[1][0], 0.7078279477857337206e-7, 1e-16,
                "iauBp00", "rb21", ref status);
            vvd(rb[1][1], 0.9999999999999969484, 1e-12,
                "iauBp00", "rb22", ref status);
            vvd(rb[1][2], 0.3306041454222136517e-7, 1e-16,
                "iauBp00", "rb23", ref status);
            vvd(rb[2][0], -0.8056217380986972157e-7, 1e-16,
                "iauBp00", "rb31", ref status);
            vvd(rb[2][1], -0.3306040883980552500e-7, 1e-16,
                "iauBp00", "rb32", ref status);
            vvd(rb[2][2], 0.9999999999999962084, 1e-12,
                "iauBp00", "rb33", ref status);

            vvd(rp[0][0], 0.9999995504864048241, 1e-12,
                "iauBp00", "rp11", ref status);
            vvd(rp[0][1], 0.8696113836207084411e-3, 1e-14,
                "iauBp00", "rp12", ref status);
            vvd(rp[0][2], 0.3778928813389333402e-3, 1e-14,
                "iauBp00", "rp13", ref status);
            vvd(rp[1][0], -0.8696113818227265968e-3, 1e-14,
                "iauBp00", "rp21", ref status);
            vvd(rp[1][1], 0.9999996218879365258, 1e-12,
                "iauBp00", "rp22", ref status);
            vvd(rp[1][2], -0.1690679263009242066e-6, 1e-14,
                "iauBp00", "rp23", ref status);
            vvd(rp[2][0], -0.3778928854764695214e-3, 1e-14,
                "iauBp00", "rp31", ref status);
            vvd(rp[2][1], -0.1595521004195286491e-6, 1e-14,
                "iauBp00", "rp32", ref status);
            vvd(rp[2][2], 0.9999999285984682756, 1e-12,
                "iauBp00", "rp33", ref status);

            vvd(rbp[0][0], 0.9999995505175087260, 1e-12,
                "iauBp00", "rbp11", ref status);
            vvd(rbp[0][1], 0.8695405883617884705e-3, 1e-14,
                "iauBp00", "rbp12", ref status);
            vvd(rbp[0][2], 0.3779734722239007105e-3, 1e-14,
                "iauBp00", "rbp13", ref status);
            vvd(rbp[1][0], -0.8695405990410863719e-3, 1e-14,
                "iauBp00", "rbp21", ref status);
            vvd(rbp[1][1], 0.9999996219494925900, 1e-12,
                "iauBp00", "rbp22", ref status);
            vvd(rbp[1][2], -0.1360775820404982209e-6, 1e-14,
                "iauBp00", "rbp23", ref status);
            vvd(rbp[2][0], -0.3779734476558184991e-3, 1e-14,
                "iauBp00", "rbp31", ref status);
            vvd(rbp[2][1], -0.1925857585832024058e-6, 1e-14,
                "iauBp00", "rbp32", ref status);
            vvd(rbp[2][2], 0.9999999285680153377, 1e-12,
                "iauBp00", "rbp33", ref status);
        }

        static void t_bp06(ref int status)
        /*
        **  - - - - - - -
        **   t _ b p 0 6
        **  - - - - - - -
        **
        **  Test iauBp06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauBp06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rb, rp, rbp;


            Sofa.iauBp06(2400000.5, 50123.9999, out rb, out rp, out rbp);

            vvd(rb[0][0], 0.9999999999999942497, 1e-12,
                "iauBp06", "rb11", ref status);
            vvd(rb[0][1], -0.7078368960971557145e-7, 1e-14,
                "iauBp06", "rb12", ref status);
            vvd(rb[0][2], 0.8056213977613185606e-7, 1e-14,
                "iauBp06", "rb13", ref status);
            vvd(rb[1][0], 0.7078368694637674333e-7, 1e-14,
                "iauBp06", "rb21", ref status);
            vvd(rb[1][1], 0.9999999999999969484, 1e-12,
                "iauBp06", "rb22", ref status);
            vvd(rb[1][2], 0.3305943742989134124e-7, 1e-14,
                "iauBp06", "rb23", ref status);
            vvd(rb[2][0], -0.8056214211620056792e-7, 1e-14,
                "iauBp06", "rb31", ref status);
            vvd(rb[2][1], -0.3305943172740586950e-7, 1e-14,
                "iauBp06", "rb32", ref status);
            vvd(rb[2][2], 0.9999999999999962084, 1e-12,
                "iauBp06", "rb33", ref status);

            vvd(rp[0][0], 0.9999995504864960278, 1e-12,
                "iauBp06", "rp11", ref status);
            vvd(rp[0][1], 0.8696112578855404832e-3, 1e-14,
                "iauBp06", "rp12", ref status);
            vvd(rp[0][2], 0.3778929293341390127e-3, 1e-14,
                "iauBp06", "rp13", ref status);
            vvd(rp[1][0], -0.8696112560510186244e-3, 1e-14,
                "iauBp06", "rp21", ref status);
            vvd(rp[1][1], 0.9999996218880458820, 1e-12,
                "iauBp06", "rp22", ref status);
            vvd(rp[1][2], -0.1691646168941896285e-6, 1e-14,
                "iauBp06", "rp23", ref status);
            vvd(rp[2][0], -0.3778929335557603418e-3, 1e-14,
                "iauBp06", "rp31", ref status);
            vvd(rp[2][1], -0.1594554040786495076e-6, 1e-14,
                "iauBp06", "rp32", ref status);
            vvd(rp[2][2], 0.9999999285984501222, 1e-12,
                "iauBp06", "rp33", ref status);

            vvd(rbp[0][0], 0.9999995505176007047, 1e-12,
                "iauBp06", "rbp11", ref status);
            vvd(rbp[0][1], 0.8695404617348208406e-3, 1e-14,
                "iauBp06", "rbp12", ref status);
            vvd(rbp[0][2], 0.3779735201865589104e-3, 1e-14,
                "iauBp06", "rbp13", ref status);
            vvd(rbp[1][0], -0.8695404723772031414e-3, 1e-14,
                "iauBp06", "rbp21", ref status);
            vvd(rbp[1][1], 0.9999996219496027161, 1e-12,
                "iauBp06", "rbp22", ref status);
            vvd(rbp[1][2], -0.1361752497080270143e-6, 1e-14,
                "iauBp06", "rbp23", ref status);
            vvd(rbp[2][0], -0.3779734957034089490e-3, 1e-14,
                "iauBp06", "rbp31", ref status);
            vvd(rbp[2][1], -0.1924880847894457113e-6, 1e-14,
                "iauBp06", "rbp32", ref status);
            vvd(rbp[2][2], 0.9999999285679971958, 1e-12,
                "iauBp06", "rbp33", ref status);
        }

        static void t_bpn2xy(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ b p n 2 x y
        **  - - - - - - - - -
        **
        **  Test iauBpn2xy function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauBpn2xy, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y;
            double[][] rbpn = { new double[3], new double[3], new double[3] };


            rbpn[0][0] = 9.999962358680738e-1;
            rbpn[0][1] = -2.516417057665452e-3;
            rbpn[0][2] = -1.093569785342370e-3;

            rbpn[1][0] = 2.516462370370876e-3;
            rbpn[1][1] = 9.999968329010883e-1;
            rbpn[1][2] = 4.006159587358310e-5;

            rbpn[2][0] = 1.093465510215479e-3;
            rbpn[2][1] = -4.281337229063151e-5;
            rbpn[2][2] = 9.999994012499173e-1;

            Sofa.iauBpn2xy(rbpn, out x, out y);

            vvd(x, 1.093465510215479e-3, 1e-12, "iauBpn2xy", "x", ref status);
            vvd(y, -4.281337229063151e-5, 1e-12, "iauBpn2xy", "y", ref status);

        }

        static void t_c2i00a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 i 0 0 a
        **  - - - - - - - - -
        **
        **  Test iauC2i00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2i00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rc2i;


            Sofa.iauC2i00a(2400000.5, 53736.0, out rc2i);

            vvd(rc2i[0][0], 0.9999998323037165557, 1e-12,
                "iauC2i00a", "11", ref status);
            vvd(rc2i[0][1], 0.5581526348992140183e-9, 1e-12,
                "iauC2i00a", "12", ref status);
            vvd(rc2i[0][2], -0.5791308477073443415e-3, 1e-12,
                "iauC2i00a", "13", ref status);

            vvd(rc2i[1][0], -0.2384266227870752452e-7, 1e-12,
                "iauC2i00a", "21", ref status);
            vvd(rc2i[1][1], 0.9999999991917405258, 1e-12,
                "iauC2i00a", "22", ref status);
            vvd(rc2i[1][2], -0.4020594955028209745e-4, 1e-12,
                "iauC2i00a", "23", ref status);

            vvd(rc2i[2][0], 0.5791308472168152904e-3, 1e-12,
                "iauC2i00a", "31", ref status);
            vvd(rc2i[2][1], 0.4020595661591500259e-4, 1e-12,
                "iauC2i00a", "32", ref status);
            vvd(rc2i[2][2], 0.9999998314954572304, 1e-12,
                "iauC2i00a", "33", ref status);

        }

        static void t_c2i00b(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 i 0 0 b
        **  - - - - - - - - -
        **
        **  Test iauC2i00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2i00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rc2i;


            Sofa.iauC2i00b(2400000.5, 53736.0, out rc2i);

            vvd(rc2i[0][0], 0.9999998323040954356, 1e-12,
                "iauC2i00b", "11", ref status);
            vvd(rc2i[0][1], 0.5581526349131823372e-9, 1e-12,
                "iauC2i00b", "12", ref status);
            vvd(rc2i[0][2], -0.5791301934855394005e-3, 1e-12,
                "iauC2i00b", "13", ref status);

            vvd(rc2i[1][0], -0.2384239285499175543e-7, 1e-12,
                "iauC2i00b", "21", ref status);
            vvd(rc2i[1][1], 0.9999999991917574043, 1e-12,
                "iauC2i00b", "22", ref status);
            vvd(rc2i[1][2], -0.4020552974819030066e-4, 1e-12,
                "iauC2i00b", "23", ref status);

            vvd(rc2i[2][0], 0.5791301929950208873e-3, 1e-12,
                "iauC2i00b", "31", ref status);
            vvd(rc2i[2][1], 0.4020553681373720832e-4, 1e-12,
                "iauC2i00b", "32", ref status);
            vvd(rc2i[2][2], 0.9999998314958529887, 1e-12,
                "iauC2i00b", "33", ref status);

        }

        static void t_c2i06a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 i 0 6 a
        **  - - - - - - - - -
        **
        **  Test iauC2i06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2i06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rc2i;


            Sofa.iauC2i06a(2400000.5, 53736.0, out rc2i);

            vvd(rc2i[0][0], 0.9999998323037159379, 1e-12,
                "iauC2i06a", "11", ref status);
            vvd(rc2i[0][1], 0.5581121329587613787e-9, 1e-12,
                "iauC2i06a", "12", ref status);
            vvd(rc2i[0][2], -0.5791308487740529749e-3, 1e-12,
                "iauC2i06a", "13", ref status);

            vvd(rc2i[1][0], -0.2384253169452306581e-7, 1e-12,
                "iauC2i06a", "21", ref status);
            vvd(rc2i[1][1], 0.9999999991917467827, 1e-12,
                "iauC2i06a", "22", ref status);
            vvd(rc2i[1][2], -0.4020579392895682558e-4, 1e-12,
                "iauC2i06a", "23", ref status);

            vvd(rc2i[2][0], 0.5791308482835292617e-3, 1e-12,
                "iauC2i06a", "31", ref status);
            vvd(rc2i[2][1], 0.4020580099454020310e-4, 1e-12,
                "iauC2i06a", "32", ref status);
            vvd(rc2i[2][2], 0.9999998314954628695, 1e-12,
                "iauC2i06a", "33", ref status);

        }

        static void t_c2ibpn(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 i b p n
        **  - - - - - - - - -
        **
        **  Test iauC2ibpn function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2ibpn, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rc2i;
            double[][] rbpn = { new double[3], new double[3], new double[3] };


            rbpn[0][0] = 9.999962358680738e-1;
            rbpn[0][1] = -2.516417057665452e-3;
            rbpn[0][2] = -1.093569785342370e-3;

            rbpn[1][0] = 2.516462370370876e-3;
            rbpn[1][1] = 9.999968329010883e-1;
            rbpn[1][2] = 4.006159587358310e-5;

            rbpn[2][0] = 1.093465510215479e-3;
            rbpn[2][1] = -4.281337229063151e-5;
            rbpn[2][2] = 9.999994012499173e-1;

            Sofa.iauC2ibpn(2400000.5, 50123.9999, rbpn, out rc2i);

            vvd(rc2i[0][0], 0.9999994021664089977, 1e-12,
                "iauC2ibpn", "11", ref status);
            vvd(rc2i[0][1], -0.3869195948017503664e-8, 1e-12,
                "iauC2ibpn", "12", ref status);
            vvd(rc2i[0][2], -0.1093465511383285076e-2, 1e-12,
                "iauC2ibpn", "13", ref status);

            vvd(rc2i[1][0], 0.5068413965715446111e-7, 1e-12,
                "iauC2ibpn", "21", ref status);
            vvd(rc2i[1][1], 0.9999999990835075686, 1e-12,
                "iauC2ibpn", "22", ref status);
            vvd(rc2i[1][2], 0.4281334246452708915e-4, 1e-12,
                "iauC2ibpn", "23", ref status);

            vvd(rc2i[2][0], 0.1093465510215479000e-2, 1e-12,
                "iauC2ibpn", "31", ref status);
            vvd(rc2i[2][1], -0.4281337229063151000e-4, 1e-12,
                "iauC2ibpn", "32", ref status);
            vvd(rc2i[2][2], 0.9999994012499173103, 1e-12,
                "iauC2ibpn", "33", ref status);

        }

        static void t_c2ixy(ref int status)
        /*
        **  - - - - - - - -
        **   t _ c 2 i x y
        **  - - - - - - - -
        **
        **  Test iauC2ixy function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2ixy, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y;
            double[][] rc2i;

            x = 0.5791308486706011000e-3;
            y = 0.4020579816732961219e-4;

            Sofa.iauC2ixy(2400000.5, 53736, x, y, out rc2i);

            vvd(rc2i[0][0], 0.9999998323037157138, 1e-12,
                "iauC2ixy", "11", ref status);
            vvd(rc2i[0][1], 0.5581526349032241205e-9, 1e-12,
                "iauC2ixy", "12", ref status);
            vvd(rc2i[0][2], -0.5791308491611263745e-3, 1e-12,
                "iauC2ixy", "13", ref status);

            vvd(rc2i[1][0], -0.2384257057469842953e-7, 1e-12,
                "iauC2ixy", "21", ref status);
            vvd(rc2i[1][1], 0.9999999991917468964, 1e-12,
                "iauC2ixy", "22", ref status);
            vvd(rc2i[1][2], -0.4020579110172324363e-4, 1e-12,
                "iauC2ixy", "23", ref status);

            vvd(rc2i[2][0], 0.5791308486706011000e-3, 1e-12,
                "iauC2ixy", "31", ref status);
            vvd(rc2i[2][1], 0.4020579816732961219e-4, 1e-12,
                "iauC2ixy", "32", ref status);
            vvd(rc2i[2][2], 0.9999998314954627590, 1e-12,
                "iauC2ixy", "33", ref status);

        }

        static void t_c2ixys(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 i x y s
        **  - - - - - - - - -
        **
        **  Test iauC2ixys function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2ixys, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y, s;
            double[][] rc2i;


            x = 0.5791308486706011000e-3;
            y = 0.4020579816732961219e-4;
            s = -0.1220040848472271978e-7;

            Sofa.iauC2ixys(x, y, s, out rc2i);

            vvd(rc2i[0][0], 0.9999998323037157138, 1e-12,
                "iauC2ixys", "11", ref status);
            vvd(rc2i[0][1], 0.5581984869168499149e-9, 1e-12,
                "iauC2ixys", "12", ref status);
            vvd(rc2i[0][2], -0.5791308491611282180e-3, 1e-12,
                "iauC2ixys", "13", ref status);

            vvd(rc2i[1][0], -0.2384261642670440317e-7, 1e-12,
                "iauC2ixys", "21", ref status);
            vvd(rc2i[1][1], 0.9999999991917468964, 1e-12,
                "iauC2ixys", "22", ref status);
            vvd(rc2i[1][2], -0.4020579110169668931e-4, 1e-12,
                "iauC2ixys", "23", ref status);

            vvd(rc2i[2][0], 0.5791308486706011000e-3, 1e-12,
                "iauC2ixys", "31", ref status);
            vvd(rc2i[2][1], 0.4020579816732961219e-4, 1e-12,
                "iauC2ixys", "32", ref status);
            vvd(rc2i[2][2], 0.9999998314954627590, 1e-12,
                "iauC2ixys", "33", ref status);

        }

        static void t_c2s(ref int status)
        /*
        **  - - - - - -
        **   t _ c 2 s
        **  - - - - - -
        **
        **  Test iauC2s function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2s, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta, phi;
            double[] p = new double[3];


            p[0] = 100.0;
            p[1] = -50.0;
            p[2] = 25.0;

            Sofa.iauC2s(p, out theta, out phi);

            vvd(theta, -0.4636476090008061162, 1e-14, "iauC2s", "theta", ref status);
            vvd(phi, 0.2199879773954594463, 1e-14, "iauC2s", "phi", ref status);

        }

        static void t_c2t00a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 t 0 0 a
        **  - - - - - - - - -
        **
        **  Test iauC2t00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2t00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double tta, ttb, uta, utb, xp, yp;
            double[][] rc2t;


            tta = 2400000.5;
            uta = 2400000.5;
            ttb = 53736.0;
            utb = 53736.0;
            xp = 2.55060238e-7;
            yp = 1.860359247e-6;

            Sofa.iauC2t00a(tta, ttb, uta, utb, xp, yp, out rc2t);

            vvd(rc2t[0][0], -0.1810332128307182668, 1e-12,
                "iauC2t00a", "11", ref status);
            vvd(rc2t[0][1], 0.9834769806938457836, 1e-12,
                "iauC2t00a", "12", ref status);
            vvd(rc2t[0][2], 0.6555535638688341725e-4, 1e-12,
                "iauC2t00a", "13", ref status);

            vvd(rc2t[1][0], -0.9834768134135984552, 1e-12,
                "iauC2t00a", "21", ref status);
            vvd(rc2t[1][1], -0.1810332203649520727, 1e-12,
                "iauC2t00a", "22", ref status);
            vvd(rc2t[1][2], 0.5749801116141056317e-3, 1e-12,
                "iauC2t00a", "23", ref status);

            vvd(rc2t[2][0], 0.5773474014081406921e-3, 1e-12,
                "iauC2t00a", "31", ref status);
            vvd(rc2t[2][1], 0.3961832391770163647e-4, 1e-12,
                "iauC2t00a", "32", ref status);
            vvd(rc2t[2][2], 0.9999998325501692289, 1e-12,
                "iauC2t00a", "33", ref status);

        }

        static void t_c2t00b(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 t 0 0 b
        **  - - - - - - - - -
        **
        **  Test iauC2t00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2t00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double tta, ttb, uta, utb, xp, yp;
            double[][] rc2t;


            tta = 2400000.5;
            uta = 2400000.5;
            ttb = 53736.0;
            utb = 53736.0;
            xp = 2.55060238e-7;
            yp = 1.860359247e-6;

            Sofa.iauC2t00b(tta, ttb, uta, utb, xp, yp, out rc2t);

            vvd(rc2t[0][0], -0.1810332128439678965, 1e-12,
                "iauC2t00b", "11", ref status);
            vvd(rc2t[0][1], 0.9834769806913872359, 1e-12,
                "iauC2t00b", "12", ref status);
            vvd(rc2t[0][2], 0.6555565082458415611e-4, 1e-12,
                "iauC2t00b", "13", ref status);

            vvd(rc2t[1][0], -0.9834768134115435923, 1e-12,
                "iauC2t00b", "21", ref status);
            vvd(rc2t[1][1], -0.1810332203784001946, 1e-12,
                "iauC2t00b", "22", ref status);
            vvd(rc2t[1][2], 0.5749793922030017230e-3, 1e-12,
                "iauC2t00b", "23", ref status);

            vvd(rc2t[2][0], 0.5773467471863534901e-3, 1e-12,
                "iauC2t00b", "31", ref status);
            vvd(rc2t[2][1], 0.3961790411549945020e-4, 1e-12,
                "iauC2t00b", "32", ref status);
            vvd(rc2t[2][2], 0.9999998325505635738, 1e-12,
                "iauC2t00b", "33", ref status);

        }

        static void t_c2t06a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 t 0 6 a
        **  - - - - - - - - -
        **
        **  Test iauC2t06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2t06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double tta, ttb, uta, utb, xp, yp;
            double[][] rc2t;


            tta = 2400000.5;
            uta = 2400000.5;
            ttb = 53736.0;
            utb = 53736.0;
            xp = 2.55060238e-7;
            yp = 1.860359247e-6;

            Sofa.iauC2t06a(tta, ttb, uta, utb, xp, yp, out rc2t);

            vvd(rc2t[0][0], -0.1810332128305897282, 1e-12,
                "iauC2t06a", "11", ref status);
            vvd(rc2t[0][1], 0.9834769806938592296, 1e-12,
                "iauC2t06a", "12", ref status);
            vvd(rc2t[0][2], 0.6555550962998436505e-4, 1e-12,
                "iauC2t06a", "13", ref status);

            vvd(rc2t[1][0], -0.9834768134136214897, 1e-12,
                "iauC2t06a", "21", ref status);
            vvd(rc2t[1][1], -0.1810332203649130832, 1e-12,
                "iauC2t06a", "22", ref status);
            vvd(rc2t[1][2], 0.5749800844905594110e-3, 1e-12,
                "iauC2t06a", "23", ref status);

            vvd(rc2t[2][0], 0.5773474024748545878e-3, 1e-12,
                "iauC2t06a", "31", ref status);
            vvd(rc2t[2][1], 0.3961816829632690581e-4, 1e-12,
                "iauC2t06a", "32", ref status);
            vvd(rc2t[2][2], 0.9999998325501747785, 1e-12,
                "iauC2t06a", "33", ref status);

        }

        static void t_c2tcio(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 t c i o
        **  - - - - - - - - -
        **
        **  Test iauC2tcio function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2tcio, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double era;
            double[][] rc2i = { new double[3], new double[3], new double[3] };
            double[][] rpom = { new double[3], new double[3], new double[3] };
            double[][] rc2t;


            rc2i[0][0] = 0.9999998323037164738;
            rc2i[0][1] = 0.5581526271714303683e-9;
            rc2i[0][2] = -0.5791308477073443903e-3;

            rc2i[1][0] = -0.2384266227524722273e-7;
            rc2i[1][1] = 0.9999999991917404296;
            rc2i[1][2] = -0.4020594955030704125e-4;

            rc2i[2][0] = 0.5791308472168153320e-3;
            rc2i[2][1] = 0.4020595661593994396e-4;
            rc2i[2][2] = 0.9999998314954572365;

            era = 1.75283325530307;

            rpom[0][0] = 0.9999999999999674705;
            rpom[0][1] = -0.1367174580728847031e-10;
            rpom[0][2] = 0.2550602379999972723e-6;

            rpom[1][0] = 0.1414624947957029721e-10;
            rpom[1][1] = 0.9999999999982694954;
            rpom[1][2] = -0.1860359246998866338e-5;

            rpom[2][0] = -0.2550602379741215275e-6;
            rpom[2][1] = 0.1860359247002413923e-5;
            rpom[2][2] = 0.9999999999982369658;


            Sofa.iauC2tcio(rc2i, era, rpom, out rc2t);

            vvd(rc2t[0][0], -0.1810332128307110439, 1e-12,
                "iauC2tcio", "11", ref status);
            vvd(rc2t[0][1], 0.9834769806938470149, 1e-12,
                "iauC2tcio", "12", ref status);
            vvd(rc2t[0][2], 0.6555535638685466874e-4, 1e-12,
                "iauC2tcio", "13", ref status);

            vvd(rc2t[1][0], -0.9834768134135996657, 1e-12,
                "iauC2tcio", "21", ref status);
            vvd(rc2t[1][1], -0.1810332203649448367, 1e-12,
                "iauC2tcio", "22", ref status);
            vvd(rc2t[1][2], 0.5749801116141106528e-3, 1e-12,
                "iauC2tcio", "23", ref status);

            vvd(rc2t[2][0], 0.5773474014081407076e-3, 1e-12,
                "iauC2tcio", "31", ref status);
            vvd(rc2t[2][1], 0.3961832391772658944e-4, 1e-12,
                "iauC2tcio", "32", ref status);
            vvd(rc2t[2][2], 0.9999998325501691969, 1e-12,
                "iauC2tcio", "33", ref status);

        }

        static void t_c2teqx(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c 2 t e q x
        **  - - - - - - - - -
        **
        **  Test iauC2teqx function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2teqx, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double gst;
            double[][] rbpn = { new double[3], new double[3], new double[3] };
            double[][] rpom = { new double[3], new double[3], new double[3] };
            double[][] rc2t;

            rbpn[0][0] = 0.9999989440476103608;
            rbpn[0][1] = -0.1332881761240011518e-2;
            rbpn[0][2] = -0.5790767434730085097e-3;

            rbpn[1][0] = 0.1332858254308954453e-2;
            rbpn[1][1] = 0.9999991109044505944;
            rbpn[1][2] = -0.4097782710401555759e-4;

            rbpn[2][0] = 0.5791308472168153320e-3;
            rbpn[2][1] = 0.4020595661593994396e-4;
            rbpn[2][2] = 0.9999998314954572365;

            gst = 1.754166138040730516;

            rpom[0][0] = 0.9999999999999674705;
            rpom[0][1] = -0.1367174580728847031e-10;
            rpom[0][2] = 0.2550602379999972723e-6;

            rpom[1][0] = 0.1414624947957029721e-10;
            rpom[1][1] = 0.9999999999982694954;
            rpom[1][2] = -0.1860359246998866338e-5;

            rpom[2][0] = -0.2550602379741215275e-6;
            rpom[2][1] = 0.1860359247002413923e-5;
            rpom[2][2] = 0.9999999999982369658;

            Sofa.iauC2teqx(rbpn, gst, rpom, out rc2t);

            vvd(rc2t[0][0], -0.1810332128528685730, 1e-12,
                "iauC2teqx", "11", ref status);
            vvd(rc2t[0][1], 0.9834769806897685071, 1e-12,
                "iauC2teqx", "12", ref status);
            vvd(rc2t[0][2], 0.6555535639982634449e-4, 1e-12,
                "iauC2teqx", "13", ref status);

            vvd(rc2t[1][0], -0.9834768134095211257, 1e-12,
                "iauC2teqx", "21", ref status);
            vvd(rc2t[1][1], -0.1810332203871023800, 1e-12,
                "iauC2teqx", "22", ref status);
            vvd(rc2t[1][2], 0.5749801116126438962e-3, 1e-12,
                "iauC2teqx", "23", ref status);

            vvd(rc2t[2][0], 0.5773474014081539467e-3, 1e-12,
                "iauC2teqx", "31", ref status);
            vvd(rc2t[2][1], 0.3961832391768640871e-4, 1e-12,
                "iauC2teqx", "32", ref status);
            vvd(rc2t[2][2], 0.9999998325501691969, 1e-12,
                "iauC2teqx", "33", ref status);

        }

        static void t_c2tpe(ref int status)
        /*
        **  - - - - - - - -
        **   t _ c 2 t p e
        **  - - - - - - - -
        **
        **  Test iauC2tpe function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2tpe, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double tta, ttb, uta, utb, dpsi, deps, xp, yp;
            double[][] rc2t;


            tta = 2400000.5;
            uta = 2400000.5;
            ttb = 53736.0;
            utb = 53736.0;
            deps = 0.4090789763356509900;
            dpsi = -0.9630909107115582393e-5;
            xp = 2.55060238e-7;
            yp = 1.860359247e-6;

            Sofa.iauC2tpe(tta, ttb, uta, utb, dpsi, deps, xp, yp, out rc2t);

            vvd(rc2t[0][0], -0.1813677995763029394, 1e-12,
                "iauC2tpe", "11", ref status);
            vvd(rc2t[0][1], 0.9023482206891683275, 1e-12,
                "iauC2tpe", "12", ref status);
            vvd(rc2t[0][2], -0.3909902938641085751, 1e-12,
                "iauC2tpe", "13", ref status);

            vvd(rc2t[1][0], -0.9834147641476804807, 1e-12,
                "iauC2tpe", "21", ref status);
            vvd(rc2t[1][1], -0.1659883635434995121, 1e-12,
                "iauC2tpe", "22", ref status);
            vvd(rc2t[1][2], 0.7309763898042819705e-1, 1e-12,
                "iauC2tpe", "23", ref status);

            vvd(rc2t[2][0], 0.1059685430673215247e-2, 1e-12,
                "iauC2tpe", "31", ref status);
            vvd(rc2t[2][1], 0.3977631855605078674, 1e-12,
                "iauC2tpe", "32", ref status);
            vvd(rc2t[2][2], 0.9174875068792735362, 1e-12,
                "iauC2tpe", "33", ref status);

        }

        static void t_c2txy(ref int status)
        /*
        **  - - - - - - - -
        **   t _ c 2 t x y
        **  - - - - - - - -
        **
        **  Test iauC2txy function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauC2txy, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double tta, ttb, uta, utb, x, y, xp, yp;
            double[][] rc2t;


            tta = 2400000.5;
            uta = 2400000.5;
            ttb = 53736.0;
            utb = 53736.0;
            x = 0.5791308486706011000e-3;
            y = 0.4020579816732961219e-4;
            xp = 2.55060238e-7;
            yp = 1.860359247e-6;

            Sofa.iauC2txy(tta, ttb, uta, utb, x, y, xp, yp, out rc2t);

            vvd(rc2t[0][0], -0.1810332128306279253, 1e-12,
                "iauC2txy", "11", ref status);
            vvd(rc2t[0][1], 0.9834769806938520084, 1e-12,
                "iauC2txy", "12", ref status);
            vvd(rc2t[0][2], 0.6555551248057665829e-4, 1e-12,
                "iauC2txy", "13", ref status);

            vvd(rc2t[1][0], -0.9834768134136142314, 1e-12,
                "iauC2txy", "21", ref status);
            vvd(rc2t[1][1], -0.1810332203649529312, 1e-12,
                "iauC2txy", "22", ref status);
            vvd(rc2t[1][2], 0.5749800843594139912e-3, 1e-12,
                "iauC2txy", "23", ref status);

            vvd(rc2t[2][0], 0.5773474028619264494e-3, 1e-12,
                "iauC2txy", "31", ref status);
            vvd(rc2t[2][1], 0.3961816546911624260e-4, 1e-12,
                "iauC2txy", "32", ref status);
            vvd(rc2t[2][2], 0.9999998325501746670, 1e-12,
                "iauC2txy", "33", ref status);

        }

        static void t_cal2jd(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ c a l 2 j d
        **  - - - - - - - - -
        **
        **  Test iauCal2jd function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauCal2jd, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            int j;
            double djm0, djm;


            j = Sofa.iauCal2jd(2003, 06, 01, out djm0, out djm);

            vvd(djm0, 2400000.5, 0.0, "iauCal2jd", "djm0", ref status);
            vvd(djm, 52791.0, 0.0, "iauCal2jd", "djm", ref status);

            viv(j, 0, "iauCal2jd", "j", ref status);

        }

        static void t_cp(ref int status)
        /*
        **  - - - - -
        **   t _ c p
        **  - - - - -
        **
        **  Test iauCp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauCp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p = new double[3];
            double[] c = new double[3];


            p[0] = 0.3;
            p[1] = 1.2;
            p[2] = -2.5;

            Sofa.iauCp(p, c);

            vvd(c[0], 0.3, 0.0, "iauCp", "1", ref status);
            vvd(c[1], 1.2, 0.0, "iauCp", "2", ref status);
            vvd(c[2], -2.5, 0.0, "iauCp", "3", ref status);
        }

        static void t_cpv(ref int status)
        /*
        **  - - - - - -
        **   t _ c p v
        **  - - - - - -
        **
        **  Test iauCpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauCpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] pv = { new double[3], new double[3] };
            double[][] c = { new double[3], new double[3] };


            pv[0][0] = 0.3;
            pv[0][1] = 1.2;
            pv[0][2] = -2.5;

            pv[1][0] = -0.5;
            pv[1][1] = 3.1;
            pv[1][2] = 0.9;

            Sofa.iauCpv(pv, c);

            vvd(c[0][0], 0.3, 0.0, "iauCpv", "p1", ref status);
            vvd(c[0][1], 1.2, 0.0, "iauCpv", "p2", ref status);
            vvd(c[0][2], -2.5, 0.0, "iauCpv", "p3", ref status);

            vvd(c[1][0], -0.5, 0.0, "iauCpv", "v1", ref status);
            vvd(c[1][1], 3.1, 0.0, "iauCpv", "v2", ref status);
            vvd(c[1][2], 0.9, 0.0, "iauCpv", "v3", ref status);

        }

        static void t_cr(ref int status)
        /*
        **  - - - - -
        **   t _ c r
        **  - - - - -
        **
        **  Test iauCr function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauCr, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r = { new double[3], new double[3], new double[3] };
            double[][] c = { new double[3], new double[3], new double[3] };


            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            Sofa.iauCr(r, c);

            vvd(c[0][0], 2.0, 0.0, "iauCr", "11", ref status);
            vvd(c[0][1], 3.0, 0.0, "iauCr", "12", ref status);
            vvd(c[0][2], 2.0, 0.0, "iauCr", "13", ref status);

            vvd(c[1][0], 3.0, 0.0, "iauCr", "21", ref status);
            vvd(c[1][1], 2.0, 0.0, "iauCr", "22", ref status);
            vvd(c[1][2], 3.0, 0.0, "iauCr", "23", ref status);

            vvd(c[2][0], 3.0, 0.0, "iauCr", "31", ref status);
            vvd(c[2][1], 4.0, 0.0, "iauCr", "32", ref status);
            vvd(c[2][2], 5.0, 0.0, "iauCr", "33", ref status);
        }

        static void t_d2dtf(ref int status)
        /*
        **  - - - - - - - -
        **   t _ d 2 d t f
        **  - - - - - - - -
        **
        **  Test iauD2dtf function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauD2dtf, viv
        **
        **  This revision:  2013 August 7
        */
        {
            int j, iy, im, id;
            int[] ihmsf;


            j = Sofa.iauD2dtf("UTC", 5, 2400000.5, 49533.99999, out iy, out im, out id, out ihmsf);

            viv(iy, 1994, "iauD2dtf", "y", ref status);
            viv(im, 6, "iauD2dtf", "mo", ref status);
            viv(id, 30, "iauD2dtf", "d", ref status);
            viv(ihmsf[0], 23, "iauD2dtf", "h", ref status);
            viv(ihmsf[1], 59, "iauD2dtf", "m", ref status);
            viv(ihmsf[2], 60, "iauD2dtf", "s", ref status);
            viv(ihmsf[3], 13599, "iauD2dtf", "f", ref status);
            viv(j, 0, "iauD2dtf", "j", ref status);

        }

        static void t_d2tf(ref int status)
        /*
        **  - - - - - - -
        **   t _ d 2 t f
        **  - - - - - - -
        **
        **  Test iauD2tf function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauD2tf, viv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            int[] ihmsf = new int[4];
            char s;


            Sofa.iauD2tf(4, -0.987654321, out s, out ihmsf);

            viv((int)s, '-', "iauD2tf", "s", ref status);

            viv(ihmsf[0], 23, "iauD2tf", "0", ref status);
            viv(ihmsf[1], 42, "iauD2tf", "1", ref status);
            viv(ihmsf[2], 13, "iauD2tf", "2", ref status);
            viv(ihmsf[3], 3333, "iauD2tf", "3", ref status);

        }

        static void t_dat(ref int status)
        /*
        **  - - - - - -
        **   t _ d a t
        **  - - - - - -
        **
        **  Test iauDat function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauDat, vvd, viv
        **
        **  This revision:  2016 July 11
        */
        {
            int j;
            double deltat;


            j = Sofa.iauDat(2003, 6, 1, 0.0, out deltat);

            vvd(deltat, 32.0, 0.0, "iauDat", "d1", ref status);
            viv(j, 0, "iauDat", "j1", ref status);

            j = Sofa.iauDat(2008, 1, 17, 0.0, out deltat);

            vvd(deltat, 33.0, 0.0, "iauDat", "d2", ref status);
            viv(j, 0, "iauDat", "j2", ref status);

            j = Sofa.iauDat(2017, 9, 1, 0.0, out deltat);

            vvd(deltat, 37.0, 0.0, "iauDat", "d3", ref status);
            viv(j, 0, "iauDat", "j3", ref status);

        }

        static void t_dtdb(ref int status)
        /*
        **  - - - - - - -
        **   t _ d t d b
        **  - - - - - - -
        **
        **  Test iauDtdb function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauDtdb, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dtdb;


            dtdb = Sofa.iauDtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0);

            vvd(dtdb, -0.1280368005936998991e-2, 1e-15, "iauDtdb", "", ref status);

        }

        static void t_dtf2d(ref int status)
        /*
        **  - - - - - - - -
        **   t _ d t f 2 d
        **  - - - - - - - -
        **
        **  Test iauDtf2d function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauDtf2d, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double u1, u2;
            int j;


            j = Sofa.iauDtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599, out u1, out u2);

            vvd(u1 + u2, 2449534.49999, 1e-6, "iauDtf2d", "u", ref status);
            viv(j, 0, "iauDtf2d", "j", ref status);

        }

        static void t_eceq06(ref int status)
        /*
        **  - - - - -
        **   t _ e c e q 0 6
        **  - - - - -
        **
        **  Test iauEceq06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEceq06, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double date1, date2, dl, db, dr, dd;


            date1 = 2456165.5;
            date2 = 0.401182685;
            dl = 5.1;
            db = -0.9;

            Sofa.iauEceq06(date1, date2, dl, db, out dr, out dd);

            vvd(dr, 5.533459733613627767, 1e-14, "iauEceq06", "dr", ref status);
            vvd(dd, -1.246542932554480576, 1e-14, "iauEceq06", "dd", ref status);

        }

        static void t_ecm06(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e c m 0 6
        **  - - - - - - - -
        **
        **  Test iauEcm06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEcm06, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double date1, date2;
            double[][] rm;


            date1 = 2456165.5;
            date2 = 0.401182685;

            Sofa.iauEcm06(date1, date2, out rm);

            vvd(rm[0][0], 0.9999952427708701137, 1e-14,
                "iauEcm06", "rm11", ref status);
            vvd(rm[0][1], -0.2829062057663042347e-2, 1e-14,
                "iauEcm06", "rm12", ref status);
            vvd(rm[0][2], -0.1229163741100017629e-2, 1e-14,
                "iauEcm06", "rm13", ref status);
            vvd(rm[1][0], 0.3084546876908653562e-2, 1e-14,
                "iauEcm06", "rm21", ref status);
            vvd(rm[1][1], 0.9174891871550392514, 1e-14,
                "iauEcm06", "rm22", ref status);
            vvd(rm[1][2], 0.3977487611849338124, 1e-14,
                "iauEcm06", "rm23", ref status);
            vvd(rm[2][0], 0.2488512951527405928e-5, 1e-14,
                "iauEcm06", "rm31", ref status);
            vvd(rm[2][1], -0.3977506604161195467, 1e-14,
                "iauEcm06", "rm32", ref status);
            vvd(rm[2][2], 0.9174935488232863071, 1e-14,
                "iauEcm06", "rm33", ref status);

        }

        static void t_ee00(ref int status)
        /*
        **  - - - - - - -
        **   t _ e e 0 0
        **  - - - - - - -
        **
        **  Test iauEe00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEe00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double epsa, dpsi, ee;


            epsa = 0.4090789763356509900;
            dpsi = -0.9630909107115582393e-5;

            ee = Sofa.iauEe00(2400000.5, 53736.0, epsa, dpsi);

            vvd(ee, -0.8834193235367965479e-5, 1e-18, "iauEe00", "", ref status);

        }

        static void t_ee00a(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e e 0 0 a
        **  - - - - - - - -
        **
        **  Test iauEe00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEe00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double ee;


            ee = Sofa.iauEe00a(2400000.5, 53736.0);

            vvd(ee, -0.8834192459222588227e-5, 1e-18, "iauEe00a", "", ref status);

        }

        static void t_ee00b(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e e 0 0 b
        **  - - - - - - - -
        **
        **  Test iauEe00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEe00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double ee;


            ee = Sofa.iauEe00b(2400000.5, 53736.0);

            vvd(ee, -0.8835700060003032831e-5, 1e-18, "iauEe00b", "", ref status);

        }

        static void t_ee06a(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e e 0 6 a
        **  - - - - - - - -
        **
        **  Test iauEe06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEe06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double ee;


            ee = Sofa.iauEe06a(2400000.5, 53736.0);

            vvd(ee, -0.8834195072043790156e-5, 1e-15, "iauEe06a", "", ref status);
        }

        static void t_eect00(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ e e c t 0 0
        **  - - - - - - - - -
        **
        **  Test iauEect00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEect00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double eect;


            eect = Sofa.iauEect00(2400000.5, 53736.0);

            vvd(eect, 0.2046085004885125264e-8, 1e-20, "iauEect00", "", ref status);

        }

        static void t_eform(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e f o r m
        **  - - - - - - - -
        **
        **  Test iauEform function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEform, viv, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            int j;
            double a, f;

            j = Sofa.iauEform(0, out a, out f);

            viv(j, -1, "iauEform", "j0", ref status);

            j = Sofa.iauEform(Sofa.ReferenceEllipsoids.WGS84, out a, out f);

            viv(j, 0, "iauEform", "j1", ref status);
            vvd(a, 6378137.0, 1e-10, "iauEform", "a1", ref status);
            vvd(f, 0.3352810664747480720e-2, 1e-18, "iauEform", "f1", ref status);

            j = Sofa.iauEform(Sofa.ReferenceEllipsoids.GRS80, out a, out f);

            viv(j, 0, "iauEform", "j2", ref status);
            vvd(a, 6378137.0, 1e-10, "iauEform", "a2", ref status);
            vvd(f, 0.3352810681182318935e-2, 1e-18, "iauEform", "f2", ref status);

            j = Sofa.iauEform(Sofa.ReferenceEllipsoids.WGS72, out a, out f);

            viv(j, 0, "iauEform", "j2", ref status);
            vvd(a, 6378135.0, 1e-10, "iauEform", "a3", ref status);
            vvd(f, 0.3352779454167504862e-2, 1e-18, "iauEform", "f3", ref status);

            //j = Sofa.iauEform(4, out a, out f);
            //viv(j, -1, "iauEform", "j3", ref status);
        }

        static void t_eo06a(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e o 0 6 a
        **  - - - - - - - -
        **
        **  Test iauEo06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEo06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double eo;


            eo = Sofa.iauEo06a(2400000.5, 53736.0);

            vvd(eo, -0.1332882371941833644e-2, 1e-15, "iauEo06a", "", ref status);

        }

        static void t_eors(ref int status)
        /*
        **  - - - - - - -
        **   t _ e o r s
        **  - - - - - - -
        **
        **  Test iauEors function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEors, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rnpb = { new double[3], new double[3], new double[3] };
            double s, eo;


            rnpb[0][0] = 0.9999989440476103608;
            rnpb[0][1] = -0.1332881761240011518e-2;
            rnpb[0][2] = -0.5790767434730085097e-3;

            rnpb[1][0] = 0.1332858254308954453e-2;
            rnpb[1][1] = 0.9999991109044505944;
            rnpb[1][2] = -0.4097782710401555759e-4;

            rnpb[2][0] = 0.5791308472168153320e-3;
            rnpb[2][1] = 0.4020595661593994396e-4;
            rnpb[2][2] = 0.9999998314954572365;

            s = -0.1220040848472271978e-7;

            eo = Sofa.iauEors(rnpb, s);

            vvd(eo, -0.1332882715130744606e-2, 1e-14, "iauEors", "", ref status);

        }

        static void t_epb(ref int status)
        /*
        **  - - - - - -
        **   t _ e p b
        **  - - - - - -
        **
        **  Test iauEpb function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEpb, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double epb;


            epb = Sofa.iauEpb(2415019.8135, 30103.18648);

            vvd(epb, 1982.418424159278580, 1e-12, "iauEpb", "", ref status);

        }

        static void t_epb2jd(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ e p b 2 j d
        **  - - - - - - - - -
        **
        **  Test iauEpb2jd function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEpb2jd, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double epb, djm0, djm;


            epb = 1957.3;

            Sofa.iauEpb2jd(epb, out djm0, out djm);

            vvd(djm0, 2400000.5, 1e-9, "iauEpb2jd", "djm0", ref status);
            vvd(djm, 35948.1915101513, 1e-9, "iauEpb2jd", "mjd", ref status);

        }

        static void t_epj(ref int status)
        /*
        **  - - - - - -
        **   t _ e p j
        **  - - - - - -
        **
        **  Test iauEpj function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEpj, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double epj;


            epj = Sofa.iauEpj(2451545, -7392.5);

            vvd(epj, 1979.760438056125941, 1e-12, "iauEpj", "", ref status);

        }

        static void t_epj2jd(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ e p j 2 j d
        **  - - - - - - - - -
        **
        **  Test iauEpj2jd function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEpj2jd, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double epj, djm0, djm;


            epj = 1996.8;

            Sofa.iauEpj2jd(epj, out djm0, out djm);

            vvd(djm0, 2400000.5, 1e-9, "iauEpj2jd", "djm0", ref status);
            vvd(djm, 50375.7, 1e-9, "iauEpj2jd", "mjd", ref status);

        }

        static void t_epv00(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e p v 0 0
        **  - - - - - - - -
        **
        **  Test iauEpv00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called: iauEpv00, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] pvh, pvb;
            int j;


            j = Sofa.iauEpv00(2400000.5, 53411.52501161, out pvh, out pvb);

            vvd(pvh[0][0], -0.7757238809297706813, 1e-14,
                "iauEpv00", "ph(x)", ref status);
            vvd(pvh[0][1], 0.5598052241363340596, 1e-14,
                "iauEpv00", "ph(y)", ref status);
            vvd(pvh[0][2], 0.2426998466481686993, 1e-14,
                "iauEpv00", "ph(z)", ref status);

            vvd(pvh[1][0], -0.1091891824147313846e-1, 1e-15,
                "iauEpv00", "vh(x)", ref status);
            vvd(pvh[1][1], -0.1247187268440845008e-1, 1e-15,
                "iauEpv00", "vh(y)", ref status);
            vvd(pvh[1][2], -0.5407569418065039061e-2, 1e-15,
                "iauEpv00", "vh(z)", ref status);

            vvd(pvb[0][0], -0.7714104440491111971, 1e-14,
                "iauEpv00", "pb(x)", ref status);
            vvd(pvb[0][1], 0.5598412061824171323, 1e-14,
                "iauEpv00", "pb(y)", ref status);
            vvd(pvb[0][2], 0.2425996277722452400, 1e-14,
                "iauEpv00", "pb(z)", ref status);

            vvd(pvb[1][0], -0.1091874268116823295e-1, 1e-15,
                "iauEpv00", "vb(x)", ref status);
            vvd(pvb[1][1], -0.1246525461732861538e-1, 1e-15,
                "iauEpv00", "vb(y)", ref status);
            vvd(pvb[1][2], -0.5404773180966231279e-2, 1e-15,
                "iauEpv00", "vb(z)", ref status);

            viv(j, 0, "iauEpv00", "j", ref status);

        }

        static void t_eqec06(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ e q e c 0 6
        **  - - - - - - - - -
        **
        **  Test iauEqec06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEqec06, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double date1, date2, dr, dd, dl, db;


            date1 = 1234.5;
            date2 = 2440000.5;
            dr = 1.234;
            dd = 0.987;

            Sofa.iauEqec06(date1, date2, dr, dd, out dl, out db);

            vvd(dl, 1.342509918994654619, 1e-14, "iauEqec06", "dl", ref status);
            vvd(db, 0.5926215259704608132, 1e-14, "iauEqec06", "db", ref status);

        }

        static void t_eqeq94(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ e q e q 9 4
        **  - - - - - - - - -
        **
        **  Test iauEqeq94 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEqeq94, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double eqeq;


            eqeq = Sofa.iauEqeq94(2400000.5, 41234.0);

            vvd(eqeq, 0.5357758254609256894e-4, 1e-17, "iauEqeq94", "", ref status);

        }

        static void t_era00(ref int status)
        /*
        **  - - - - - - - -
        **   t _ e r a 0 0
        **  - - - - - - - -
        **
        **  Test iauEra00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauEra00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double era00;


            era00 = Sofa.iauEra00(2400000.5, 54388.0);

            vvd(era00, 0.4022837240028158102, 1e-12, "iauEra00", "", ref status);

        }

        static void t_fad03(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f a d 0 3
        **  - - - - - - - -
        **
        **  Test iauFad03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFad03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFad03(0.80), 1.946709205396925672, 1e-12,
                "iauFad03", "", ref status);
        }

        static void t_fae03(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f a e 0 3
        **  - - - - - - - -
        **
        **  Test iauFae03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFae03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFae03(0.80), 1.744713738913081846, 1e-12,
                "iauFae03", "", ref status);
        }

        static void t_faf03(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f a f 0 3
        **  - - - - - - - -
        **
        **  Test iauFaf03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFaf03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFaf03(0.80), 0.2597711366745499518, 1e-12,
                "iauFaf03", "", ref status);
        }

        static void t_faju03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a j u 0 3
        **  - - - - - - - - -
        **
        **  Test iauFaju03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFaju03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFaju03(0.80), 5.275711665202481138, 1e-12,
                "iauFaju03", "", ref status);
        }

        static void t_fal03(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f a l 0 3
        **  - - - - - - - -
        **
        **  Test iauFal03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFal03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFal03(0.80), 5.132369751108684150, 1e-12,
                "iauFal03", "", ref status);
        }

        static void t_falp03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a l p 0 3
        **  - - - - - - - - -
        **
        **  Test iauFalp03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFalp03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFalp03(0.80), 6.226797973505507345, 1e-12,
               "iauFalp03", "", ref status);
        }

        static void t_fama03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a m a 0 3
        **  - - - - - - - - -
        **
        **  Test iauFama03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFama03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFama03(0.80), 3.275506840277781492, 1e-12,
                "iauFama03", "", ref status);
        }

        static void t_fame03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a m e 0 3
        **  - - - - - - - - -
        **
        **  Test iauFame03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFame03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFame03(0.80), 5.417338184297289661, 1e-12,
                "iauFame03", "", ref status);
        }

        static void t_fane03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a n e 0 3
        **  - - - - - - - - -
        **
        **  Test iauFane03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFane03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFane03(0.80), 2.079343830860413523, 1e-12,
                "iauFane03", "", ref status);
        }

        static void t_faom03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a o m 0 3
        **  - - - - - - - - -
        **
        **  Test iauFaom03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFaom03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFaom03(0.80), -5.973618440951302183, 1e-12,
                "iauFaom03", "", ref status);
        }

        static void t_fapa03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a p a 0 3
        **  - - - - - - - - -
        **
        **  Test iauFapa03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFapa03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFapa03(0.80), 0.1950884762240000000e-1, 1e-12,
                "iauFapa03", "", ref status);
        }

        static void t_fasa03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a s a 0 3
        **  - - - - - - - - -
        **
        **  Test iauFasa03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFasa03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFasa03(0.80), 5.371574539440827046, 1e-12,
                "iauFasa03", "", ref status);
        }

        static void t_faur03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a u r 0 3
        **  - - - - - - - - -
        **
        **  Test iauFaur03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFaur03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFaur03(0.80), 5.180636450180413523, 1e-12,
                "iauFaur03", "", ref status);
        }

        static void t_fave03(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f a v e 0 3
        **  - - - - - - - - -
        **
        **  Test iauFave03 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFave03, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauFave03(0.80), 3.424900460533758000, 1e-12,
                "iauFave03", "", ref status);
        }

        static void t_fk425(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f k 4 2 5
        **  - - - - - - - -
        **
        **  Test iauFk425 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFk425, vvd
        **
        **  This revision:  2018 December 6
        */
        {
            double r1950, d1950, dr1950, dd1950, p1950, v1950,
                   r2000, d2000, dr2000, dd2000, p2000, v2000;


            r1950 = 0.07626899753879587532;
            d1950 = -1.137405378399605780;
            dr1950 = 0.1973749217849087460e-4;
            dd1950 = 0.5659714913272723189e-5;
            p1950 = 0.134;
            v1950 = 8.7;

            Sofa.iauFk425(r1950, d1950, dr1950, dd1950, p1950, v1950,
                     out r2000, out d2000, out dr2000, out dd2000, out p2000, out v2000);

            vvd(r2000, 0.08757989933556446040, 1e-14,
                "iauFk425", "r2000", ref status);
            vvd(d2000, -1.132279113042091895, 1e-12,
                "iauFk425", "d2000", ref status);
            vvd(dr2000, 0.1953670614474396139e-4, 1e-17,
                "iauFk425", "dr2000", ref status);
            vvd(dd2000, 0.5637686678659640164e-5, 1e-18,
                "iauFk425", "dd2000", ref status);
            vvd(p2000, 0.1339919950582767871, 1e-13, "iauFk425", "p2000", ref status);
            vvd(v2000, 8.736999669183529069, 1e-12, "iauFk425", "v2000", ref status);

        }

        static void t_fk45z(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f k 4 5 z
        **  - - - - - - - -
        **
        **  Test iauFk45z function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFk45z, vvd
        **
        **  This revision:  2018 December 6
        */
        {
            double r1950, d1950, bepoch, r2000, d2000;


            r1950 = 0.01602284975382960982;
            d1950 = -0.1164347929099906024;
            bepoch = 1954.677617625256806;

            Sofa.iauFk45z(r1950, d1950, bepoch, out r2000, out d2000);

            vvd(r2000, 0.02719295911606862303, 1e-15,
                "iauFk45z", "r2000", ref status);
            vvd(d2000, -0.1115766001565926892, 1e-13,
                "iauFk45z", "d2000", ref status);

        }

        static void t_fk524(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f k 5 2 4
        **  - - - - - - - -
        **
        **  Test iauFk524 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFk524, vvd
        **
        **  This revision:  2018 December 6
        */
        {
            double r2000, d2000, dr2000, dd2000, p2000, v2000,
                   r1950, d1950, dr1950, dd1950, p1950, v1950;


            r2000 = 0.8723503576487275595;
            d2000 = -0.7517076365138887672;
            dr2000 = 0.2019447755430472323e-4;
            dd2000 = 0.3541563940505160433e-5;
            p2000 = 0.1559;
            v2000 = 86.87;

            Sofa.iauFk524(r2000, d2000, dr2000, dd2000, p2000, v2000,
                     out r1950, out d1950, out dr1950, out dd1950, out p1950, out v1950);

            vvd(r1950, 0.8636359659799603487, 1e-13,
                "iauFk524", "r1950", ref status);
            vvd(d1950, -0.7550281733160843059, 1e-13,
                "iauFk524", "d1950", ref status);
            vvd(dr1950, 0.2023628192747172486e-4, 1e-17,
                "iauFk524", "dr1950", ref status);
            vvd(dd1950, 0.3624459754935334718e-5, 1e-18,
                "iauFk524", "dd1950", ref status);
            vvd(p1950, 0.1560079963299390241, 1e-13,
                "iauFk524", "p1950", ref status);
            vvd(v1950, 86.79606353469163751, 1e-11, "iauFk524", "v1950", ref status);

        }

        static void t_fk52h(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f k 5 2 h
        **  - - - - - - - -
        **
        **  Test iauFk52h function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFk52h, vvd
        **
        **  This revision:  2021 January 5
        */
        {
            double r5, d5, dr5, dd5, px5, rv5, rh, dh, drh, ddh, pxh, rvh;


            r5 = 1.76779433;
            d5 = -0.2917517103;
            dr5 = -1.91851572e-7;
            dd5 = -5.8468475e-6;
            px5 = 0.379210;
            rv5 = -7.6;

            Sofa.iauFk52h(r5, d5, dr5, dd5, px5, rv5,
                     out rh, out dh, out drh, out ddh, out pxh, out rvh);

            vvd(rh, 1.767794226299947632, 1e-14,
                "iauFk52h", "ra", ref status);
            vvd(dh, -0.2917516070530391757, 1e-14,
                "iauFk52h", "dec", ref status);
            vvd(drh, -0.1961874125605721270e-6, 1e-19,
                "iauFk52h", "dr5", ref status);
            vvd(ddh, -0.58459905176693911e-5, 1e-19,
                "iauFk52h", "dd5", ref status);
            vvd(pxh, 0.37921, 1e-14,
                "iauFk52h", "px", ref status);
            vvd(rvh, -7.6000000940000254, 1e-11,
                "iauFk52h", "rv", ref status);

        }

        static void t_fk54z(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f k 5 4 z
        **  - - - - - - - -
        **
        **  Test iauFk54z function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFk54z, vvd
        **
        **  This revision:  2018 December 6
        */
        {
            double r2000, d2000, bepoch, r1950, d1950, dr1950, dd1950;


            r2000 = 0.02719026625066316119;
            d2000 = -0.1115815170738754813;
            bepoch = 1954.677308160316374;

            Sofa.iauFk54z(r2000, d2000, bepoch, out r1950, out d1950, out dr1950, out dd1950);

            vvd(r1950, 0.01602015588390065476, 1e-14,
                "iauFk54z", "r1950", ref status);
            vvd(d1950, -0.1164397101110765346, 1e-13,
                "iauFk54z", "d1950", ref status);
            vvd(dr1950, -0.1175712648471090704e-7, 1e-20,
                "iauFk54z", "dr1950", ref status);
            vvd(dd1950, 0.2108109051316431056e-7, 1e-20,
                "iauFk54z", "dd1950", ref status);

        }

        static void t_fk5hip(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ f k 5 h i p
        **  - - - - - - - - -
        **
        **  Test iauFk5hip function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFk5hip, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r5h;
            double[] s5h;


            Sofa.iauFk5hip(out r5h, out s5h);

            vvd(r5h[0][0], 0.9999999999999928638, 1e-14,
                "iauFk5hip", "11", ref status);
            vvd(r5h[0][1], 0.1110223351022919694e-6, 1e-17,
                "iauFk5hip", "12", ref status);
            vvd(r5h[0][2], 0.4411803962536558154e-7, 1e-17,
                "iauFk5hip", "13", ref status);
            vvd(r5h[1][0], -0.1110223308458746430e-6, 1e-17,
                "iauFk5hip", "21", ref status);
            vvd(r5h[1][1], 0.9999999999999891830, 1e-14,
                "iauFk5hip", "22", ref status);
            vvd(r5h[1][2], -0.9647792498984142358e-7, 1e-17,
                "iauFk5hip", "23", ref status);
            vvd(r5h[2][0], -0.4411805033656962252e-7, 1e-17,
                "iauFk5hip", "31", ref status);
            vvd(r5h[2][1], 0.9647792009175314354e-7, 1e-17,
                "iauFk5hip", "32", ref status);
            vvd(r5h[2][2], 0.9999999999999943728, 1e-14,
                "iauFk5hip", "33", ref status);
            vvd(s5h[0], -0.1454441043328607981e-8, 1e-17,
                "iauFk5hip", "s1", ref status);
            vvd(s5h[1], 0.2908882086657215962e-8, 1e-17,
                "iauFk5hip", "s2", ref status);
            vvd(s5h[2], 0.3393695767766751955e-8, 1e-17,
                "iauFk5hip", "s3", ref status);

        }

        static void t_fk5hz(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f k 5 h z
        **  - - - - - - - -
        **
        **  Test iauFk5hz function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFk5hz, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double r5, d5, rh, dh;


            r5 = 1.76779433;
            d5 = -0.2917517103;

            Sofa.iauFk5hz(r5, d5, 2400000.5, 54479.0, out rh, out dh);

            vvd(rh, 1.767794191464423978, 1e-12, "iauFk5hz", "ra", ref status);
            vvd(dh, -0.2917516001679884419, 1e-12, "iauFk5hz", "dec", ref status);

        }

        static void t_fw2m(ref int status)
        /*
        **  - - - - - - -
        **   t _ f w 2 m
        **  - - - - - - -
        **
        **  Test iauFw2m function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFw2m, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double gamb, phib, psi, eps;
            double[][] r;


            gamb = -0.2243387670997992368e-5;
            phib = 0.4091014602391312982;
            psi = -0.9501954178013015092e-3;
            eps = 0.4091014316587367472;

            Sofa.iauFw2m(gamb, phib, psi, eps, out r);

            vvd(r[0][0], 0.9999995505176007047, 1e-12,
                "iauFw2m", "11", ref status);
            vvd(r[0][1], 0.8695404617348192957e-3, 1e-12,
                "iauFw2m", "12", ref status);
            vvd(r[0][2], 0.3779735201865582571e-3, 1e-12,
                "iauFw2m", "13", ref status);

            vvd(r[1][0], -0.8695404723772016038e-3, 1e-12,
                "iauFw2m", "21", ref status);
            vvd(r[1][1], 0.9999996219496027161, 1e-12,
                "iauFw2m", "22", ref status);
            vvd(r[1][2], -0.1361752496887100026e-6, 1e-12,
                "iauFw2m", "23", ref status);

            vvd(r[2][0], -0.3779734957034082790e-3, 1e-12,
                "iauFw2m", "31", ref status);
            vvd(r[2][1], -0.1924880848087615651e-6, 1e-12,
                "iauFw2m", "32", ref status);
            vvd(r[2][2], 0.9999999285679971958, 1e-12,
                "iauFw2m", "33", ref status);

        }

        static void t_fw2xy(ref int status)
        /*
        **  - - - - - - - -
        **   t _ f w 2 x y
        **  - - - - - - - -
        **
        **  Test iauFw2xy function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauFw2xy, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double gamb, phib, psi, eps, x, y;


            gamb = -0.2243387670997992368e-5;
            phib = 0.4091014602391312982;
            psi = -0.9501954178013015092e-3;
            eps = 0.4091014316587367472;

            Sofa.iauFw2xy(gamb, phib, psi, eps, out x, out y);

            vvd(x, -0.3779734957034082790e-3, 1e-14, "iauFw2xy", "x", ref status);
            vvd(y, -0.1924880848087615651e-6, 1e-14, "iauFw2xy", "y", ref status);

        }

        static void t_g2icrs(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g 2 i c r s
        **  - - - - - - - - -
        **
        **  Test iauG2icrs function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauG2icrs, vvd
        **
        **  This revision:  2015 January 30
        */
        {
            double dl, db, dr, dd;


            dl = 5.5850536063818546461558105;
            db = -0.7853981633974483096156608;
            Sofa.iauG2icrs(dl, db, out dr, out dd);
            vvd(dr, 5.9338074302227188048671, 1e-14, "iauG2icrs", "R", ref status);
            vvd(dd, -1.1784870613579944551541, 1e-14, "iauG2icrs", "D", ref status);
        }

        static void t_gc2gd(ref int status)
        /*
        **  - - - - - - - -
        **   t _ g c 2 g d
        **  - - - - - - - -
        **
        **  Test iauGc2gd function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGc2gd, viv, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            int j;
            double[] xyz = { 2e6, 3e6, 5.244e6 };
            double e, p, h;

            j = Sofa.iauGc2gd(0, xyz, out e, out p, out h);

            viv(j, -1, "iauGc2gd", "j0", ref status);

            j = Sofa.iauGc2gd(Sofa.ReferenceEllipsoids.WGS84, xyz, out e, out p, out h);

            viv(j, 0, "iauGc2gd", "j1", ref status);
            vvd(e, 0.9827937232473290680, 1e-14, "iauGc2gd", "e1", ref status);
            vvd(p, 0.97160184819075459, 1e-14, "iauGc2gd", "p1", ref status);
            vvd(h, 331.4172461426059892, 1e-8, "iauGc2gd", "h1", ref status);

            j = Sofa.iauGc2gd(Sofa.ReferenceEllipsoids.GRS80, xyz, out e, out p, out h);

            viv(j, 0, "iauGc2gd", "j2", ref status);
            vvd(e, 0.9827937232473290680, 1e-14, "iauGc2gd", "e2", ref status);
            vvd(p, 0.97160184820607853, 1e-14, "iauGc2gd", "p2", ref status);
            vvd(h, 331.41731754844348, 1e-8, "iauGc2gd", "h2", ref status);

            j = Sofa.iauGc2gd(Sofa.ReferenceEllipsoids.WGS72, xyz, out e, out p, out h);

            viv(j, 0, "iauGc2gd", "j3", ref status);
            vvd(e, 0.9827937232473290680, 1e-14, "iauGc2gd", "e3", ref status);
            vvd(p, 0.9716018181101511937, 1e-14, "iauGc2gd", "p3", ref status);
            vvd(h, 333.2770726130318123, 1e-8, "iauGc2gd", "h3", ref status);

            //j = Sofa.iauGc2gd(4, xyz, out e, out p, out h);

            //viv(j, -1, "iauGc2gd", "j4", ref status);
        }

        static void t_gc2gde(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g c 2 g d e
        **  - - - - - - - - -
        **
        **  Test iauGc2gde function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGc2gde, viv, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            int j;
            double a = 6378136.0, f = 0.0033528;
            double[] xyz = { 2e6, 3e6, 5.244e6 };
            double e, p, h;

            j = Sofa.iauGc2gde(a, f, xyz, out e, out p, out h);

            viv(j, 0, "iauGc2gde", "j", ref status);
            vvd(e, 0.9827937232473290680, 1e-14, "iauGc2gde", "e", ref status);
            vvd(p, 0.9716018377570411532, 1e-14, "iauGc2gde", "p", ref status);
            vvd(h, 332.36862495764397, 1e-8, "iauGc2gde", "h", ref status);
        }

        static void t_gd2gc(ref int status)
        /*
        **  - - - - - - - -
        **   t _ g d 2 g c
        **  - - - - - - - -
        **
        **  Test iauGd2gc function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGd2gc, viv, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            int j;
            double e = 3.1, p = -0.5, h = 2500.0;
            double[] xyz;

            j = Sofa.iauGd2gc(0, e, p, h, out xyz);

            viv(j, -1, "iauGd2gc", "j0", ref status);

            j = Sofa.iauGd2gc(Sofa.ReferenceEllipsoids.WGS84, e, p, h, out xyz);

            viv(j, 0, "iauGd2gc", "j1", ref status);
            vvd(xyz[0], -5599000.5577049947, 1e-7, "iauGd2gc", "1/1", ref status);
            vvd(xyz[1], 233011.67223479203, 1e-7, "iauGd2gc", "2/1", ref status);
            vvd(xyz[2], -3040909.4706983363, 1e-7, "iauGd2gc", "3/1", ref status);

            j = Sofa.iauGd2gc(Sofa.ReferenceEllipsoids.GRS80, e, p, h, out xyz);

            viv(j, 0, "iauGd2gc", "j2", ref status);
            vvd(xyz[0], -5599000.5577260984, 1e-7, "iauGd2gc", "1/2", ref status);
            vvd(xyz[1], 233011.6722356702949, 1e-7, "iauGd2gc", "2/2", ref status);
            vvd(xyz[2], -3040909.4706095476, 1e-7, "iauGd2gc", "3/2", ref status);

            j = Sofa.iauGd2gc(Sofa.ReferenceEllipsoids.WGS72, e, p, h, out xyz);

            viv(j, 0, "iauGd2gc", "j3", ref status);
            vvd(xyz[0], -5598998.7626301490, 1e-7, "iauGd2gc", "1/3", ref status);
            vvd(xyz[1], 233011.5975297822211, 1e-7, "iauGd2gc", "2/3", ref status);
            vvd(xyz[2], -3040908.6861467111, 1e-7, "iauGd2gc", "3/3", ref status);

            //j = Sofa.iauGd2gc(4, e, p, h, xyz);

            //viv(j, -1, "iauGd2gc", "j4", ref status);
        }

        static void t_gd2gce(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g d 2 g c e
        **  - - - - - - - - -
        **
        **  Test iauGd2gce function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGd2gce, viv, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            int j;
            double a = 6378136.0, f = 0.0033528;
            double e = 3.1, p = -0.5, h = 2500.0;
            double[] xyz;

            j = Sofa.iauGd2gce(a, f, e, p, h, out xyz);

            viv(j, 0, "iauGd2gce", "j", ref status);
            vvd(xyz[0], -5598999.6665116328, 1e-7, "iauGd2gce", "1", ref status);
            vvd(xyz[1], 233011.6351463057189, 1e-7, "iauGd2gce", "2", ref status);
            vvd(xyz[2], -3040909.0517314132, 1e-7, "iauGd2gce", "3", ref status);
        }

        static void t_gmst00(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g m s t 0 0
        **  - - - - - - - - -
        **
        **  Test iauGmst00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGmst00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;


            theta = Sofa.iauGmst00(2400000.5, 53736.0, 2400000.5, 53736.0);

            vvd(theta, 1.754174972210740592, 1e-12, "iauGmst00", "", ref status);

        }

        static void t_gmst06(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g m s t 0 6
        **  - - - - - - - - -
        **
        **  Test iauGmst06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGmst06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;


            theta = Sofa.iauGmst06(2400000.5, 53736.0, 2400000.5, 53736.0);

            vvd(theta, 1.754174971870091203, 1e-12, "iauGmst06", "", ref status);

        }

        static void t_gmst82(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g m s t 8 2
        **  - - - - - - - - -
        **
        **  Test iauGmst82 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGmst82, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;


            theta = Sofa.iauGmst82(2400000.5, 53736.0);

            vvd(theta, 1.754174981860675096, 1e-12, "iauGmst82", "", ref status);

        }

        static void t_gst00a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g s t 0 0 a
        **  - - - - - - - - -
        **
        **  Test iauGst00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGst00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;


            theta = Sofa.iauGst00a(2400000.5, 53736.0, 2400000.5, 53736.0);

            vvd(theta, 1.754166138018281369, 1e-12, "iauGst00a", "", ref status);

        }

        static void t_gst00b(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g s t 0 0 b
        **  - - - - - - - - -
        **
        **  Test iauGst00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGst00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;


            theta = Sofa.iauGst00b(2400000.5, 53736.0);

            vvd(theta, 1.754166136510680589, 1e-12, "iauGst00b", "", ref status);

        }

        static void t_gst06(ref int status)
        /*
        **  - - - - - - - -
        **   t _ g s t 0 6
        **  - - - - - - - -
        **
        **  Test iauGst06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGst06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rnpb = { new double[3], new double[3], new double[3] };
            double theta;


            rnpb[0][0] = 0.9999989440476103608;
            rnpb[0][1] = -0.1332881761240011518e-2;
            rnpb[0][2] = -0.5790767434730085097e-3;

            rnpb[1][0] = 0.1332858254308954453e-2;
            rnpb[1][1] = 0.9999991109044505944;
            rnpb[1][2] = -0.4097782710401555759e-4;

            rnpb[2][0] = 0.5791308472168153320e-3;
            rnpb[2][1] = 0.4020595661593994396e-4;
            rnpb[2][2] = 0.9999998314954572365;

            theta = Sofa.iauGst06(2400000.5, 53736.0, 2400000.5, 53736.0, rnpb);

            vvd(theta, 1.754166138018167568, 1e-12, "iauGst06", "", ref status);

        }

        static void t_gst06a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ g s t 0 6 a
        **  - - - - - - - - -
        **
        **  Test iauGst06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGst06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;


            theta = Sofa.iauGst06a(2400000.5, 53736.0, 2400000.5, 53736.0);

            vvd(theta, 1.754166137675019159, 1e-12, "iauGst06a", "", ref status);

        }

        static void t_gst94(ref int status)
        /*
        **  - - - - - - - -
        **   t _ g s t 9 4
        **  - - - - - - - -
        **
        **  Test iauGst94 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauGst94, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;


            theta = Sofa.iauGst94(2400000.5, 53736.0);

            vvd(theta, 1.754166136020645203, 1e-12, "iauGst94", "", ref status);

        }

        static void t_icrs2g(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ i c r s 2 g
        **  - - - - - - - - -
        **
        **  Test iauIcrs2g function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauIcrs2g, vvd
        **
        **  This revision:  2015 January 30
        */
        {
            double dr, dd, dl, db;

            dr = 5.9338074302227188048671087;
            dd = -1.1784870613579944551540570;
            Sofa.iauIcrs2g(dr, dd, out dl, out db);
            vvd(dl, 5.5850536063818546461558, 1e-14, "iauIcrs2g", "L", ref status);
            vvd(db, -0.7853981633974483096157, 1e-14, "iauIcrs2g", "B", ref status);
        }

        static void t_h2fk5(ref int status)
        /*
        **  - - - - - - - -
        **   t _ h 2 f k 5
        **  - - - - - - - -
        **
        **  Test iauH2fk5 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauH2fk5, vvd
        **
        **  This revision:  2017 January 3
        */
        {
            double rh, dh, drh, ddh, pxh, rvh, r5, d5, dr5, dd5, px5, rv5;


            rh = 1.767794352;
            dh = -0.2917512594;
            drh = -2.76413026e-6;
            ddh = -5.92994449e-6;
            pxh = 0.379210;
            rvh = -7.6;

            Sofa.iauH2fk5(rh, dh, drh, ddh, pxh, rvh,
                     out r5, out d5, out dr5, out dd5, out px5, out rv5);

            vvd(r5, 1.767794455700065506, 1e-13,
                "iauH2fk5", "ra", ref status);
            vvd(d5, -0.2917513626469638890, 1e-13,
                "iauH2fk5", "dec", ref status);
            vvd(dr5, -0.27597945024511204e-5, 1e-18,
                "iauH2fk5", "dr5", ref status);
            vvd(dd5, -0.59308014093262838e-5, 1e-18,
                "iauH2fk5", "dd5", ref status);
            vvd(px5, 0.37921, 1e-13,
                "iauH2fk5", "px", ref status);
            vvd(rv5, -7.6000001309071126, 1e-11,
                "iauH2fk5", "rv", ref status);

        }

        static void t_hd2ae(ref int status)
        /*
        **  - - - - - - - -
        **   t _ h d 2 a e
        **  - - - - - - - -
        **
        **  Test iauHd2ae function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauHd2ae and vvd
        **
        **  This revision:  2017 October 21
        */
        {
            double h, d, p, a, e;


            h = 1.1;
            d = 1.2;
            p = 0.3;

            Sofa.iauHd2ae(h, d, p, out a, out e);

            vvd(a, 5.916889243730066194, 1e-13, "iauHd2ae", "a", ref status);
            vvd(e, 0.4472186304990486228, 1e-14, "iauHd2ae", "e", ref status);

        }

        static void t_hd2pa(ref int status)
        /*
        **  - - - - - - - -
        **   t _ h d 2 p a
        **  - - - - - - - -
        **
        **  Test iauHd2pa function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauHd2pa and vvd
        **
        **  This revision:  2017 October 21
        */
        {
            double h, d, p, q;


            h = 1.1;
            d = 1.2;
            p = 0.3;

            q = Sofa.iauHd2pa(h, d, p);

            vvd(q, 1.906227428001995580, 1e-13, "iauHd2pa", "q", ref status);

        }

        static void t_hfk5z(ref int status)
        /*
        **  - - - - - - - -
        **   t _ h f k 5 z
        **  - - - - - - - -
        **
        **  Test iauHfk5z function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauHfk5z, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double rh, dh, r5, d5, dr5, dd5;



            rh = 1.767794352;
            dh = -0.2917512594;

            Sofa.iauHfk5z(rh, dh, 2400000.5, 54479.0, out r5, out d5, out dr5, out dd5);

            vvd(r5, 1.767794490535581026, 1e-13,
                "iauHfk5z", "ra", ref status);
            vvd(d5, -0.2917513695320114258, 1e-14,
                "iauHfk5z", "dec", ref status);
            vvd(dr5, 0.4335890983539243029e-8, 1e-22,
                "iauHfk5z", "dr5", ref status);
            vvd(dd5, -0.8569648841237745902e-9, 1e-23,
                "iauHfk5z", "dd5", ref status);

        }

        static void t_ir(ref int status)
        /*
        **  - - - - -
        **   t _ i r
        **  - - - - -
        **
        **  Test iauIr function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauIr, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r = { new double[3], new double[3], new double[3] };


            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            Sofa.iauIr(out r);

            vvd(r[0][0], 1.0, 0.0, "iauIr", "11", ref status);
            vvd(r[0][1], 0.0, 0.0, "iauIr", "12", ref status);
            vvd(r[0][2], 0.0, 0.0, "iauIr", "13", ref status);

            vvd(r[1][0], 0.0, 0.0, "iauIr", "21", ref status);
            vvd(r[1][1], 1.0, 0.0, "iauIr", "22", ref status);
            vvd(r[1][2], 0.0, 0.0, "iauIr", "23", ref status);

            vvd(r[2][0], 0.0, 0.0, "iauIr", "31", ref status);
            vvd(r[2][1], 0.0, 0.0, "iauIr", "32", ref status);
            vvd(r[2][2], 1.0, 0.0, "iauIr", "33", ref status);

        }

        static void t_jd2cal(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ j d 2 c a l
        **  - - - - - - - - -
        **
        **  Test iauJd2cal function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauJd2cal, viv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dj1, dj2, fd;
            int iy, im, id, j;


            dj1 = 2400000.5;
            dj2 = 50123.9999;

            j = Sofa.iauJd2cal(dj1, dj2, out iy, out im, out id, out fd);

            viv(iy, 1996, "iauJd2cal", "y", ref status);
            viv(im, 2, "iauJd2cal", "m", ref status);
            viv(id, 10, "iauJd2cal", "d", ref status);
            vvd(fd, 0.9999, 1e-7, "iauJd2cal", "fd", ref status);
            viv(j, 0, "iauJd2cal", "j", ref status);

        }

        static void t_jdcalf(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ j d c a l f
        **  - - - - - - - - -
        **
        **  Test iauJdcalf function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauJdcalf, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double dj1, dj2;
            int[] iydmf;
            int j;


            dj1 = 2400000.5;
            dj2 = 50123.9999;

            j = Sofa.iauJdcalf(4, dj1, dj2, out iydmf);

            viv(iydmf[0], 1996, "iauJdcalf", "y", ref status);
            viv(iydmf[1], 2, "iauJdcalf", "m", ref status);
            viv(iydmf[2], 10, "iauJdcalf", "d", ref status);
            viv(iydmf[3], 9999, "iauJdcalf", "f", ref status);

            viv(j, 0, "iauJdcalf", "j", ref status);

        }

        static void t_ld(ref int status)
        /*
        **  - - - - -
        **   t _ l d
        **  - - - - -
        **
        **  Test iauLd function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLd, vvd
        *
        **  This revision:  2013 October 2
        */
        {
            double bm, em, dlim;
            double[] p = new double[3];
            double[] q = new double[3];
            double[] e = new double[3];
            double[] p1 = new double[3];

            bm = 0.00028574;
            p[0] = -0.763276255;
            p[1] = -0.608633767;
            p[2] = -0.216735543;
            q[0] = -0.763276255;
            q[1] = -0.608633767;
            q[2] = -0.216735543;
            e[0] = 0.76700421;
            e[1] = 0.605629598;
            e[2] = 0.211937094;
            em = 8.91276983;
            dlim = 3e-10;

            Sofa.iauLd(bm, p, q, e, em, dlim, ref p1);

            vvd(p1[0], -0.7632762548968159627, 1e-12,
                        "iauLd", "1", ref status);
            vvd(p1[1], -0.6086337670823762701, 1e-12,
                        "iauLd", "2", ref status);
            vvd(p1[2], -0.2167355431320546947, 1e-12,
                        "iauLd", "3", ref status);

        }

        static void t_ldn(ref int status)
        /*
        **  - - - - - -
        **   t _ l d n
        **  - - - - - -
        **
        **  Test iauLdn function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLdn, vvd
        **
        **  This revision:  2013 October 2
        */
        {
            int n;
            iauLDBODY[] b = new iauLDBODY[] { new iauLDBODY(), new iauLDBODY(), new iauLDBODY() };
            double[] ob = new double[3];
            double[] sc = new double[3];
            double[] sn = new double[3];


            n = 3;
            b[0].bm = 0.00028574;
            b[0].dl = 3e-10;
            b[0].pv[0][0] = -7.81014427;
            b[0].pv[0][1] = -5.60956681;
            b[0].pv[0][2] = -1.98079819;
            b[0].pv[1][0] = 0.0030723249;
            b[0].pv[1][1] = -0.00406995477;
            b[0].pv[1][2] = -0.00181335842;
            b[1].bm = 0.00095435;
            b[1].dl = 3e-9;
            b[1].pv[0][0] = 0.738098796;
            b[1].pv[0][1] = 4.63658692;
            b[1].pv[0][2] = 1.9693136;
            b[1].pv[1][0] = -0.00755816922;
            b[1].pv[1][1] = 0.00126913722;
            b[1].pv[1][2] = 0.000727999001;
            b[2].bm = 1.0;
            b[2].dl = 6e-6;
            b[2].pv[0][0] = -0.000712174377;
            b[2].pv[0][1] = -0.00230478303;
            b[2].pv[0][2] = -0.00105865966;
            b[2].pv[1][0] = 6.29235213e-6;
            b[2].pv[1][1] = -3.30888387e-7;
            b[2].pv[1][2] = -2.96486623e-7;
            ob[0] = -0.974170437;
            ob[1] = -0.2115201;
            ob[2] = -0.0917583114;
            sc[0] = -0.763276255;
            sc[1] = -0.608633767;
            sc[2] = -0.216735543;

            Sofa.iauLdn(n, b, ob, sc, ref sn);

            vvd(sn[0], -0.7632762579693333866, 1e-12,
                        "iauLdn", "1", ref status);
            vvd(sn[1], -0.6086337636093002660, 1e-12,
                        "iauLdn", "2", ref status);
            vvd(sn[2], -0.2167355420646328159, 1e-12,
                        "iauLdn", "3", ref status);

        }

        static void t_ldsun(ref int status)
        /*
        **  - - - - - - - -
        **   t _ l d s u n
        **  - - - - - - - -
        **
        **  Test iauLdsun function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLdsun, vvd
        **
        **  This revision:  2013 October 2
        */
        {
            double em;
            double[] p = new double[3];
            double[] e = new double[3];
            double[] p1 = new double[3];


            p[0] = -0.763276255;
            p[1] = -0.608633767;
            p[2] = -0.216735543;
            e[0] = -0.973644023;
            e[1] = -0.20925523;
            e[2] = -0.0907169552;
            em = 0.999809214;

            Sofa.iauLdsun(p, e, em, ref p1);

            vvd(p1[0], -0.7632762580731413169, 1e-12,
                        "iauLdsun", "1", ref status);
            vvd(p1[1], -0.6086337635262647900, 1e-12,
                        "iauLdsun", "2", ref status);
            vvd(p1[2], -0.2167355419322321302, 1e-12,
                        "iauLdsun", "3", ref status);

        }

        static void t_lteceq(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ l t e c e q
        **  - - - - - - - - -
        **
        **  Test iauLteceq function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLteceq, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double epj, dl, db, dr, dd;


            epj = 2500.0;
            dl = 1.5;
            db = 0.6;

            Sofa.iauLteceq(epj, dl, db, out dr, out dd);

            vvd(dr, 1.275156021861921167, 1e-14, "iauLteceq", "dr", ref status);
            vvd(dd, 0.9966573543519204791, 1e-14, "iauLteceq", "dd", ref status);

        }

        static void t_ltecm(ref int status)
        /*
        **  - - - - - - - -
        **   t _ l t e c m
        **  - - - - - - - -
        **
        **  Test iauLtecm function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLtecm, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double epj;
            double[][] rm;


            epj = -3000.0;

            Sofa.iauLtecm(epj, out rm);

            vvd(rm[0][0], 0.3564105644859788825, 1e-14,
                "iauLtecm", "rm11", ref status);
            vvd(rm[0][1], 0.8530575738617682284, 1e-14,
                "iauLtecm", "rm12", ref status);
            vvd(rm[0][2], 0.3811355207795060435, 1e-14,
                "iauLtecm", "rm13", ref status);
            vvd(rm[1][0], -0.9343283469640709942, 1e-14,
                "iauLtecm", "rm21", ref status);
            vvd(rm[1][1], 0.3247830597681745976, 1e-14,
                "iauLtecm", "rm22", ref status);
            vvd(rm[1][2], 0.1467872751535940865, 1e-14,
                "iauLtecm", "rm23", ref status);
            vvd(rm[2][0], 0.1431636191201167793e-2, 1e-14,
                "iauLtecm", "rm31", ref status);
            vvd(rm[2][1], -0.4084222566960599342, 1e-14,
                "iauLtecm", "rm32", ref status);
            vvd(rm[2][2], 0.9127919865189030899, 1e-14,
                "iauLtecm", "rm33", ref status);

        }

        static void t_lteqec(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ l t e q e c
        **  - - - - - - - - -
        **
        **  Test iauLteqec function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLteqec, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double epj, dr, dd, dl, db;


            epj = -1500.0;
            dr = 1.234;
            dd = 0.987;

            Sofa.iauLteqec(epj, dr, dd, out dl, out db);

            vvd(dl, 0.5039483649047114859, 1e-14, "iauLteqec", "dl", ref status);
            vvd(db, 0.5848534459726224882, 1e-14, "iauLteqec", "db", ref status);

        }

        static void t_ltp(ref int status)
        /*
        **  - - - - - -
        **   t _ l t p
        **  - - - - - -
        **
        **  Test iauLtp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLtp, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double epj;
            double[][] rp;


            epj = 1666.666;

            Sofa.iauLtp(epj, out rp);

            vvd(rp[0][0], 0.9967044141159213819, 1e-14,
                "iauLtp", "rp11", ref status);
            vvd(rp[0][1], 0.7437801893193210840e-1, 1e-14,
                "iauLtp", "rp12", ref status);
            vvd(rp[0][2], 0.3237624409345603401e-1, 1e-14,
                "iauLtp", "rp13", ref status);
            vvd(rp[1][0], -0.7437802731819618167e-1, 1e-14,
                "iauLtp", "rp21", ref status);
            vvd(rp[1][1], 0.9972293894454533070, 1e-14,
                "iauLtp", "rp22", ref status);
            vvd(rp[1][2], -0.1205768842723593346e-2, 1e-14,
                "iauLtp", "rp23", ref status);
            vvd(rp[2][0], -0.3237622482766575399e-1, 1e-14,
                "iauLtp", "rp31", ref status);
            vvd(rp[2][1], -0.1206286039697609008e-2, 1e-14,
                "iauLtp", "rp32", ref status);
            vvd(rp[2][2], 0.9994750246704010914, 1e-14,
                "iauLtp", "rp33", ref status);

        }

        static void t_ltpb(ref int status)
        /*
        **  - - - - - - -
        **   t _ l t p b
        **  - - - - - - -
        **
        **  Test iauLtpb function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLtpb, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double epj;
            double[][] rpb;


            epj = 1666.666;

            Sofa.iauLtpb(epj, out rpb);

            vvd(rpb[0][0], 0.9967044167723271851, 1e-14,
                "iauLtpb", "rpb11", ref status);
            vvd(rpb[0][1], 0.7437794731203340345e-1, 1e-14,
                "iauLtpb", "rpb12", ref status);
            vvd(rpb[0][2], 0.3237632684841625547e-1, 1e-14,
                "iauLtpb", "rpb13", ref status);
            vvd(rpb[1][0], -0.7437795663437177152e-1, 1e-14,
                "iauLtpb", "rpb21", ref status);
            vvd(rpb[1][1], 0.9972293947500013666, 1e-14,
                "iauLtpb", "rpb22", ref status);
            vvd(rpb[1][2], -0.1205741865911243235e-2, 1e-14,
                "iauLtpb", "rpb23", ref status);
            vvd(rpb[2][0], -0.3237630543224664992e-1, 1e-14,
                "iauLtpb", "rpb31", ref status);
            vvd(rpb[2][1], -0.1206316791076485295e-2, 1e-14,
                "iauLtpb", "rpb32", ref status);
            vvd(rpb[2][2], 0.9994750220222438819, 1e-14,
                "iauLtpb", "rpb33", ref status);

        }

        static void t_ltpecl(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ l t p e c l
        **  - - - - - - - - -
        **
        **  Test iauLtpecl function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLtpecl, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double epj;
            double[] vec;


            epj = -1500.0;

            Sofa.iauLtpecl(epj, out vec);

            vvd(vec[0], 0.4768625676477096525e-3, 1e-14,
                "iauLtpecl", "vec1", ref status);
            vvd(vec[1], -0.4052259533091875112, 1e-14,
                "iauLtpecl", "vec2", ref status);
            vvd(vec[2], 0.9142164401096448012, 1e-14,
                "iauLtpecl", "vec3", ref status);

        }

        static void t_ltpequ(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ l t p e q u
        **  - - - - - - - - -
        **
        **  Test iauLtpequ function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauLtpequ, vvd
        **
        **  This revision:  2016 March 12
        */
        {
            double epj;
            double[] veq;


            epj = -2500.0;

            Sofa.iauLtpequ(epj, out veq);

            vvd(veq[0], -0.3586652560237326659, 1e-14,
                "iauLtpequ", "veq1", ref status);
            vvd(veq[1], -0.1996978910771128475, 1e-14,
                "iauLtpequ", "veq2", ref status);
            vvd(veq[2], 0.9118552442250819624, 1e-14,
                "iauLtpequ", "veq3", ref status);

        }

        static void t_moon98(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ m o o n 9 8
        **  - - - - - - - - -
        **
        **  Test iauMoon98 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauMoon98, vvd, viv
        **
        **  This revision:  2021 April 12
        */
        {
            double[][] pv;


            Sofa.iauMoon98(2400000.5, 43999.9, out pv);

            vvd(pv[0][0], -0.2601295959971044180e-2, 1e-11,
                "iauMoon98", "x 4", ref status);
            vvd(pv[0][1], 0.6139750944302742189e-3, 1e-11,
                "iauMoon98", "y 4", ref status);
            vvd(pv[0][2], 0.2640794528229828909e-3, 1e-11,
                "iauMoon98", "z 4", ref status);

            vvd(pv[1][0], -0.1244321506649895021e-3, 1e-11,
                "iauMoon98", "xd 4", ref status);
            vvd(pv[1][1], -0.5219076942678119398e-3, 1e-11,
                "iauMoon98", "yd 4", ref status);
            vvd(pv[1][2], -0.1716132214378462047e-3, 1e-11,
                "iauMoon98", "zd 4", ref status);

        }

        static void t_num00a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ n u m 0 0 a
        **  - - - - - - - - -
        **
        **  Test iauNum00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNum00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rmatn;


            Sofa.iauNum00a(2400000.5, 53736.0, out rmatn);

            vvd(rmatn[0][0], 0.9999999999536227949, 1e-12,
                "iauNum00a", "11", ref status);
            vvd(rmatn[0][1], 0.8836238544090873336e-5, 1e-12,
                "iauNum00a", "12", ref status);
            vvd(rmatn[0][2], 0.3830835237722400669e-5, 1e-12,
                "iauNum00a", "13", ref status);

            vvd(rmatn[1][0], -0.8836082880798569274e-5, 1e-12,
                "iauNum00a", "21", ref status);
            vvd(rmatn[1][1], 0.9999999991354655028, 1e-12,
                "iauNum00a", "22", ref status);
            vvd(rmatn[1][2], -0.4063240865362499850e-4, 1e-12,
                "iauNum00a", "23", ref status);

            vvd(rmatn[2][0], -0.3831194272065995866e-5, 1e-12,
                "iauNum00a", "31", ref status);
            vvd(rmatn[2][1], 0.4063237480216291775e-4, 1e-12,
                "iauNum00a", "32", ref status);
            vvd(rmatn[2][2], 0.9999999991671660338, 1e-12,
                "iauNum00a", "33", ref status);

        }

        static void t_num00b(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ n u m 0 0 b
        **  - - - - - - - - -
        **
        **  Test iauNum00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNum00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rmatn;

            Sofa.iauNum00b(2400000.5, 53736, out rmatn);

            vvd(rmatn[0][0], 0.9999999999536069682, 1e-12,
                "iauNum00b", "11", ref status);
            vvd(rmatn[0][1], 0.8837746144871248011e-5, 1e-12,
                "iauNum00b", "12", ref status);
            vvd(rmatn[0][2], 0.3831488838252202945e-5, 1e-12,
                "iauNum00b", "13", ref status);

            vvd(rmatn[1][0], -0.8837590456632304720e-5, 1e-12,
                "iauNum00b", "21", ref status);
            vvd(rmatn[1][1], 0.9999999991354692733, 1e-12,
                "iauNum00b", "22", ref status);
            vvd(rmatn[1][2], -0.4063198798559591654e-4, 1e-12,
                "iauNum00b", "23", ref status);

            vvd(rmatn[2][0], -0.3831847930134941271e-5, 1e-12,
                "iauNum00b", "31", ref status);
            vvd(rmatn[2][1], 0.4063195412258168380e-4, 1e-12,
                "iauNum00b", "32", ref status);
            vvd(rmatn[2][2], 0.9999999991671806225, 1e-12,
                "iauNum00b", "33", ref status);

        }

        static void t_num06a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ n u m 0 6 a
        **  - - - - - - - - -
        **
        **  Test iauNum06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNum06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rmatn;

            Sofa.iauNum06a(2400000.5, 53736, out rmatn);

            vvd(rmatn[0][0], 0.9999999999536227668, 1e-12,
                "iauNum06a", "11", ref status);
            vvd(rmatn[0][1], 0.8836241998111535233e-5, 1e-12,
                "iauNum06a", "12", ref status);
            vvd(rmatn[0][2], 0.3830834608415287707e-5, 1e-12,
                "iauNum06a", "13", ref status);

            vvd(rmatn[1][0], -0.8836086334870740138e-5, 1e-12,
                "iauNum06a", "21", ref status);
            vvd(rmatn[1][1], 0.9999999991354657474, 1e-12,
                "iauNum06a", "22", ref status);
            vvd(rmatn[1][2], -0.4063240188248455065e-4, 1e-12,
                "iauNum06a", "23", ref status);

            vvd(rmatn[2][0], -0.3831193642839398128e-5, 1e-12,
                "iauNum06a", "31", ref status);
            vvd(rmatn[2][1], 0.4063236803101479770e-4, 1e-12,
                "iauNum06a", "32", ref status);
            vvd(rmatn[2][2], 0.9999999991671663114, 1e-12,
                "iauNum06a", "33", ref status);

        }

        static void t_numat(ref int status)
        /*
        **  - - - - - - - -
        **   t _ n u m a t
        **  - - - - - - - -
        **
        **  Test iauNumat function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNumat, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double epsa, dpsi, deps;
            double[][] rmatn;


            epsa = 0.4090789763356509900;
            dpsi = -0.9630909107115582393e-5;
            deps = 0.4063239174001678826e-4;

            Sofa.iauNumat(epsa, dpsi, deps, out rmatn);

            vvd(rmatn[0][0], 0.9999999999536227949, 1e-12,
                "iauNumat", "11", ref status);
            vvd(rmatn[0][1], 0.8836239320236250577e-5, 1e-12,
                "iauNumat", "12", ref status);
            vvd(rmatn[0][2], 0.3830833447458251908e-5, 1e-12,
                "iauNumat", "13", ref status);

            vvd(rmatn[1][0], -0.8836083657016688588e-5, 1e-12,
                "iauNumat", "21", ref status);
            vvd(rmatn[1][1], 0.9999999991354654959, 1e-12,
                "iauNumat", "22", ref status);
            vvd(rmatn[1][2], -0.4063240865361857698e-4, 1e-12,
                "iauNumat", "23", ref status);

            vvd(rmatn[2][0], -0.3831192481833385226e-5, 1e-12,
                "iauNumat", "31", ref status);
            vvd(rmatn[2][1], 0.4063237480216934159e-4, 1e-12,
                "iauNumat", "32", ref status);
            vvd(rmatn[2][2], 0.9999999991671660407, 1e-12,
                "iauNumat", "33", ref status);

        }

        static void t_nut00a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ n u t 0 0 a
        **  - - - - - - - - -
        **
        **  Test iauNut00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNut00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps;


            Sofa.iauNut00a(2400000.5, 53736.0, out dpsi, out deps);

            vvd(dpsi, -0.9630909107115518431e-5, 1e-13,
                "iauNut00a", "dpsi", ref status);
            vvd(deps, 0.4063239174001678710e-4, 1e-13,
                "iauNut00a", "deps", ref status);

        }

        static void t_nut00b(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ n u t 0 0 b
        **  - - - - - - - - -
        **
        **  Test iauNut00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNut00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps;


            Sofa.iauNut00b(2400000.5, 53736.0, out dpsi, out deps);

            vvd(dpsi, -0.9632552291148362783e-5, 1e-13,
                "iauNut00b", "dpsi", ref status);
            vvd(deps, 0.4063197106621159367e-4, 1e-13,
                "iauNut00b", "deps", ref status);

        }

        static void t_nut06a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ n u t 0 6 a
        **  - - - - - - - - -
        **
        **  Test iauNut06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNut06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps;


            Sofa.iauNut06a(2400000.5, 53736.0, out dpsi, out deps);

            vvd(dpsi, -0.9630912025820308797e-5, 1e-13,
                "iauNut06a", "dpsi", ref status);
            vvd(deps, 0.4063238496887249798e-4, 1e-13,
                "iauNut06a", "deps", ref status);

        }

        static void t_nut80(ref int status)
        /*
        **  - - - - - - - -
        **   t _ n u t 8 0
        **  - - - - - - - -
        **
        **  Test iauNut80 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNut80, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps;


            Sofa.iauNut80(2400000.5, 53736.0, out dpsi, out deps);

            vvd(dpsi, -0.9643658353226563966e-5, 1e-13,
                "iauNut80", "dpsi", ref status);
            vvd(deps, 0.4060051006879713322e-4, 1e-13,
                "iauNut80", "deps", ref status);

        }

        static void t_nutm80(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ n u t m 8 0
        **  - - - - - - - - -
        **
        **  Test iauNutm80 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauNutm80, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rmatn;


            Sofa.iauNutm80(2400000.5, 53736.0, out rmatn);

            vvd(rmatn[0][0], 0.9999999999534999268, 1e-12,
               "iauNutm80", "11", ref status);
            vvd(rmatn[0][1], 0.8847935789636432161e-5, 1e-12,
               "iauNutm80", "12", ref status);
            vvd(rmatn[0][2], 0.3835906502164019142e-5, 1e-12,
               "iauNutm80", "13", ref status);

            vvd(rmatn[1][0], -0.8847780042583435924e-5, 1e-12,
               "iauNutm80", "21", ref status);
            vvd(rmatn[1][1], 0.9999999991366569963, 1e-12,
               "iauNutm80", "22", ref status);
            vvd(rmatn[1][2], -0.4060052702727130809e-4, 1e-12,
               "iauNutm80", "23", ref status);

            vvd(rmatn[2][0], -0.3836265729708478796e-5, 1e-12,
               "iauNutm80", "31", ref status);
            vvd(rmatn[2][1], 0.4060049308612638555e-4, 1e-12,
               "iauNutm80", "32", ref status);
            vvd(rmatn[2][2], 0.9999999991684415129, 1e-12,
               "iauNutm80", "33", ref status);

        }

        static void t_obl06(ref int status)
        /*
        **  - - - - - - - -
        **   t _ o b l 0 6
        **  - - - - - - - -
        **
        **  Test iauObl06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauObl06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauObl06(2400000.5, 54388.0), 0.4090749229387258204, 1e-14,
                "iauObl06", "", ref status);
        }

        static void t_obl80(ref int status)
        /*
        **  - - - - - - - -
        **   t _ o b l 8 0
        **  - - - - - - - -
        **
        **  Test iauObl80 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauObl80, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double eps0;


            eps0 = Sofa.iauObl80(2400000.5, 54388.0);

            vvd(eps0, 0.4090751347643816218, 1e-14, "iauObl80", "", ref status);

        }

        static void t_p06e(ref int status)
        /*
        **  - - - - - - -
        **   t _ p 0 6 e
        **  - - - - - - -
        **
        **  Test iauP06e function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauP06e, vvd
        **
        **  This revision:  2020 May 30
        */
        {
            double eps0, psia, oma, bpa, bqa, pia, bpia,
                   epsa, chia, za, zetaa, thetaa, pa, gam, phi, psi;


            Sofa.iauP06e(2400000.5, 52541.0, out eps0, out psia, out oma, out bpa,
                    out bqa, out pia, out bpia, out epsa, out chia, out za,
                    out zetaa, out thetaa, out pa, out gam, out phi, out psi);

            vvd(eps0, 0.4090926006005828715, 1e-14,
                "iauP06e", "eps0", ref status);
            vvd(psia, 0.6664369630191613431e-3, 1e-14,
                "iauP06e", "psia", ref status);
            vvd(oma, 0.4090925973783255982, 1e-14,
                "iauP06e", "oma", ref status);
            vvd(bpa, 0.5561149371265209445e-6, 1e-14,
                "iauP06e", "bpa", ref status);
            vvd(bqa, -0.6191517193290621270e-5, 1e-14,
                "iauP06e", "bqa", ref status);
            vvd(pia, 0.6216441751884382923e-5, 1e-14,
                "iauP06e", "pia", ref status);
            vvd(bpia, 3.052014180023779882, 1e-14,
                "iauP06e", "bpia", ref status);
            vvd(epsa, 0.4090864054922431688, 1e-14,
                "iauP06e", "epsa", ref status);
            vvd(chia, 0.1387703379530915364e-5, 1e-14,
                "iauP06e", "chia", ref status);
            vvd(za, 0.2921789846651790546e-3, 1e-14,
                "iauP06e", "za", ref status);
            vvd(zetaa, 0.3178773290332009310e-3, 1e-14,
                "iauP06e", "zetaa", ref status);
            vvd(thetaa, 0.2650932701657497181e-3, 1e-14,
                "iauP06e", "thetaa", ref status);
            vvd(pa, 0.6651637681381016288e-3, 1e-14,
                "iauP06e", "pa", ref status);
            vvd(gam, 0.1398077115963754987e-5, 1e-14,
                "iauP06e", "gam", ref status);
            vvd(phi, 0.4090864090837462602, 1e-14,
                "iauP06e", "phi", ref status);
            vvd(psi, 0.6664464807480920325e-3, 1e-14,
                "iauP06e", "psi", ref status);

        }

        static void t_p2pv(ref int status)
        /*
        **  - - - - - - -
        **   t _ p 2 p v
        **  - - - - - - -
        **
        **  Test iauP2pv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauP2pv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p = new double[3];
            double[][] pv = { new double[3], new double[3] };


            p[0] = 0.25;
            p[1] = 1.2;
            p[2] = 3.0;

            pv[0][0] = 0.3;
            pv[0][1] = 1.2;
            pv[0][2] = -2.5;

            pv[1][0] = -0.5;
            pv[1][1] = 3.1;
            pv[1][2] = 0.9;

            Sofa.iauP2pv(p, out pv);

            vvd(pv[0][0], 0.25, 0.0, "iauP2pv", "p1", ref status);
            vvd(pv[0][1], 1.2, 0.0, "iauP2pv", "p2", ref status);
            vvd(pv[0][2], 3.0, 0.0, "iauP2pv", "p3", ref status);

            vvd(pv[1][0], 0.0, 0.0, "iauP2pv", "v1", ref status);
            vvd(pv[1][1], 0.0, 0.0, "iauP2pv", "v2", ref status);
            vvd(pv[1][2], 0.0, 0.0, "iauP2pv", "v3", ref status);

        }

        static void t_p2s(ref int status)
        /*
        **  - - - - - -
        **   t _ p 2 s
        **  - - - - - -
        **
        **  Test iauP2s function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauP2s, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p = new double[3];
            double theta, phi, r;


            p[0] = 100.0;
            p[1] = -50.0;
            p[2] = 25.0;

            Sofa.iauP2s(p, out theta, out phi, out r);

            vvd(theta, -0.4636476090008061162, 1e-12, "iauP2s", "theta", ref status);
            vvd(phi, 0.2199879773954594463, 1e-12, "iauP2s", "phi", ref status);
            vvd(r, 114.5643923738960002, 1e-9, "iauP2s", "r", ref status);

        }

        static void t_pap(ref int status)
        /*
        **  - - - - - -
        **   t _ p a p
        **  - - - - - -
        **
        **  Test iauPap function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPap, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] a = new double[3];
            double[] b = new double[3];
            double theta;


            a[0] = 1.0;
            a[1] = 0.1;
            a[2] = 0.2;

            b[0] = -3.0;
            b[1] = 1e-3;
            b[2] = 0.2;

            theta = Sofa.iauPap(a, b);

            vvd(theta, 0.3671514267841113674, 1e-12, "iauPap", "", ref status);

        }

        static void t_pas(ref int status)
        /*
        **  - - - - - -
        **   t _ p a s
        **  - - - - - -
        **
        **  Test iauPas function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPas, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double al, ap, bl, bp, theta;


            al = 1.0;
            ap = 0.1;
            bl = 0.2;
            bp = -1.0;

            theta = Sofa.iauPas(al, ap, bl, bp);

            vvd(theta, -2.724544922932270424, 1e-12, "iauPas", "", ref status);

        }

        static void t_pb06(ref int status)
        /*
        **  - - - - - - -
        **   t _ p b 0 6
        **  - - - - - - -
        **
        **  Test iauPb06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPb06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double bzeta, bz, btheta;


            Sofa.iauPb06(2400000.5, 50123.9999, out bzeta, out bz, out btheta);

            vvd(bzeta, -0.5092634016326478238e-3, 1e-12,
                "iauPb06", "bzeta", ref status);
            vvd(bz, -0.3602772060566044413e-3, 1e-12,
                "iauPb06", "bz", ref status);
            vvd(btheta, -0.3779735537167811177e-3, 1e-12,
                "iauPb06", "btheta", ref status);

        }

        static void t_pdp(ref int status)
        /*
        **  - - - - - -
        **   t _ p d p
        **  - - - - - -
        **
        **  Test iauPdp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPdp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] a = new double[3];
            double[] b = new double[3];
            double adb;


            a[0] = 2.0;
            a[1] = 2.0;
            a[2] = 3.0;

            b[0] = 1.0;
            b[1] = 3.0;
            b[2] = 4.0;

            adb = Sofa.iauPdp(a, b);

            vvd(adb, 20, 1e-12, "iauPdp", "", ref status);

        }

        static void t_pfw06(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p f w 0 6
        **  - - - - - - - -
        **
        **  Test iauPfw06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPfw06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double gamb, phib, psib, epsa;


            Sofa.iauPfw06(2400000.5, 50123.9999, out gamb, out phib, out psib, out epsa);

            vvd(gamb, -0.2243387670997995690e-5, 1e-16,
                "iauPfw06", "gamb", ref status);
            vvd(phib, 0.4091014602391312808, 1e-12,
                "iauPfw06", "phib", ref status);
            vvd(psib, -0.9501954178013031895e-3, 1e-14,
                "iauPfw06", "psib", ref status);
            vvd(epsa, 0.4091014316587367491, 1e-12,
                "iauPfw06", "epsa", ref status);

        }

        static void t_plan94(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p l a n 9 4
        **  - - - - - - - - -
        **
        **  Test iauPlan94 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPlan94, vvd, viv
        **
        **  This revision:  2013 October 2
        */
        {
            double[][] pv;
            int j;


            j = Sofa.iauPlan94(2400000.5, 1e6, 0, out pv);

            vvd(pv[0][0], 0.0, 0.0, "iauPlan94", "x 1", ref status);
            vvd(pv[0][1], 0.0, 0.0, "iauPlan94", "y 1", ref status);
            vvd(pv[0][2], 0.0, 0.0, "iauPlan94", "z 1", ref status);

            vvd(pv[1][0], 0.0, 0.0, "iauPlan94", "xd 1", ref status);
            vvd(pv[1][1], 0.0, 0.0, "iauPlan94", "yd 1", ref status);
            vvd(pv[1][2], 0.0, 0.0, "iauPlan94", "zd 1", ref status);

            viv(j, -1, "iauPlan94", "j 1", ref status);

            j = Sofa.iauPlan94(2400000.5, 1e6, 10, out pv);

            viv(j, -1, "iauPlan94", "j 2", ref status);

            j = Sofa.iauPlan94(2400000.5, -320000, 3, out pv);

            vvd(pv[0][0], 0.9308038666832975759, 1e-11,
                "iauPlan94", "x 3", ref status);
            vvd(pv[0][1], 0.3258319040261346000, 1e-11,
                "iauPlan94", "y 3", ref status);
            vvd(pv[0][2], 0.1422794544481140560, 1e-11,
                "iauPlan94", "z 3", ref status);

            vvd(pv[1][0], -0.6429458958255170006e-2, 1e-11,
                "iauPlan94", "xd 3", ref status);
            vvd(pv[1][1], 0.1468570657704237764e-1, 1e-11,
                "iauPlan94", "yd 3", ref status);
            vvd(pv[1][2], 0.6406996426270981189e-2, 1e-11,
                "iauPlan94", "zd 3", ref status);

            viv(j, 1, "iauPlan94", "j 3", ref status);

            j = Sofa.iauPlan94(2400000.5, 43999.9, 1, out pv);

            vvd(pv[0][0], 0.2945293959257430832, 1e-11,
                "iauPlan94", "x 4", ref status);
            vvd(pv[0][1], -0.2452204176601049596, 1e-11,
                "iauPlan94", "y 4", ref status);
            vvd(pv[0][2], -0.1615427700571978153, 1e-11,
                "iauPlan94", "z 4", ref status);

            vvd(pv[1][0], 0.1413867871404614441e-1, 1e-11,
                "iauPlan94", "xd 4", ref status);
            vvd(pv[1][1], 0.1946548301104706582e-1, 1e-11,
                "iauPlan94", "yd 4", ref status);
            vvd(pv[1][2], 0.8929809783898904786e-2, 1e-11,
                "iauPlan94", "zd 4", ref status);

            viv(j, 0, "iauPlan94", "j 4", ref status);

        }

        static void t_pmat00(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p m a t 0 0
        **  - - - - - - - - -
        **
        **  Test iauPmat00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPmat00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rbp;


            Sofa.iauPmat00(2400000.5, 50123.9999, out rbp);

            vvd(rbp[0][0], 0.9999995505175087260, 1e-12,
                "iauPmat00", "11", ref status);
            vvd(rbp[0][1], 0.8695405883617884705e-3, 1e-14,
                "iauPmat00", "12", ref status);
            vvd(rbp[0][2], 0.3779734722239007105e-3, 1e-14,
                "iauPmat00", "13", ref status);

            vvd(rbp[1][0], -0.8695405990410863719e-3, 1e-14,
                "iauPmat00", "21", ref status);
            vvd(rbp[1][1], 0.9999996219494925900, 1e-12,
                "iauPmat00", "22", ref status);
            vvd(rbp[1][2], -0.1360775820404982209e-6, 1e-14,
                "iauPmat00", "23", ref status);

            vvd(rbp[2][0], -0.3779734476558184991e-3, 1e-14,
                "iauPmat00", "31", ref status);
            vvd(rbp[2][1], -0.1925857585832024058e-6, 1e-14,
                "iauPmat00", "32", ref status);
            vvd(rbp[2][2], 0.9999999285680153377, 1e-12,
                "iauPmat00", "33", ref status);

        }

        static void t_pmat06(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p m a t 0 6
        **  - - - - - - - - -
        **
        **  Test iauPmat06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPmat06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rbp;


            Sofa.iauPmat06(2400000.5, 50123.9999, out rbp);

            vvd(rbp[0][0], 0.9999995505176007047, 1e-12,
                "iauPmat06", "11", ref status);
            vvd(rbp[0][1], 0.8695404617348208406e-3, 1e-14,
                "iauPmat06", "12", ref status);
            vvd(rbp[0][2], 0.3779735201865589104e-3, 1e-14,
                "iauPmat06", "13", ref status);

            vvd(rbp[1][0], -0.8695404723772031414e-3, 1e-14,
                "iauPmat06", "21", ref status);
            vvd(rbp[1][1], 0.9999996219496027161, 1e-12,
                "iauPmat06", "22", ref status);
            vvd(rbp[1][2], -0.1361752497080270143e-6, 1e-14,
                "iauPmat06", "23", ref status);

            vvd(rbp[2][0], -0.3779734957034089490e-3, 1e-14,
                "iauPmat06", "31", ref status);
            vvd(rbp[2][1], -0.1924880847894457113e-6, 1e-14,
                "iauPmat06", "32", ref status);
            vvd(rbp[2][2], 0.9999999285679971958, 1e-12,
                "iauPmat06", "33", ref status);

        }

        static void t_pmat76(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p m a t 7 6
        **  - - - - - - - - -
        **
        **  Test iauPmat76 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPmat76, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rmatp;


            Sofa.iauPmat76(2400000.5, 50123.9999, out rmatp);

            vvd(rmatp[0][0], 0.9999995504328350733, 1e-12,
                "iauPmat76", "11", ref status);
            vvd(rmatp[0][1], 0.8696632209480960785e-3, 1e-14,
                "iauPmat76", "12", ref status);
            vvd(rmatp[0][2], 0.3779153474959888345e-3, 1e-14,
                "iauPmat76", "13", ref status);

            vvd(rmatp[1][0], -0.8696632209485112192e-3, 1e-14,
                "iauPmat76", "21", ref status);
            vvd(rmatp[1][1], 0.9999996218428560614, 1e-12,
                "iauPmat76", "22", ref status);
            vvd(rmatp[1][2], -0.1643284776111886407e-6, 1e-14,
                "iauPmat76", "23", ref status);

            vvd(rmatp[2][0], -0.3779153474950335077e-3, 1e-14,
                "iauPmat76", "31", ref status);
            vvd(rmatp[2][1], -0.1643306746147366896e-6, 1e-14,
                "iauPmat76", "32", ref status);
            vvd(rmatp[2][2], 0.9999999285899790119, 1e-12,
                "iauPmat76", "33", ref status);

        }

        static void t_pm(ref int status)
        /*
        **  - - - - -
        **   t _ p m
        **  - - - - -
        **
        **  Test iauPm function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPm, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p = new double[3];
            double r;


            p[0] = 0.3;
            p[1] = 1.2;
            p[2] = -2.5;

            r = Sofa.iauPm(p);

            vvd(r, 2.789265136196270604, 1e-12, "iauPm", "", ref status);

        }

        static void t_pmp(ref int status)
        /*
        **  - - - - - -
        **   t _ p m p
        **  - - - - - -
        **
        **  Test iauPmp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPmp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] a = new double[3];
            double[] b = new double[3];
            double[] amb = new double[3];


            a[0] = 2.0;
            a[1] = 2.0;
            a[2] = 3.0;

            b[0] = 1.0;
            b[1] = 3.0;
            b[2] = 4.0;

            Sofa.iauPmp(a, b, ref amb);

            vvd(amb[0], 1.0, 1e-12, "iauPmp", "0", ref status);
            vvd(amb[1], -1.0, 1e-12, "iauPmp", "1", ref status);
            vvd(amb[2], -1.0, 1e-12, "iauPmp", "2", ref status);

        }

        static void t_pmpx(ref int status)
        /*
        **  - - - - - - -
        **   t _ p m p x
        **  - - - - - - -
        **
        **  Test iauPmpx function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPmpx, vvd
        **
        **  This revision:  2017 March 15
        */
        {
            double rc, dc, pr, pd, px, rv, pmt;
            double[] pob = new double[3];
            double[] pco;


            rc = 1.234;
            dc = 0.789;
            pr = 1e-5;
            pd = -2e-5;
            px = 1e-2;
            rv = 10.0;
            pmt = 8.75;
            pob[0] = 0.9;
            pob[1] = 0.4;
            pob[2] = 0.1;

            Sofa.iauPmpx(rc, dc, pr, pd, px, rv, pmt, pob, out pco);

            vvd(pco[0], 0.2328137623960308438, 1e-12,
                        "iauPmpx", "1", ref status);
            vvd(pco[1], 0.6651097085397855328, 1e-12,
                        "iauPmpx", "2", ref status);
            vvd(pco[2], 0.7095257765896359837, 1e-12,
                        "iauPmpx", "3", ref status);

        }

        static void t_pmsafe(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p m s a f e
        **  - - - - - - - - -
        **
        **  Test iauPmsafe function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPmsafe, vvd, viv
        **
        **  This revision:  2017 March 15
        */
        {
            int j;
            double ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b,
                   ra2, dec2, pmr2, pmd2, px2, rv2;


            ra1 = 1.234;
            dec1 = 0.789;
            pmr1 = 1e-5;
            pmd1 = -2e-5;
            px1 = 1e-2;
            rv1 = 10.0;
            ep1a = 2400000.5;
            ep1b = 48348.5625;
            ep2a = 2400000.5;
            ep2b = 51544.5;

            j = Sofa.iauPmsafe(ra1, dec1, pmr1, pmd1, px1, rv1,
                          ep1a, ep1b, ep2a, ep2b,
                          out ra2, out dec2, out pmr2, out pmd2, out px2, out rv2);

            vvd(ra2, 1.234087484501017061, 1e-12,
                     "iauPmsafe", "ra2", ref status);
            vvd(dec2, 0.7888249982450468567, 1e-12,
                     "iauPmsafe", "dec2", ref status);
            vvd(pmr2, 0.9996457663586073988e-5, 1e-12,
                      "iauPmsafe", "pmr2", ref status);
            vvd(pmd2, -0.2000040085106754565e-4, 1e-16,
                      "iauPmsafe", "pmd2", ref status);
            vvd(px2, 0.9999997295356830666e-2, 1e-12,
                     "iauPmsafe", "px2", ref status);
            vvd(rv2, 10.38468380293920069, 1e-10,
                     "iauPmsafe", "rv2", ref status);
            viv(j, 0, "iauPmsafe", "j", ref status);

        }

        static void t_pn(ref int status)
        /*
        **  - - - - -
        **   t _ p n
        **  - - - - -
        **
        **  Test iauPn function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPn, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double r;
            double[] p = new double[3];
            double[] u = new double[3];


            p[0] = 0.3;
            p[1] = 1.2;
            p[2] = -2.5;

            Sofa.iauPn(p, out r, ref u);

            vvd(r, 2.789265136196270604, 1e-12, "iauPn", "r", ref status);

            vvd(u[0], 0.1075552109073112058, 1e-12, "iauPn", "u1", ref status);
            vvd(u[1], 0.4302208436292448232, 1e-12, "iauPn", "u2", ref status);
            vvd(u[2], -0.8962934242275933816, 1e-12, "iauPn", "u3", ref status);

        }

        static void t_pn00(ref int status)
        /*
        **  - - - - - - -
        **   t _ p n 0 0
        **  - - - - - - -
        **
        **  Test iauPn00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPn00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps, epsa;
            double[][] rb, rp, rbp, rn, rbpn;


            dpsi = -0.9632552291149335877e-5;
            deps = 0.4063197106621141414e-4;

            Sofa.iauPn00(2400000.5, 53736.0, dpsi, deps,
                    out epsa, out rb, out rp, out rbp, out rn, out rbpn);

            vvd(epsa, 0.4090791789404229916, 1e-12, "iauPn00", "epsa", ref status);

            vvd(rb[0][0], 0.9999999999999942498, 1e-12,
                "iauPn00", "rb11", ref status);
            vvd(rb[0][1], -0.7078279744199196626e-7, 1e-18,
                "iauPn00", "rb12", ref status);
            vvd(rb[0][2], 0.8056217146976134152e-7, 1e-18,
                "iauPn00", "rb13", ref status);

            vvd(rb[1][0], 0.7078279477857337206e-7, 1e-18,
                "iauPn00", "rb21", ref status);
            vvd(rb[1][1], 0.9999999999999969484, 1e-12,
                "iauPn00", "rb22", ref status);
            vvd(rb[1][2], 0.3306041454222136517e-7, 1e-18,
                "iauPn00", "rb23", ref status);

            vvd(rb[2][0], -0.8056217380986972157e-7, 1e-18,
                "iauPn00", "rb31", ref status);
            vvd(rb[2][1], -0.3306040883980552500e-7, 1e-18,
                "iauPn00", "rb32", ref status);
            vvd(rb[2][2], 0.9999999999999962084, 1e-12,
                "iauPn00", "rb33", ref status);

            vvd(rp[0][0], 0.9999989300532289018, 1e-12,
                "iauPn00", "rp11", ref status);
            vvd(rp[0][1], -0.1341647226791824349e-2, 1e-14,
                "iauPn00", "rp12", ref status);
            vvd(rp[0][2], -0.5829880927190296547e-3, 1e-14,
                "iauPn00", "rp13", ref status);

            vvd(rp[1][0], 0.1341647231069759008e-2, 1e-14,
                "iauPn00", "rp21", ref status);
            vvd(rp[1][1], 0.9999990999908750433, 1e-12,
                "iauPn00", "rp22", ref status);
            vvd(rp[1][2], -0.3837444441583715468e-6, 1e-14,
                "iauPn00", "rp23", ref status);

            vvd(rp[2][0], 0.5829880828740957684e-3, 1e-14,
                "iauPn00", "rp31", ref status);
            vvd(rp[2][1], -0.3984203267708834759e-6, 1e-14,
                "iauPn00", "rp32", ref status);
            vvd(rp[2][2], 0.9999998300623538046, 1e-12,
                "iauPn00", "rp33", ref status);

            vvd(rbp[0][0], 0.9999989300052243993, 1e-12,
                "iauPn00", "rbp11", ref status);
            vvd(rbp[0][1], -0.1341717990239703727e-2, 1e-14,
                "iauPn00", "rbp12", ref status);
            vvd(rbp[0][2], -0.5829075749891684053e-3, 1e-14,
                "iauPn00", "rbp13", ref status);

            vvd(rbp[1][0], 0.1341718013831739992e-2, 1e-14,
                "iauPn00", "rbp21", ref status);
            vvd(rbp[1][1], 0.9999990998959191343, 1e-12,
                "iauPn00", "rbp22", ref status);
            vvd(rbp[1][2], -0.3505759733565421170e-6, 1e-14,
                "iauPn00", "rbp23", ref status);

            vvd(rbp[2][0], 0.5829075206857717883e-3, 1e-14,
                "iauPn00", "rbp31", ref status);
            vvd(rbp[2][1], -0.4315219955198608970e-6, 1e-14,
                "iauPn00", "rbp32", ref status);
            vvd(rbp[2][2], 0.9999998301093036269, 1e-12,
                "iauPn00", "rbp33", ref status);

            vvd(rn[0][0], 0.9999999999536069682, 1e-12,
                "iauPn00", "rn11", ref status);
            vvd(rn[0][1], 0.8837746144872140812e-5, 1e-16,
                "iauPn00", "rn12", ref status);
            vvd(rn[0][2], 0.3831488838252590008e-5, 1e-16,
                "iauPn00", "rn13", ref status);

            vvd(rn[1][0], -0.8837590456633197506e-5, 1e-16,
                "iauPn00", "rn21", ref status);
            vvd(rn[1][1], 0.9999999991354692733, 1e-12,
                "iauPn00", "rn22", ref status);
            vvd(rn[1][2], -0.4063198798559573702e-4, 1e-16,
                "iauPn00", "rn23", ref status);

            vvd(rn[2][0], -0.3831847930135328368e-5, 1e-16,
                "iauPn00", "rn31", ref status);
            vvd(rn[2][1], 0.4063195412258150427e-4, 1e-16,
                "iauPn00", "rn32", ref status);
            vvd(rn[2][2], 0.9999999991671806225, 1e-12,
                "iauPn00", "rn33", ref status);

            vvd(rbpn[0][0], 0.9999989440499982806, 1e-12,
                "iauPn00", "rbpn11", ref status);
            vvd(rbpn[0][1], -0.1332880253640848301e-2, 1e-14,
                "iauPn00", "rbpn12", ref status);
            vvd(rbpn[0][2], -0.5790760898731087295e-3, 1e-14,
                "iauPn00", "rbpn13", ref status);

            vvd(rbpn[1][0], 0.1332856746979948745e-2, 1e-14,
                "iauPn00", "rbpn21", ref status);
            vvd(rbpn[1][1], 0.9999991109064768883, 1e-12,
                "iauPn00", "rbpn22", ref status);
            vvd(rbpn[1][2], -0.4097740555723063806e-4, 1e-14,
                "iauPn00", "rbpn23", ref status);

            vvd(rbpn[2][0], 0.5791301929950205000e-3, 1e-14,
                "iauPn00", "rbpn31", ref status);
            vvd(rbpn[2][1], 0.4020553681373702931e-4, 1e-14,
                "iauPn00", "rbpn32", ref status);
            vvd(rbpn[2][2], 0.9999998314958529887, 1e-12,
                "iauPn00", "rbpn33", ref status);

        }

        static void t_pn00a(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p n 0 0 a
        **  - - - - - - - -
        **
        **  Test iauPn00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPn00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps, epsa;
            double[][] rb, rp, rbp, rn, rbpn;


            Sofa.iauPn00a(2400000.5, 53736.0,
                     out dpsi, out deps, out epsa, out rb, out rp, out rbp, out rn, out rbpn);

            vvd(dpsi, -0.9630909107115518431e-5, 1e-12,
                "iauPn00a", "dpsi", ref status);
            vvd(deps, 0.4063239174001678710e-4, 1e-12,
                "iauPn00a", "deps", ref status);
            vvd(epsa, 0.4090791789404229916, 1e-12, "iauPn00a", "epsa", ref status);

            vvd(rb[0][0], 0.9999999999999942498, 1e-12,
                "iauPn00a", "rb11", ref status);
            vvd(rb[0][1], -0.7078279744199196626e-7, 1e-16,
                "iauPn00a", "rb12", ref status);
            vvd(rb[0][2], 0.8056217146976134152e-7, 1e-16,
                "iauPn00a", "rb13", ref status);

            vvd(rb[1][0], 0.7078279477857337206e-7, 1e-16,
                "iauPn00a", "rb21", ref status);
            vvd(rb[1][1], 0.9999999999999969484, 1e-12,
                "iauPn00a", "rb22", ref status);
            vvd(rb[1][2], 0.3306041454222136517e-7, 1e-16,
                "iauPn00a", "rb23", ref status);

            vvd(rb[2][0], -0.8056217380986972157e-7, 1e-16,
                "iauPn00a", "rb31", ref status);
            vvd(rb[2][1], -0.3306040883980552500e-7, 1e-16,
                "iauPn00a", "rb32", ref status);
            vvd(rb[2][2], 0.9999999999999962084, 1e-12,
                "iauPn00a", "rb33", ref status);

            vvd(rp[0][0], 0.9999989300532289018, 1e-12,
                "iauPn00a", "rp11", ref status);
            vvd(rp[0][1], -0.1341647226791824349e-2, 1e-14,
                "iauPn00a", "rp12", ref status);
            vvd(rp[0][2], -0.5829880927190296547e-3, 1e-14,
                "iauPn00a", "rp13", ref status);

            vvd(rp[1][0], 0.1341647231069759008e-2, 1e-14,
                "iauPn00a", "rp21", ref status);
            vvd(rp[1][1], 0.9999990999908750433, 1e-12,
                "iauPn00a", "rp22", ref status);
            vvd(rp[1][2], -0.3837444441583715468e-6, 1e-14,
                "iauPn00a", "rp23", ref status);

            vvd(rp[2][0], 0.5829880828740957684e-3, 1e-14,
                "iauPn00a", "rp31", ref status);
            vvd(rp[2][1], -0.3984203267708834759e-6, 1e-14,
                "iauPn00a", "rp32", ref status);
            vvd(rp[2][2], 0.9999998300623538046, 1e-12,
                "iauPn00a", "rp33", ref status);

            vvd(rbp[0][0], 0.9999989300052243993, 1e-12,
                "iauPn00a", "rbp11", ref status);
            vvd(rbp[0][1], -0.1341717990239703727e-2, 1e-14,
                "iauPn00a", "rbp12", ref status);
            vvd(rbp[0][2], -0.5829075749891684053e-3, 1e-14,
                "iauPn00a", "rbp13", ref status);

            vvd(rbp[1][0], 0.1341718013831739992e-2, 1e-14,
                "iauPn00a", "rbp21", ref status);
            vvd(rbp[1][1], 0.9999990998959191343, 1e-12,
                "iauPn00a", "rbp22", ref status);
            vvd(rbp[1][2], -0.3505759733565421170e-6, 1e-14,
                "iauPn00a", "rbp23", ref status);

            vvd(rbp[2][0], 0.5829075206857717883e-3, 1e-14,
                "iauPn00a", "rbp31", ref status);
            vvd(rbp[2][1], -0.4315219955198608970e-6, 1e-14,
                "iauPn00a", "rbp32", ref status);
            vvd(rbp[2][2], 0.9999998301093036269, 1e-12,
                "iauPn00a", "rbp33", ref status);

            vvd(rn[0][0], 0.9999999999536227949, 1e-12,
                "iauPn00a", "rn11", ref status);
            vvd(rn[0][1], 0.8836238544090873336e-5, 1e-14,
                "iauPn00a", "rn12", ref status);
            vvd(rn[0][2], 0.3830835237722400669e-5, 1e-14,
                "iauPn00a", "rn13", ref status);

            vvd(rn[1][0], -0.8836082880798569274e-5, 1e-14,
                "iauPn00a", "rn21", ref status);
            vvd(rn[1][1], 0.9999999991354655028, 1e-12,
                "iauPn00a", "rn22", ref status);
            vvd(rn[1][2], -0.4063240865362499850e-4, 1e-14,
                "iauPn00a", "rn23", ref status);

            vvd(rn[2][0], -0.3831194272065995866e-5, 1e-14,
                "iauPn00a", "rn31", ref status);
            vvd(rn[2][1], 0.4063237480216291775e-4, 1e-14,
                "iauPn00a", "rn32", ref status);
            vvd(rn[2][2], 0.9999999991671660338, 1e-12,
                "iauPn00a", "rn33", ref status);

            vvd(rbpn[0][0], 0.9999989440476103435, 1e-12,
                "iauPn00a", "rbpn11", ref status);
            vvd(rbpn[0][1], -0.1332881761240011763e-2, 1e-14,
                "iauPn00a", "rbpn12", ref status);
            vvd(rbpn[0][2], -0.5790767434730085751e-3, 1e-14,
                "iauPn00a", "rbpn13", ref status);

            vvd(rbpn[1][0], 0.1332858254308954658e-2, 1e-14,
                "iauPn00a", "rbpn21", ref status);
            vvd(rbpn[1][1], 0.9999991109044505577, 1e-12,
                "iauPn00a", "rbpn22", ref status);
            vvd(rbpn[1][2], -0.4097782710396580452e-4, 1e-14,
                "iauPn00a", "rbpn23", ref status);

            vvd(rbpn[2][0], 0.5791308472168152904e-3, 1e-14,
                "iauPn00a", "rbpn31", ref status);
            vvd(rbpn[2][1], 0.4020595661591500259e-4, 1e-14,
                "iauPn00a", "rbpn32", ref status);
            vvd(rbpn[2][2], 0.9999998314954572304, 1e-12,
                "iauPn00a", "rbpn33", ref status);

        }

        static void t_pn00b(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p n 0 0 b
        **  - - - - - - - -
        **
        **  Test iauPn00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPn00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps, epsa;
            double[][] rb, rp, rbp, rn, rbpn;


            Sofa.iauPn00b(2400000.5, 53736.0, out dpsi, out deps, out epsa,
                     out rb, out rp, out rbp, out rn, out rbpn);

            vvd(dpsi, -0.9632552291148362783e-5, 1e-12,
                "iauPn00b", "dpsi", ref status);
            vvd(deps, 0.4063197106621159367e-4, 1e-12,
                "iauPn00b", "deps", ref status);
            vvd(epsa, 0.4090791789404229916, 1e-12, "iauPn00b", "epsa", ref status);

            vvd(rb[0][0], 0.9999999999999942498, 1e-12,
               "iauPn00b", "rb11", ref status);
            vvd(rb[0][1], -0.7078279744199196626e-7, 1e-16,
               "iauPn00b", "rb12", ref status);
            vvd(rb[0][2], 0.8056217146976134152e-7, 1e-16,
               "iauPn00b", "rb13", ref status);

            vvd(rb[1][0], 0.7078279477857337206e-7, 1e-16,
               "iauPn00b", "rb21", ref status);
            vvd(rb[1][1], 0.9999999999999969484, 1e-12,
               "iauPn00b", "rb22", ref status);
            vvd(rb[1][2], 0.3306041454222136517e-7, 1e-16,
               "iauPn00b", "rb23", ref status);

            vvd(rb[2][0], -0.8056217380986972157e-7, 1e-16,
               "iauPn00b", "rb31", ref status);
            vvd(rb[2][1], -0.3306040883980552500e-7, 1e-16,
               "iauPn00b", "rb32", ref status);
            vvd(rb[2][2], 0.9999999999999962084, 1e-12,
               "iauPn00b", "rb33", ref status);

            vvd(rp[0][0], 0.9999989300532289018, 1e-12,
               "iauPn00b", "rp11", ref status);
            vvd(rp[0][1], -0.1341647226791824349e-2, 1e-14,
               "iauPn00b", "rp12", ref status);
            vvd(rp[0][2], -0.5829880927190296547e-3, 1e-14,
               "iauPn00b", "rp13", ref status);

            vvd(rp[1][0], 0.1341647231069759008e-2, 1e-14,
               "iauPn00b", "rp21", ref status);
            vvd(rp[1][1], 0.9999990999908750433, 1e-12,
               "iauPn00b", "rp22", ref status);
            vvd(rp[1][2], -0.3837444441583715468e-6, 1e-14,
               "iauPn00b", "rp23", ref status);

            vvd(rp[2][0], 0.5829880828740957684e-3, 1e-14,
               "iauPn00b", "rp31", ref status);
            vvd(rp[2][1], -0.3984203267708834759e-6, 1e-14,
               "iauPn00b", "rp32", ref status);
            vvd(rp[2][2], 0.9999998300623538046, 1e-12,
               "iauPn00b", "rp33", ref status);

            vvd(rbp[0][0], 0.9999989300052243993, 1e-12,
               "iauPn00b", "rbp11", ref status);
            vvd(rbp[0][1], -0.1341717990239703727e-2, 1e-14,
               "iauPn00b", "rbp12", ref status);
            vvd(rbp[0][2], -0.5829075749891684053e-3, 1e-14,
               "iauPn00b", "rbp13", ref status);

            vvd(rbp[1][0], 0.1341718013831739992e-2, 1e-14,
               "iauPn00b", "rbp21", ref status);
            vvd(rbp[1][1], 0.9999990998959191343, 1e-12,
               "iauPn00b", "rbp22", ref status);
            vvd(rbp[1][2], -0.3505759733565421170e-6, 1e-14,
               "iauPn00b", "rbp23", ref status);

            vvd(rbp[2][0], 0.5829075206857717883e-3, 1e-14,
               "iauPn00b", "rbp31", ref status);
            vvd(rbp[2][1], -0.4315219955198608970e-6, 1e-14,
               "iauPn00b", "rbp32", ref status);
            vvd(rbp[2][2], 0.9999998301093036269, 1e-12,
               "iauPn00b", "rbp33", ref status);

            vvd(rn[0][0], 0.9999999999536069682, 1e-12,
               "iauPn00b", "rn11", ref status);
            vvd(rn[0][1], 0.8837746144871248011e-5, 1e-14,
               "iauPn00b", "rn12", ref status);
            vvd(rn[0][2], 0.3831488838252202945e-5, 1e-14,
               "iauPn00b", "rn13", ref status);

            vvd(rn[1][0], -0.8837590456632304720e-5, 1e-14,
               "iauPn00b", "rn21", ref status);
            vvd(rn[1][1], 0.9999999991354692733, 1e-12,
               "iauPn00b", "rn22", ref status);
            vvd(rn[1][2], -0.4063198798559591654e-4, 1e-14,
               "iauPn00b", "rn23", ref status);

            vvd(rn[2][0], -0.3831847930134941271e-5, 1e-14,
               "iauPn00b", "rn31", ref status);
            vvd(rn[2][1], 0.4063195412258168380e-4, 1e-14,
               "iauPn00b", "rn32", ref status);
            vvd(rn[2][2], 0.9999999991671806225, 1e-12,
               "iauPn00b", "rn33", ref status);

            vvd(rbpn[0][0], 0.9999989440499982806, 1e-12,
               "iauPn00b", "rbpn11", ref status);
            vvd(rbpn[0][1], -0.1332880253640849194e-2, 1e-14,
               "iauPn00b", "rbpn12", ref status);
            vvd(rbpn[0][2], -0.5790760898731091166e-3, 1e-14,
               "iauPn00b", "rbpn13", ref status);

            vvd(rbpn[1][0], 0.1332856746979949638e-2, 1e-14,
               "iauPn00b", "rbpn21", ref status);
            vvd(rbpn[1][1], 0.9999991109064768883, 1e-12,
               "iauPn00b", "rbpn22", ref status);
            vvd(rbpn[1][2], -0.4097740555723081811e-4, 1e-14,
               "iauPn00b", "rbpn23", ref status);

            vvd(rbpn[2][0], 0.5791301929950208873e-3, 1e-14,
               "iauPn00b", "rbpn31", ref status);
            vvd(rbpn[2][1], 0.4020553681373720832e-4, 1e-14,
               "iauPn00b", "rbpn32", ref status);
            vvd(rbpn[2][2], 0.9999998314958529887, 1e-12,
               "iauPn00b", "rbpn33", ref status);

        }

        static void t_pn06a(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p n 0 6 a
        **  - - - - - - - -
        **
        **  Test iauPn06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPn06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps, epsa;
            double[][] rb, rp, rbp, rn, rbpn;


            Sofa.iauPn06a(2400000.5, 53736.0, out dpsi, out deps, out epsa,
                     out rb, out rp, out rbp, out rn, out rbpn);

            vvd(dpsi, -0.9630912025820308797e-5, 1e-12,
                "iauPn06a", "dpsi", ref status);
            vvd(deps, 0.4063238496887249798e-4, 1e-12,
                "iauPn06a", "deps", ref status);
            vvd(epsa, 0.4090789763356509926, 1e-12, "iauPn06a", "epsa", ref status);

            vvd(rb[0][0], 0.9999999999999942497, 1e-12,
                "iauPn06a", "rb11", ref status);
            vvd(rb[0][1], -0.7078368960971557145e-7, 1e-14,
                "iauPn06a", "rb12", ref status);
            vvd(rb[0][2], 0.8056213977613185606e-7, 1e-14,
                "iauPn06a", "rb13", ref status);

            vvd(rb[1][0], 0.7078368694637674333e-7, 1e-14,
                "iauPn06a", "rb21", ref status);
            vvd(rb[1][1], 0.9999999999999969484, 1e-12,
                "iauPn06a", "rb22", ref status);
            vvd(rb[1][2], 0.3305943742989134124e-7, 1e-14,
                "iauPn06a", "rb23", ref status);

            vvd(rb[2][0], -0.8056214211620056792e-7, 1e-14,
                "iauPn06a", "rb31", ref status);
            vvd(rb[2][1], -0.3305943172740586950e-7, 1e-14,
                "iauPn06a", "rb32", ref status);
            vvd(rb[2][2], 0.9999999999999962084, 1e-12,
                "iauPn06a", "rb33", ref status);

            vvd(rp[0][0], 0.9999989300536854831, 1e-12,
                "iauPn06a", "rp11", ref status);
            vvd(rp[0][1], -0.1341646886204443795e-2, 1e-14,
                "iauPn06a", "rp12", ref status);
            vvd(rp[0][2], -0.5829880933488627759e-3, 1e-14,
                "iauPn06a", "rp13", ref status);

            vvd(rp[1][0], 0.1341646890569782183e-2, 1e-14,
                "iauPn06a", "rp21", ref status);
            vvd(rp[1][1], 0.9999990999913319321, 1e-12,
                "iauPn06a", "rp22", ref status);
            vvd(rp[1][2], -0.3835944216374477457e-6, 1e-14,
                "iauPn06a", "rp23", ref status);

            vvd(rp[2][0], 0.5829880833027867368e-3, 1e-14,
                "iauPn06a", "rp31", ref status);
            vvd(rp[2][1], -0.3985701514686976112e-6, 1e-14,
                "iauPn06a", "rp32", ref status);
            vvd(rp[2][2], 0.9999998300623534950, 1e-12,
                "iauPn06a", "rp33", ref status);

            vvd(rbp[0][0], 0.9999989300056797893, 1e-12,
                "iauPn06a", "rbp11", ref status);
            vvd(rbp[0][1], -0.1341717650545059598e-2, 1e-14,
                "iauPn06a", "rbp12", ref status);
            vvd(rbp[0][2], -0.5829075756493728856e-3, 1e-14,
                "iauPn06a", "rbp13", ref status);

            vvd(rbp[1][0], 0.1341717674223918101e-2, 1e-14,
                "iauPn06a", "rbp21", ref status);
            vvd(rbp[1][1], 0.9999990998963748448, 1e-12,
                "iauPn06a", "rbp22", ref status);
            vvd(rbp[1][2], -0.3504269280170069029e-6, 1e-14,
                "iauPn06a", "rbp23", ref status);

            vvd(rbp[2][0], 0.5829075211461454599e-3, 1e-14,
                "iauPn06a", "rbp31", ref status);
            vvd(rbp[2][1], -0.4316708436255949093e-6, 1e-14,
                "iauPn06a", "rbp32", ref status);
            vvd(rbp[2][2], 0.9999998301093032943, 1e-12,
                "iauPn06a", "rbp33", ref status);

            vvd(rn[0][0], 0.9999999999536227668, 1e-12,
                "iauPn06a", "rn11", ref status);
            vvd(rn[0][1], 0.8836241998111535233e-5, 1e-14,
                "iauPn06a", "rn12", ref status);
            vvd(rn[0][2], 0.3830834608415287707e-5, 1e-14,
                "iauPn06a", "rn13", ref status);

            vvd(rn[1][0], -0.8836086334870740138e-5, 1e-14,
                "iauPn06a", "rn21", ref status);
            vvd(rn[1][1], 0.9999999991354657474, 1e-12,
                "iauPn06a", "rn22", ref status);
            vvd(rn[1][2], -0.4063240188248455065e-4, 1e-14,
                "iauPn06a", "rn23", ref status);

            vvd(rn[2][0], -0.3831193642839398128e-5, 1e-14,
                "iauPn06a", "rn31", ref status);
            vvd(rn[2][1], 0.4063236803101479770e-4, 1e-14,
                "iauPn06a", "rn32", ref status);
            vvd(rn[2][2], 0.9999999991671663114, 1e-12,
                "iauPn06a", "rn33", ref status);

            vvd(rbpn[0][0], 0.9999989440480669738, 1e-12,
                "iauPn06a", "rbpn11", ref status);
            vvd(rbpn[0][1], -0.1332881418091915973e-2, 1e-14,
                "iauPn06a", "rbpn12", ref status);
            vvd(rbpn[0][2], -0.5790767447612042565e-3, 1e-14,
                "iauPn06a", "rbpn13", ref status);

            vvd(rbpn[1][0], 0.1332857911250989133e-2, 1e-14,
                "iauPn06a", "rbpn21", ref status);
            vvd(rbpn[1][1], 0.9999991109049141908, 1e-12,
                "iauPn06a", "rbpn22", ref status);
            vvd(rbpn[1][2], -0.4097767128546784878e-4, 1e-14,
                "iauPn06a", "rbpn23", ref status);

            vvd(rbpn[2][0], 0.5791308482835292617e-3, 1e-14,
                "iauPn06a", "rbpn31", ref status);
            vvd(rbpn[2][1], 0.4020580099454020310e-4, 1e-14,
                "iauPn06a", "rbpn32", ref status);
            vvd(rbpn[2][2], 0.9999998314954628695, 1e-12,
                "iauPn06a", "rbpn33", ref status);

        }

        static void t_pn06(ref int status)
        /*
        **  - - - - - - -
        **   t _ p n 0 6
        **  - - - - - - -
        **
        **  Test iauPn06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPn06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsi, deps, epsa;
            double[][] rb, rp, rbp, rn, rbpn;


            dpsi = -0.9632552291149335877e-5;
            deps = 0.4063197106621141414e-4;

            Sofa.iauPn06(2400000.5, 53736.0, dpsi, deps,
                    out epsa, out rb, out rp, out rbp, out rn, out rbpn);

            vvd(epsa, 0.4090789763356509926, 1e-12, "iauPn06", "epsa", ref status);

            vvd(rb[0][0], 0.9999999999999942497, 1e-12,
                "iauPn06", "rb11", ref status);
            vvd(rb[0][1], -0.7078368960971557145e-7, 1e-14,
                "iauPn06", "rb12", ref status);
            vvd(rb[0][2], 0.8056213977613185606e-7, 1e-14,
                "iauPn06", "rb13", ref status);

            vvd(rb[1][0], 0.7078368694637674333e-7, 1e-14,
                "iauPn06", "rb21", ref status);
            vvd(rb[1][1], 0.9999999999999969484, 1e-12,
                "iauPn06", "rb22", ref status);
            vvd(rb[1][2], 0.3305943742989134124e-7, 1e-14,
                "iauPn06", "rb23", ref status);

            vvd(rb[2][0], -0.8056214211620056792e-7, 1e-14,
                "iauPn06", "rb31", ref status);
            vvd(rb[2][1], -0.3305943172740586950e-7, 1e-14,
                "iauPn06", "rb32", ref status);
            vvd(rb[2][2], 0.9999999999999962084, 1e-12,
                "iauPn06", "rb33", ref status);

            vvd(rp[0][0], 0.9999989300536854831, 1e-12,
                "iauPn06", "rp11", ref status);
            vvd(rp[0][1], -0.1341646886204443795e-2, 1e-14,
                "iauPn06", "rp12", ref status);
            vvd(rp[0][2], -0.5829880933488627759e-3, 1e-14,
                "iauPn06", "rp13", ref status);

            vvd(rp[1][0], 0.1341646890569782183e-2, 1e-14,
                "iauPn06", "rp21", ref status);
            vvd(rp[1][1], 0.9999990999913319321, 1e-12,
                "iauPn06", "rp22", ref status);
            vvd(rp[1][2], -0.3835944216374477457e-6, 1e-14,
                "iauPn06", "rp23", ref status);

            vvd(rp[2][0], 0.5829880833027867368e-3, 1e-14,
                "iauPn06", "rp31", ref status);
            vvd(rp[2][1], -0.3985701514686976112e-6, 1e-14,
                "iauPn06", "rp32", ref status);
            vvd(rp[2][2], 0.9999998300623534950, 1e-12,
                "iauPn06", "rp33", ref status);

            vvd(rbp[0][0], 0.9999989300056797893, 1e-12,
                "iauPn06", "rbp11", ref status);
            vvd(rbp[0][1], -0.1341717650545059598e-2, 1e-14,
                "iauPn06", "rbp12", ref status);
            vvd(rbp[0][2], -0.5829075756493728856e-3, 1e-14,
                "iauPn06", "rbp13", ref status);

            vvd(rbp[1][0], 0.1341717674223918101e-2, 1e-14,
                "iauPn06", "rbp21", ref status);
            vvd(rbp[1][1], 0.9999990998963748448, 1e-12,
                "iauPn06", "rbp22", ref status);
            vvd(rbp[1][2], -0.3504269280170069029e-6, 1e-14,
                "iauPn06", "rbp23", ref status);

            vvd(rbp[2][0], 0.5829075211461454599e-3, 1e-14,
                "iauPn06", "rbp31", ref status);
            vvd(rbp[2][1], -0.4316708436255949093e-6, 1e-14,
                "iauPn06", "rbp32", ref status);
            vvd(rbp[2][2], 0.9999998301093032943, 1e-12,
                "iauPn06", "rbp33", ref status);

            vvd(rn[0][0], 0.9999999999536069682, 1e-12,
                "iauPn06", "rn11", ref status);
            vvd(rn[0][1], 0.8837746921149881914e-5, 1e-14,
                "iauPn06", "rn12", ref status);
            vvd(rn[0][2], 0.3831487047682968703e-5, 1e-14,
                "iauPn06", "rn13", ref status);

            vvd(rn[1][0], -0.8837591232983692340e-5, 1e-14,
                "iauPn06", "rn21", ref status);
            vvd(rn[1][1], 0.9999999991354692664, 1e-12,
                "iauPn06", "rn22", ref status);
            vvd(rn[1][2], -0.4063198798558931215e-4, 1e-14,
                "iauPn06", "rn23", ref status);

            vvd(rn[2][0], -0.3831846139597250235e-5, 1e-14,
                "iauPn06", "rn31", ref status);
            vvd(rn[2][1], 0.4063195412258792914e-4, 1e-14,
                "iauPn06", "rn32", ref status);
            vvd(rn[2][2], 0.9999999991671806293, 1e-12,
                "iauPn06", "rn33", ref status);

            vvd(rbpn[0][0], 0.9999989440504506688, 1e-12,
                "iauPn06", "rbpn11", ref status);
            vvd(rbpn[0][1], -0.1332879913170492655e-2, 1e-14,
                "iauPn06", "rbpn12", ref status);
            vvd(rbpn[0][2], -0.5790760923225655753e-3, 1e-14,
                "iauPn06", "rbpn13", ref status);

            vvd(rbpn[1][0], 0.1332856406595754748e-2, 1e-14,
                "iauPn06", "rbpn21", ref status);
            vvd(rbpn[1][1], 0.9999991109069366795, 1e-12,
                "iauPn06", "rbpn22", ref status);
            vvd(rbpn[1][2], -0.4097725651142641812e-4, 1e-14,
                "iauPn06", "rbpn23", ref status);

            vvd(rbpn[2][0], 0.5791301952321296716e-3, 1e-14,
                "iauPn06", "rbpn31", ref status);
            vvd(rbpn[2][1], 0.4020538796195230577e-4, 1e-14,
                "iauPn06", "rbpn32", ref status);
            vvd(rbpn[2][2], 0.9999998314958576778, 1e-12,
                "iauPn06", "rbpn33", ref status);

        }

        static void t_pnm00a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p n m 0 0 a
        **  - - - - - - - - -
        **
        **  Test iauPnm00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPnm00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rbpn;


            Sofa.iauPnm00a(2400000.5, 50123.9999, out rbpn);

            vvd(rbpn[0][0], 0.9999995832793134257, 1e-12,
                "iauPnm00a", "11", ref status);
            vvd(rbpn[0][1], 0.8372384254137809439e-3, 1e-14,
                "iauPnm00a", "12", ref status);
            vvd(rbpn[0][2], 0.3639684306407150645e-3, 1e-14,
                "iauPnm00a", "13", ref status);

            vvd(rbpn[1][0], -0.8372535226570394543e-3, 1e-14,
                "iauPnm00a", "21", ref status);
            vvd(rbpn[1][1], 0.9999996486491582471, 1e-12,
                "iauPnm00a", "22", ref status);
            vvd(rbpn[1][2], 0.4132915262664072381e-4, 1e-14,
                "iauPnm00a", "23", ref status);

            vvd(rbpn[2][0], -0.3639337004054317729e-3, 1e-14,
                "iauPnm00a", "31", ref status);
            vvd(rbpn[2][1], -0.4163386925461775873e-4, 1e-14,
                "iauPnm00a", "32", ref status);
            vvd(rbpn[2][2], 0.9999999329094390695, 1e-12,
                "iauPnm00a", "33", ref status);

        }

        static void t_pnm00b(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p n m 0 0 b
        **  - - - - - - - - -
        **
        **  Test iauPnm00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPnm00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rbpn;


            Sofa.iauPnm00b(2400000.5, 50123.9999, out rbpn);

            vvd(rbpn[0][0], 0.9999995832776208280, 1e-12,
                "iauPnm00b", "11", ref status);
            vvd(rbpn[0][1], 0.8372401264429654837e-3, 1e-14,
                "iauPnm00b", "12", ref status);
            vvd(rbpn[0][2], 0.3639691681450271771e-3, 1e-14,
                "iauPnm00b", "13", ref status);

            vvd(rbpn[1][0], -0.8372552234147137424e-3, 1e-14,
                "iauPnm00b", "21", ref status);
            vvd(rbpn[1][1], 0.9999996486477686123, 1e-12,
                "iauPnm00b", "22", ref status);
            vvd(rbpn[1][2], 0.4132832190946052890e-4, 1e-14,
                "iauPnm00b", "23", ref status);

            vvd(rbpn[2][0], -0.3639344385341866407e-3, 1e-14,
                "iauPnm00b", "31", ref status);
            vvd(rbpn[2][1], -0.4163303977421522785e-4, 1e-14,
                "iauPnm00b", "32", ref status);
            vvd(rbpn[2][2], 0.9999999329092049734, 1e-12,
                "iauPnm00b", "33", ref status);

        }

        static void t_pnm06a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p n m 0 6 a
        **  - - - - - - - - -
        **
        **  Test iauPnm06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPnm06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rbpn;


            Sofa.iauPnm06a(2400000.5, 50123.9999, out rbpn);

            vvd(rbpn[0][0], 0.9999995832794205484, 1e-12,
                "iauPnm06a", "11", ref status);
            vvd(rbpn[0][1], 0.8372382772630962111e-3, 1e-14,
                "iauPnm06a", "12", ref status);
            vvd(rbpn[0][2], 0.3639684771140623099e-3, 1e-14,
                "iauPnm06a", "13", ref status);

            vvd(rbpn[1][0], -0.8372533744743683605e-3, 1e-14,
                "iauPnm06a", "21", ref status);
            vvd(rbpn[1][1], 0.9999996486492861646, 1e-12,
                "iauPnm06a", "22", ref status);
            vvd(rbpn[1][2], 0.4132905944611019498e-4, 1e-14,
                "iauPnm06a", "23", ref status);

            vvd(rbpn[2][0], -0.3639337469629464969e-3, 1e-14,
                "iauPnm06a", "31", ref status);
            vvd(rbpn[2][1], -0.4163377605910663999e-4, 1e-14,
                "iauPnm06a", "32", ref status);
            vvd(rbpn[2][2], 0.9999999329094260057, 1e-12,
                "iauPnm06a", "33", ref status);

        }

        static void t_pnm80(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p n m 8 0
        **  - - - - - - - -
        **
        **  Test iauPnm80 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPnm80, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] rmatpn;


            Sofa.iauPnm80(2400000.5, 50123.9999, out rmatpn);

            vvd(rmatpn[0][0], 0.9999995831934611169, 1e-12,
                "iauPnm80", "11", ref status);
            vvd(rmatpn[0][1], 0.8373654045728124011e-3, 1e-14,
                "iauPnm80", "12", ref status);
            vvd(rmatpn[0][2], 0.3639121916933106191e-3, 1e-14,
                "iauPnm80", "13", ref status);

            vvd(rmatpn[1][0], -0.8373804896118301316e-3, 1e-14,
                "iauPnm80", "21", ref status);
            vvd(rmatpn[1][1], 0.9999996485439674092, 1e-12,
                "iauPnm80", "22", ref status);
            vvd(rmatpn[1][2], 0.4130202510421549752e-4, 1e-14,
                "iauPnm80", "23", ref status);

            vvd(rmatpn[2][0], -0.3638774789072144473e-3, 1e-14,
                "iauPnm80", "31", ref status);
            vvd(rmatpn[2][1], -0.4160674085851722359e-4, 1e-14,
                "iauPnm80", "32", ref status);
            vvd(rmatpn[2][2], 0.9999999329310274805, 1e-12,
                "iauPnm80", "33", ref status);

        }

        static void t_pom00(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p o m 0 0
        **  - - - - - - - -
        **
        **  Test iauPom00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPom00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double xp, yp, sp;
            double[][] rpom;


            xp = 2.55060238e-7;
            yp = 1.860359247e-6;
            sp = -0.1367174580728891460e-10;

            Sofa.iauPom00(xp, yp, sp, out rpom);

            vvd(rpom[0][0], 0.9999999999999674721, 1e-12,
                "iauPom00", "11", ref status);
            vvd(rpom[0][1], -0.1367174580728846989e-10, 1e-16,
                "iauPom00", "12", ref status);
            vvd(rpom[0][2], 0.2550602379999972345e-6, 1e-16,
                "iauPom00", "13", ref status);

            vvd(rpom[1][0], 0.1414624947957029801e-10, 1e-16,
                "iauPom00", "21", ref status);
            vvd(rpom[1][1], 0.9999999999982695317, 1e-12,
                "iauPom00", "22", ref status);
            vvd(rpom[1][2], -0.1860359246998866389e-5, 1e-16,
                "iauPom00", "23", ref status);

            vvd(rpom[2][0], -0.2550602379741215021e-6, 1e-16,
                "iauPom00", "31", ref status);
            vvd(rpom[2][1], 0.1860359247002414021e-5, 1e-16,
                "iauPom00", "32", ref status);
            vvd(rpom[2][2], 0.9999999999982370039, 1e-12,
                "iauPom00", "33", ref status);

        }

        static void t_ppp(ref int status)
        /*
        **  - - - - - -
        **   t _ p p p
        **  - - - - - -
        **
        **  Test iauPpp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPpp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] a = new double[3];
            double[] b = new double[3];
            double[] apb = new double[3];


            a[0] = 2.0;
            a[1] = 2.0;
            a[2] = 3.0;

            b[0] = 1.0;
            b[1] = 3.0;
            b[2] = 4.0;

            Sofa.iauPpp(a, b, ref apb);

            vvd(apb[0], 3.0, 1e-12, "iauPpp", "0", ref status);
            vvd(apb[1], 5.0, 1e-12, "iauPpp", "1", ref status);
            vvd(apb[2], 7.0, 1e-12, "iauPpp", "2", ref status);

        }

        static void t_ppsp(ref int status)
        /*
        **  - - - - - - -
        **   t _ p p s p
        **  - - - - - - -
        **
        **  Test iauPpsp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPpsp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s;
            double[] a = new double[3];
            double[] b = new double[3];
            double[] apsb = new double[3];


            a[0] = 2.0;
            a[1] = 2.0;
            a[2] = 3.0;

            s = 5.0;

            b[0] = 1.0;
            b[1] = 3.0;
            b[2] = 4.0;

            Sofa.iauPpsp(a, s, b, ref apsb);

            vvd(apsb[0], 7.0, 1e-12, "iauPpsp", "0", ref status);
            vvd(apsb[1], 17.0, 1e-12, "iauPpsp", "1", ref status);
            vvd(apsb[2], 23.0, 1e-12, "iauPpsp", "2", ref status);

        }

        static void t_pr00(ref int status)
        /*
        **  - - - - - - -
        **   t _ p r 0 0
        **  - - - - - - -
        **
        **  Test iauPr00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPr00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double dpsipr, depspr;

            Sofa.iauPr00(2400000.5, 53736, out dpsipr, out depspr);

            vvd(dpsipr, -0.8716465172668347629e-7, 1e-22,
               "iauPr00", "dpsipr", ref status);
            vvd(depspr, -0.7342018386722813087e-8, 1e-22,
               "iauPr00", "depspr", ref status);

        }

        static void t_prec76(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p r e c 7 6
        **  - - - - - - - - -
        **
        **  Test iauPrec76 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPrec76, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double ep01, ep02, ep11, ep12, zeta, z, theta;


            ep01 = 2400000.5;
            ep02 = 33282.0;
            ep11 = 2400000.5;
            ep12 = 51544.0;

            Sofa.iauPrec76(ep01, ep02, ep11, ep12, out zeta, out z, out theta);

            vvd(zeta, 0.5588961642000161243e-2, 1e-12,
                "iauPrec76", "zeta", ref status);
            vvd(z, 0.5589922365870680624e-2, 1e-12,
                "iauPrec76", "z", ref status);
            vvd(theta, 0.4858945471687296760e-2, 1e-12,
                "iauPrec76", "theta", ref status);

        }

        static void t_pv2p(ref int status)
        /*
        **  - - - - - - -
        **   t _ p v 2 p
        **  - - - - - - -
        **
        **  Test iauPv2p function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPv2p, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p = new double[3];
            double[][] pv = { new double[3], new double[3] };

            pv[0][0] = 0.3;
            pv[0][1] = 1.2;
            pv[0][2] = -2.5;

            pv[1][0] = -0.5;
            pv[1][1] = 3.1;
            pv[1][2] = 0.9;

            Sofa.iauPv2p(pv, out p);

            vvd(p[0], 0.3, 0.0, "iauPv2p", "1", ref status);
            vvd(p[1], 1.2, 0.0, "iauPv2p", "2", ref status);
            vvd(p[2], -2.5, 0.0, "iauPv2p", "3", ref status);

        }

        static void t_pv2s(ref int status)
        /*
        **  - - - - - - -
        **   t _ p v 2 s
        **  - - - - - - -
        **
        **  Test iauPv2s function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPv2s, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta, phi, r, td, pd, rd;
            double[][] pv = { new double[3], new double[3] };


            pv[0][0] = -0.4514964673880165;
            pv[0][1] = 0.03093394277342585;
            pv[0][2] = 0.05594668105108779;

            pv[1][0] = 1.292270850663260e-5;
            pv[1][1] = 2.652814182060692e-6;
            pv[1][2] = 2.568431853930293e-6;

            Sofa.iauPv2s(pv, out theta, out phi, out r, out td, out pd, out rd);

            vvd(theta, 3.073185307179586515, 1e-12, "iauPv2s", "theta", ref status);
            vvd(phi, 0.1229999999999999992, 1e-12, "iauPv2s", "phi", ref status);
            vvd(r, 0.4559999999999999757, 1e-12, "iauPv2s", "r", ref status);
            vvd(td, -0.7800000000000000364e-5, 1e-16, "iauPv2s", "td", ref status);
            vvd(pd, 0.9010000000000001639e-5, 1e-16, "iauPv2s", "pd", ref status);
            vvd(rd, -0.1229999999999999832e-4, 1e-16, "iauPv2s", "rd", ref status);

        }

        static void t_pvdpv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p v d p v
        **  - - - - - - - -
        **
        **  Test iauPvdpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvdpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] adb = new double[3];
            double[][] a = { new double[3], new double[3] };
            double[][] b = { new double[3], new double[3] };


            a[0][0] = 2.0;
            a[0][1] = 2.0;
            a[0][2] = 3.0;

            a[1][0] = 6.0;
            a[1][1] = 0.0;
            a[1][2] = 4.0;

            b[0][0] = 1.0;
            b[0][1] = 3.0;
            b[0][2] = 4.0;

            b[1][0] = 0.0;
            b[1][1] = 2.0;
            b[1][2] = 8.0;

            Sofa.iauPvdpv(a, b, out adb);

            vvd(adb[0], 20.0, 1e-12, "iauPvdpv", "1", ref status);
            vvd(adb[1], 50.0, 1e-12, "iauPvdpv", "2", ref status);

        }

        static void t_pvm(ref int status)
        /*
        **  - - - - - -
        **   t _ p v m
        **  - - - - - -
        **
        **  Test iauPvm function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvm, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double r, s;
            double[][] pv = { new double[3], new double[3] };


            pv[0][0] = 0.3;
            pv[0][1] = 1.2;
            pv[0][2] = -2.5;

            pv[1][0] = 0.45;
            pv[1][1] = -0.25;
            pv[1][2] = 1.1;

            Sofa.iauPvm(pv, out r, out s);

            vvd(r, 2.789265136196270604, 1e-12, "iauPvm", "r", ref status);
            vvd(s, 1.214495780149111922, 1e-12, "iauPvm", "s", ref status);

        }

        static void t_pvmpv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p v m p v
        **  - - - - - - - -
        **
        **  Test iauPvmpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvmpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] a = { new double[3], new double[3] };
            double[][] b = { new double[3], new double[3] };
            double[][] amb = { new double[3], new double[3] };


            a[0][0] = 2.0;
            a[0][1] = 2.0;
            a[0][2] = 3.0;

            a[1][0] = 5.0;
            a[1][1] = 6.0;
            a[1][2] = 3.0;

            b[0][0] = 1.0;
            b[0][1] = 3.0;
            b[0][2] = 4.0;

            b[1][0] = 3.0;
            b[1][1] = 2.0;
            b[1][2] = 1.0;

            Sofa.iauPvmpv(a, b, ref amb);

            vvd(amb[0][0], 1.0, 1e-12, "iauPvmpv", "11", ref status);
            vvd(amb[0][1], -1.0, 1e-12, "iauPvmpv", "21", ref status);
            vvd(amb[0][2], -1.0, 1e-12, "iauPvmpv", "31", ref status);

            vvd(amb[1][0], 2.0, 1e-12, "iauPvmpv", "12", ref status);
            vvd(amb[1][1], 4.0, 1e-12, "iauPvmpv", "22", ref status);
            vvd(amb[1][2], 2.0, 1e-12, "iauPvmpv", "32", ref status);

        }

        static void t_pvppv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p v p p v
        **  - - - - - - - -
        **
        **  Test iauPvppv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvppv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] a = { new double[3], new double[3] };
            double[][] b = { new double[3], new double[3] };
            double[][] apb = { new double[3], new double[3] };


            a[0][0] = 2.0;
            a[0][1] = 2.0;
            a[0][2] = 3.0;

            a[1][0] = 5.0;
            a[1][1] = 6.0;
            a[1][2] = 3.0;

            b[0][0] = 1.0;
            b[0][1] = 3.0;
            b[0][2] = 4.0;

            b[1][0] = 3.0;
            b[1][1] = 2.0;
            b[1][2] = 1.0;

            Sofa.iauPvppv(a, b, ref apb);

            vvd(apb[0][0], 3.0, 1e-12, "iauPvppv", "p1", ref status);
            vvd(apb[0][1], 5.0, 1e-12, "iauPvppv", "p2", ref status);
            vvd(apb[0][2], 7.0, 1e-12, "iauPvppv", "p3", ref status);

            vvd(apb[1][0], 8.0, 1e-12, "iauPvppv", "v1", ref status);
            vvd(apb[1][1], 8.0, 1e-12, "iauPvppv", "v2", ref status);
            vvd(apb[1][2], 4.0, 1e-12, "iauPvppv", "v3", ref status);

        }

        static void t_pvstar(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ p v s t a r
        **  - - - - - - - - -
        **
        **  Test iauPvstar function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvstar, vvd, viv
        **
        **  This revision:  2017 March 15
        */
        {
            double ra, dec, pmr, pmd, px, rv;
            double[][] pv = { new double[3], new double[3] };
            int j;


            pv[0][0] = 126668.5912743160601;
            pv[0][1] = 2136.792716839935195;
            pv[0][2] = -245251.2339876830091;

            pv[1][0] = -0.4051854035740712739e-2;
            pv[1][1] = -0.6253919754866173866e-2;
            pv[1][2] = 0.1189353719774107189e-1;

            j = Sofa.iauPvstar(pv, out ra, out dec, out pmr, out pmd, out px, out rv);

            vvd(ra, 0.1686756e-1, 1e-12, "iauPvstar", "ra", ref status);
            vvd(dec, -1.093989828, 1e-12, "iauPvstar", "dec", ref status);
            vvd(pmr, -0.1783235160000472788e-4, 1e-16, "iauPvstar", "pmr", ref status);
            vvd(pmd, 0.2336024047000619347e-5, 1e-16, "iauPvstar", "pmd", ref status);
            vvd(px, 0.74723, 1e-12, "iauPvstar", "px", ref status);
            vvd(rv, -21.60000010107306010, 1e-11, "iauPvstar", "rv", ref status);

            viv(j, 0, "iauPvstar", "j", ref status);

        }

        static void t_pvtob(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p v t o b
        **  - - - - - - - -
        **
        **  Test iauPvtob function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvtob, vvd
        **
        **  This revision:  2013 October 2
        */
        {
            double elong, phi, hm, xp, yp, sp, theta;
            double[][] pv;


            elong = 2.0;
            phi = 0.5;
            hm = 3000.0;
            xp = 1e-6;
            yp = -0.5e-6;
            sp = 1e-8;
            theta = 5.0;

            Sofa.iauPvtob(elong, phi, hm, xp, yp, sp, theta, out pv);

            vvd(pv[0][0], 4225081.367071159207, 1e-5,
                          "iauPvtob", "p(1)", ref status);
            vvd(pv[0][1], 3681943.215856198144, 1e-5,
                          "iauPvtob", "p(2)", ref status);
            vvd(pv[0][2], 3041149.399241260785, 1e-5,
                          "iauPvtob", "p(3)", ref status);
            vvd(pv[1][0], -268.4915389365998787, 1e-9,
                          "iauPvtob", "v(1)", ref status);
            vvd(pv[1][1], 308.0977983288903123, 1e-9,
                          "iauPvtob", "v(2)", ref status);
            vvd(pv[1][2], 0, 0,
                          "iauPvtob", "v(3)", ref status);

        }

        static void t_pvu(ref int status)
        /*
        **  - - - - - -
        **   t _ p v u
        **  - - - - - -
        **
        **  Test iauPvu function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvu, vvd
        **
        **  This revision:  2021 January 5
        */
        {
            double[][] pv = { new double[3], new double[3] };
            double[][] upv = { new double[3], new double[3] };


            pv[0][0] = 126668.5912743160734;
            pv[0][1] = 2136.792716839935565;
            pv[0][2] = -245251.2339876830229;

            pv[1][0] = -0.4051854035740713039e-2;
            pv[1][1] = -0.6253919754866175788e-2;
            pv[1][2] = 0.1189353719774107615e-1;

            Sofa.iauPvu(2920.0, pv, ref upv);

            vvd(upv[0][0], 126656.7598605317105, 1e-6,
                "iauPvu", "p1", ref status);
            vvd(upv[0][1], 2118.531271155726332, 1e-8,
                "iauPvu", "p2", ref status);
            vvd(upv[0][2], -245216.5048590656190, 1e-6,
                "iauPvu", "p3", ref status);

            vvd(upv[1][0], -0.4051854035740713039e-2, 1e-12,
                "iauPvu", "v1", ref status);
            vvd(upv[1][1], -0.6253919754866175788e-2, 1e-12,
                "iauPvu", "v2", ref status);
            vvd(upv[1][2], 0.1189353719774107615e-1, 1e-12,
                "iauPvu", "v3", ref status);

        }

        static void t_pvup(ref int status)
        /*
        **  - - - - - - -
        **   t _ p v u p
        **  - - - - - - -
        **
        **  Test iauPvup function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvup, vvd
        **
        **  This revision:  2021 January 5
        */
        {
            double[] p;
            double[][] pv = { new double[3], new double[3] };


            pv[0][0] = 126668.5912743160734;
            pv[0][1] = 2136.792716839935565;
            pv[0][2] = -245251.2339876830229;

            pv[1][0] = -0.4051854035740713039e-2;
            pv[1][1] = -0.6253919754866175788e-2;
            pv[1][2] = 0.1189353719774107615e-1;

            Sofa.iauPvup(2920.0, pv, out p);

            vvd(p[0], 126656.7598605317105, 1e-6, "iauPvup", "1", ref status);
            vvd(p[1], 2118.531271155726332, 1e-8, "iauPvup", "2", ref status);
            vvd(p[2], -245216.5048590656190, 1e-6, "iauPvup", "3", ref status);

        }

        static void t_pvxpv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ p v x p v
        **  - - - - - - - -
        **
        **  Test iauPvxpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPvxpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] a = { new double[3], new double[3] };
            double[][] b = { new double[3], new double[3] };
            double[][] axb = { new double[3], new double[3] };

            a[0][0] = 2.0;
            a[0][1] = 2.0;
            a[0][2] = 3.0;

            a[1][0] = 6.0;
            a[1][1] = 0.0;
            a[1][2] = 4.0;

            b[0][0] = 1.0;
            b[0][1] = 3.0;
            b[0][2] = 4.0;

            b[1][0] = 0.0;
            b[1][1] = 2.0;
            b[1][2] = 8.0;

            Sofa.iauPvxpv(a, b, ref axb);

            vvd(axb[0][0], -1.0, 1e-12, "iauPvxpv", "p1", ref status);
            vvd(axb[0][1], -5.0, 1e-12, "iauPvxpv", "p2", ref status);
            vvd(axb[0][2], 4.0, 1e-12, "iauPvxpv", "p3", ref status);

            vvd(axb[1][0], -2.0, 1e-12, "iauPvxpv", "v1", ref status);
            vvd(axb[1][1], -36.0, 1e-12, "iauPvxpv", "v2", ref status);
            vvd(axb[1][2], 22.0, 1e-12, "iauPvxpv", "v3", ref status);

        }

        static void t_pxp(ref int status)
        /*
        **  - - - - - -
        **   t _ p x p
        **  - - - - - -
        **
        **  Test iauPxp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauPxp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] a = new double[3];
            double[] b = new double[3];
            double[] axb = new double[3];

            a[0] = 2.0;
            a[1] = 2.0;
            a[2] = 3.0;

            b[0] = 1.0;
            b[1] = 3.0;
            b[2] = 4.0;

            Sofa.iauPxp(a, b, ref axb);

            vvd(axb[0], -1.0, 1e-12, "iauPxp", "1", ref status);
            vvd(axb[1], -5.0, 1e-12, "iauPxp", "2", ref status);
            vvd(axb[2], 4.0, 1e-12, "iauPxp", "3", ref status);

        }

        static void t_refco(ref int status)
        /*
        **  - - - - - - - -
        **   t _ r e f c o
        **  - - - - - - - -
        **
        **  Test iauRefco function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRefco, vvd
        **
        **  This revision:  2013 October 2
        */
        {
            double phpa, tc, rh, wl, refa, refb;


            phpa = 800.0;
            tc = 10.0;
            rh = 0.9;
            wl = 0.4;

            Sofa.iauRefco(phpa, tc, rh, wl, out refa, out refb);

            vvd(refa, 0.2264949956241415009e-3, 1e-15,
                      "iauRefco", "refa", ref status);
            vvd(refb, -0.2598658261729343970e-6, 1e-18,
                      "iauRefco", "refb", ref status);

        }

        static void t_rm2v(ref int status)
        /*
        **  - - - - - - -
        **   t _ r m 2 v
        **  - - - - - - -
        **
        **  Test iauRm2v function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRm2v, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] w;
            double[][] r = { new double[3], new double[3], new double[3] };


            r[0][0] = 0.00;
            r[0][1] = -0.80;
            r[0][2] = -0.60;

            r[1][0] = 0.80;
            r[1][1] = -0.36;
            r[1][2] = 0.48;

            r[2][0] = 0.60;
            r[2][1] = 0.48;
            r[2][2] = -0.64;

            Sofa.iauRm2v(r, out w);

            vvd(w[0], 0.0, 1e-12, "iauRm2v", "1", ref status);
            vvd(w[1], 1.413716694115406957, 1e-12, "iauRm2v", "2", ref status);
            vvd(w[2], -1.884955592153875943, 1e-12, "iauRm2v", "3", ref status);

        }

        static void t_rv2m(ref int status)
        /*
        **  - - - - - - -
        **   t _ r v 2 m
        **  - - - - - - -
        **
        **  Test iauRv2m function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRv2m, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] w = new double[3];
            double[][] r;


            w[0] = 0.0;
            w[1] = 1.41371669;
            w[2] = -1.88495559;

            Sofa.iauRv2m(w, out r);

            vvd(r[0][0], -0.7071067782221119905, 1e-14, "iauRv2m", "11", ref status);
            vvd(r[0][1], -0.5656854276809129651, 1e-14, "iauRv2m", "12", ref status);
            vvd(r[0][2], -0.4242640700104211225, 1e-14, "iauRv2m", "13", ref status);

            vvd(r[1][0], 0.5656854276809129651, 1e-14, "iauRv2m", "21", ref status);
            vvd(r[1][1], -0.0925483394532274246, 1e-14, "iauRv2m", "22", ref status);
            vvd(r[1][2], -0.8194112531408833269, 1e-14, "iauRv2m", "23", ref status);

            vvd(r[2][0], 0.4242640700104211225, 1e-14, "iauRv2m", "31", ref status);
            vvd(r[2][1], -0.8194112531408833269, 1e-14, "iauRv2m", "32", ref status);
            vvd(r[2][2], 0.3854415612311154341, 1e-14, "iauRv2m", "33", ref status);

        }

        static void t_rx(ref int status)
        /*
        **  - - - - -
        **   t _ r x
        **  - - - - -
        **
        **  Test iauRx function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRx, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double phi;
            double[][] r = { new double[3], new double[3], new double[3] };


            phi = 0.3456789;

            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            Sofa.iauRx(phi, r);

            vvd(r[0][0], 2.0, 0.0, "iauRx", "11", ref status);
            vvd(r[0][1], 3.0, 0.0, "iauRx", "12", ref status);
            vvd(r[0][2], 2.0, 0.0, "iauRx", "13", ref status);

            vvd(r[1][0], 3.839043388235612460, 1e-12, "iauRx", "21", ref status);
            vvd(r[1][1], 3.237033249594111899, 1e-12, "iauRx", "22", ref status);
            vvd(r[1][2], 4.516714379005982719, 1e-12, "iauRx", "23", ref status);

            vvd(r[2][0], 1.806030415924501684, 1e-12, "iauRx", "31", ref status);
            vvd(r[2][1], 3.085711545336372503, 1e-12, "iauRx", "32", ref status);
            vvd(r[2][2], 3.687721683977873065, 1e-12, "iauRx", "33", ref status);

        }

        static void t_rxp(ref int status)
        /*
        **  - - - - - -
        **   t _ r x p
        **  - - - - - -
        **
        **  Test iauRxp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRxp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r = { new double[3], new double[3], new double[3] };
            double[] p = new double[3];
            double[] rp = new double[3];


            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            p[0] = 0.2;
            p[1] = 1.5;
            p[2] = 0.1;

            Sofa.iauRxp(r, p, ref rp);

            vvd(rp[0], 5.1, 1e-12, "iauRxp", "1", ref status);
            vvd(rp[1], 3.9, 1e-12, "iauRxp", "2", ref status);
            vvd(rp[2], 7.1, 1e-12, "iauRxp", "3", ref status);

        }

        static void t_rxpv(ref int status)
        /*
        **  - - - - - - -
        **   t _ r x p v
        **  - - - - - - -
        **
        **  Test iauRxpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRxpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r = { new double[3], new double[3], new double[3] };
            double[][] pv = { new double[3], new double[3] };
            double[][] rpv = { new double[3], new double[3] };


            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            pv[0][0] = 0.2;
            pv[0][1] = 1.5;
            pv[0][2] = 0.1;

            pv[1][0] = 1.5;
            pv[1][1] = 0.2;
            pv[1][2] = 0.1;

            Sofa.iauRxpv(r, pv, ref rpv);

            vvd(rpv[0][0], 5.1, 1e-12, "iauRxpv", "11", ref status);
            vvd(rpv[1][0], 3.8, 1e-12, "iauRxpv", "12", ref status);

            vvd(rpv[0][1], 3.9, 1e-12, "iauRxpv", "21", ref status);
            vvd(rpv[1][1], 5.2, 1e-12, "iauRxpv", "22", ref status);

            vvd(rpv[0][2], 7.1, 1e-12, "iauRxpv", "31", ref status);
            vvd(rpv[1][2], 5.8, 1e-12, "iauRxpv", "32", ref status);

        }

        static void t_rxr(ref int status)
        /*
        **  - - - - - -
        **   t _ r x r
        **  - - - - - -
        **
        **  Test iauRxr function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRxr, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] a = { new double[3], new double[3], new double[3] };
            double[][] b = { new double[3], new double[3], new double[3] };
            double[][] atb = { new double[3], new double[3], new double[3] };


            a[0][0] = 2.0;
            a[0][1] = 3.0;
            a[0][2] = 2.0;

            a[1][0] = 3.0;
            a[1][1] = 2.0;
            a[1][2] = 3.0;

            a[2][0] = 3.0;
            a[2][1] = 4.0;
            a[2][2] = 5.0;

            b[0][0] = 1.0;
            b[0][1] = 2.0;
            b[0][2] = 2.0;

            b[1][0] = 4.0;
            b[1][1] = 1.0;
            b[1][2] = 1.0;

            b[2][0] = 3.0;
            b[2][1] = 0.0;
            b[2][2] = 1.0;

            Sofa.iauRxr(a, b, ref atb);

            vvd(atb[0][0], 20.0, 1e-12, "iauRxr", "11", ref status);
            vvd(atb[0][1], 7.0, 1e-12, "iauRxr", "12", ref status);
            vvd(atb[0][2], 9.0, 1e-12, "iauRxr", "13", ref status);

            vvd(atb[1][0], 20.0, 1e-12, "iauRxr", "21", ref status);
            vvd(atb[1][1], 8.0, 1e-12, "iauRxr", "22", ref status);
            vvd(atb[1][2], 11.0, 1e-12, "iauRxr", "23", ref status);

            vvd(atb[2][0], 34.0, 1e-12, "iauRxr", "31", ref status);
            vvd(atb[2][1], 10.0, 1e-12, "iauRxr", "32", ref status);
            vvd(atb[2][2], 15.0, 1e-12, "iauRxr", "33", ref status);

        }

        static void t_ry(ref int status)
        /*
        **  - - - - -
        **   t _ r y
        **  - - - - -
        **
        **  Test iauRy function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRy, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double theta;
            double[][] r = { new double[3], new double[3], new double[3] };


            theta = 0.3456789;

            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            Sofa.iauRy(theta, r);

            vvd(r[0][0], 0.8651847818978159930, 1e-12, "iauRy", "11", ref status);
            vvd(r[0][1], 1.467194920539316554, 1e-12, "iauRy", "12", ref status);
            vvd(r[0][2], 0.1875137911274457342, 1e-12, "iauRy", "13", ref status);

            vvd(r[1][0], 3, 1e-12, "iauRy", "21", ref status);
            vvd(r[1][1], 2, 1e-12, "iauRy", "22", ref status);
            vvd(r[1][2], 3, 1e-12, "iauRy", "23", ref status);

            vvd(r[2][0], 3.500207892850427330, 1e-12, "iauRy", "31", ref status);
            vvd(r[2][1], 4.779889022262298150, 1e-12, "iauRy", "32", ref status);
            vvd(r[2][2], 5.381899160903798712, 1e-12, "iauRy", "33", ref status);

        }

        static void t_rz(ref int status)
        /*
        **  - - - - -
        **   t _ r z
        **  - - - - -
        **
        **  Test iauRz function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauRz, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double psi;
            double[][] r = { new double[3], new double[3], new double[3] };


            psi = 0.3456789;

            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            Sofa.iauRz(psi, r);

            vvd(r[0][0], 2.898197754208926769, 1e-12, "iauRz", "11", ref status);
            vvd(r[0][1], 3.500207892850427330, 1e-12, "iauRz", "12", ref status);
            vvd(r[0][2], 2.898197754208926769, 1e-12, "iauRz", "13", ref status);

            vvd(r[1][0], 2.144865911309686813, 1e-12, "iauRz", "21", ref status);
            vvd(r[1][1], 0.865184781897815993, 1e-12, "iauRz", "22", ref status);
            vvd(r[1][2], 2.144865911309686813, 1e-12, "iauRz", "23", ref status);

            vvd(r[2][0], 3.0, 1e-12, "iauRz", "31", ref status);
            vvd(r[2][1], 4.0, 1e-12, "iauRz", "32", ref status);
            vvd(r[2][2], 5.0, 1e-12, "iauRz", "33", ref status);

        }

        static void t_s00a(ref int status)
        /*
        **  - - - - - - -
        **   t _ s 0 0 a
        **  - - - - - - -
        **
        **  Test iauS00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s;


            s = Sofa.iauS00a(2400000.5, 52541.0);

            vvd(s, -0.1340684448919163584e-7, 1e-18, "iauS00a", "", ref status);

        }

        static void t_s00b(ref int status)
        /*
        **  - - - - - - -
        **   t _ s 0 0 b
        **  - - - - - - -
        **
        **  Test iauS00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s;


            s = Sofa.iauS00b(2400000.5, 52541.0);

            vvd(s, -0.1340695782951026584e-7, 1e-18, "iauS00b", "", ref status);

        }

        static void t_s00(ref int status)
        /*
        **  - - - - - -
        **   t _ s 0 0
        **  - - - - - -
        **
        **  Test iauS00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y, s;


            x = 0.5791308486706011000e-3;
            y = 0.4020579816732961219e-4;

            s = Sofa.iauS00(2400000.5, 53736.0, x, y);

            vvd(s, -0.1220036263270905693e-7, 1e-18, "iauS00", "", ref status);

        }

        static void t_s06a(ref int status)
        /*
        **  - - - - - - -
        **   t _ s 0 6 a
        **  - - - - - - -
        **
        **  Test iauS06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s;


            s = Sofa.iauS06a(2400000.5, 52541.0);

            vvd(s, -0.1340680437291812383e-7, 1e-18, "iauS06a", "", ref status);

        }

        static void t_s06(ref int status)
        /*
        **  - - - - - -
        **   t _ s 0 6
        **  - - - - - -
        **
        **  Test iauS06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y, s;


            x = 0.5791308486706011000e-3;
            y = 0.4020579816732961219e-4;

            s = Sofa.iauS06(2400000.5, 53736.0, x, y);

            vvd(s, -0.1220032213076463117e-7, 1e-18, "iauS06", "", ref status);

        }

        static void t_s2c(ref int status)
        /*
        **  - - - - - -
        **   t _ s 2 c
        **  - - - - - -
        **
        **  Test iauS2c function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS2c, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] c;


            Sofa.iauS2c(3.0123, -0.999, out c);

            vvd(c[0], -0.5366267667260523906, 1e-12, "iauS2c", "1", ref status);
            vvd(c[1], 0.0697711109765145365, 1e-12, "iauS2c", "2", ref status);
            vvd(c[2], -0.8409302618566214041, 1e-12, "iauS2c", "3", ref status);

        }

        static void t_s2p(ref int status)
        /*
        **  - - - - - -
        **   t _ s 2 p
        **  - - - - - -
        **
        **  Test iauS2p function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS2p, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p;


            Sofa.iauS2p(-3.21, 0.123, 0.456, out p);

            vvd(p[0], -0.4514964673880165228, 1e-12, "iauS2p", "x", ref status);
            vvd(p[1], 0.0309339427734258688, 1e-12, "iauS2p", "y", ref status);
            vvd(p[2], 0.0559466810510877933, 1e-12, "iauS2p", "z", ref status);

        }

        static void t_s2pv(ref int status)
        /*
        **  - - - - - - -
        **   t _ s 2 p v
        **  - - - - - - -
        **
        **  Test iauS2pv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS2pv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] pv;


            Sofa.iauS2pv(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5, out pv);

            vvd(pv[0][0], -0.4514964673880165228, 1e-12, "iauS2pv", "x", ref status);
            vvd(pv[0][1], 0.0309339427734258688, 1e-12, "iauS2pv", "y", ref status);
            vvd(pv[0][2], 0.0559466810510877933, 1e-12, "iauS2pv", "z", ref status);

            vvd(pv[1][0], 0.1292270850663260170e-4, 1e-16,
                "iauS2pv", "vx", ref status);
            vvd(pv[1][1], 0.2652814182060691422e-5, 1e-16,
                "iauS2pv", "vy", ref status);
            vvd(pv[1][2], 0.2568431853930292259e-5, 1e-16,
                "iauS2pv", "vz", ref status);

        }

        static void t_s2xpv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ s 2 x p v
        **  - - - - - - - -
        **
        **  Test iauS2xpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauS2xpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s1, s2;
            double[][] pv = { new double[3], new double[3] };
            double[][] spv = { new double[3], new double[3] };


            s1 = 2.0;
            s2 = 3.0;

            pv[0][0] = 0.3;
            pv[0][1] = 1.2;
            pv[0][2] = -2.5;

            pv[1][0] = 0.5;
            pv[1][1] = 2.3;
            pv[1][2] = -0.4;

            Sofa.iauS2xpv(s1, s2, pv, ref spv);

            vvd(spv[0][0], 0.6, 1e-12, "iauS2xpv", "p1", ref status);
            vvd(spv[0][1], 2.4, 1e-12, "iauS2xpv", "p2", ref status);
            vvd(spv[0][2], -5.0, 1e-12, "iauS2xpv", "p3", ref status);

            vvd(spv[1][0], 1.5, 1e-12, "iauS2xpv", "v1", ref status);
            vvd(spv[1][1], 6.9, 1e-12, "iauS2xpv", "v2", ref status);
            vvd(spv[1][2], -1.2, 1e-12, "iauS2xpv", "v3", ref status);

        }

        static void t_sepp(ref int status)
        /*
        **  - - - - - - -
        **   t _ s e p p
        **  - - - - - - -
        **
        **  Test iauSepp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauSepp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s;
            double[] a = new double[3];
            double[] b = new double[3];


            a[0] = 1.0;
            a[1] = 0.1;
            a[2] = 0.2;

            b[0] = -3.0;
            b[1] = 1e-3;
            b[2] = 0.2;

            s = Sofa.iauSepp(a, b);

            vvd(s, 2.860391919024660768, 1e-12, "iauSepp", "", ref status);

        }

        static void t_seps(ref int status)
        /*
        **  - - - - - - -
        **   t _ s e p s
        **  - - - - - - -
        **
        **  Test iauSeps function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauSeps, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double al, ap, bl, bp, s;


            al = 1.0;
            ap = 0.1;

            bl = 0.2;
            bp = -3.0;

            s = Sofa.iauSeps(al, ap, bl, bp);

            vvd(s, 2.346722016996998842, 1e-14, "iauSeps", "", ref status);

        }

        static void t_sp00(ref int status)
        /*
        **  - - - - - - -
        **   t _ s p 0 0
        **  - - - - - - -
        **
        **  Test iauSp00 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauSp00, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            vvd(Sofa.iauSp00(2400000.5, 52541.0),
                -0.6216698469981019309e-11, 1e-12, "iauSp00", "", ref status);

        }

        static void t_starpm(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ s t a r p m
        **  - - - - - - - - -
        **
        **  Test iauStarpm function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauStarpm, vvd, viv
        **
        **  This revision:  2017 March 15
        */
        {
            double ra1, dec1, pmr1, pmd1, px1, rv1;
            double ra2, dec2, pmr2, pmd2, px2, rv2;
            int j;


            ra1 = 0.01686756;
            dec1 = -1.093989828;
            pmr1 = -1.78323516e-5;
            pmd1 = 2.336024047e-6;
            px1 = 0.74723;
            rv1 = -21.6;

            j = Sofa.iauStarpm(ra1, dec1, pmr1, pmd1, px1, rv1,
                          2400000.5, 50083.0, 2400000.5, 53736.0,
                          out ra2, out dec2, out pmr2, out pmd2, out px2, out rv2);

            vvd(ra2, 0.01668919069414256149, 1e-13,
                "iauStarpm", "ra", ref status);
            vvd(dec2, -1.093966454217127897, 1e-13,
                "iauStarpm", "dec", ref status);
            vvd(pmr2, -0.1783662682153176524e-4, 1e-17,
                "iauStarpm", "pmr", ref status);
            vvd(pmd2, 0.2338092915983989595e-5, 1e-17,
                "iauStarpm", "pmd", ref status);
            vvd(px2, 0.7473533835317719243, 1e-13,
                "iauStarpm", "px", ref status);
            vvd(rv2, -21.59905170476417175, 1e-11,
                "iauStarpm", "rv", ref status);

            viv(j, 0, "iauStarpm", "j", ref status);

        }

        static void t_starpv(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ s t a r p v
        **  - - - - - - - - -
        **
        **  Test iauStarpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauStarpv, vvd, viv
        **
        **  This revision:  2017 March 15
        */
        {
            double ra, dec, pmr, pmd, px, rv;
            double[][] pv = { new double[3], new double[3] };
            int j;


            ra = 0.01686756;
            dec = -1.093989828;
            pmr = -1.78323516e-5;
            pmd = 2.336024047e-6;
            px = 0.74723;
            rv = -21.6;

            j = Sofa.iauStarpv(ra, dec, pmr, pmd, px, rv, out pv);

            vvd(pv[0][0], 126668.5912743160601, 1e-10,
                "iauStarpv", "11", ref status);
            vvd(pv[0][1], 2136.792716839935195, 1e-12,
                "iauStarpv", "12", ref status);
            vvd(pv[0][2], -245251.2339876830091, 1e-10,
                "iauStarpv", "13", ref status);

            vvd(pv[1][0], -0.4051854008955659551e-2, 1e-13,
                "iauStarpv", "21", ref status);
            vvd(pv[1][1], -0.6253919754414777970e-2, 1e-15,
                "iauStarpv", "22", ref status);
            vvd(pv[1][2], 0.1189353714588109341e-1, 1e-13,
                "iauStarpv", "23", ref status);

            viv(j, 0, "iauStarpv", "j", ref status);

        }

        static void t_sxp(ref int status)
        /*
        **  - - - - - -
        **   t _ s x p
        **  - - - - - -
        **
        **  Test iauSxp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauSxp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s;
            double[] p = new double[3];
            double[] sp = new double[3];


            s = 2.0;

            p[0] = 0.3;
            p[1] = 1.2;
            p[2] = -2.5;

            Sofa.iauSxp(s, p, ref sp);

            vvd(sp[0], 0.6, 0.0, "iauSxp", "1", ref status);
            vvd(sp[1], 2.4, 0.0, "iauSxp", "2", ref status);
            vvd(sp[2], -5.0, 0.0, "iauSxp", "3", ref status);

        }


        static void t_sxpv(ref int status)
        /*
        **  - - - - - - -
        **   t _ s x p v
        **  - - - - - - -
        **
        **  Test iauSxpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauSxpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double s;
            double[][] pv = { new double[3], new double[3] };
            double[][] spv = { new double[3], new double[3] };


            s = 2.0;

            pv[0][0] = 0.3;
            pv[0][1] = 1.2;
            pv[0][2] = -2.5;

            pv[1][0] = 0.5;
            pv[1][1] = 3.2;
            pv[1][2] = -0.7;

            Sofa.iauSxpv(s, pv, ref spv);

            vvd(spv[0][0], 0.6, 0.0, "iauSxpv", "p1", ref status);
            vvd(spv[0][1], 2.4, 0.0, "iauSxpv", "p2", ref status);
            vvd(spv[0][2], -5.0, 0.0, "iauSxpv", "p3", ref status);

            vvd(spv[1][0], 1.0, 0.0, "iauSxpv", "v1", ref status);
            vvd(spv[1][1], 6.4, 0.0, "iauSxpv", "v2", ref status);
            vvd(spv[1][2], -1.4, 0.0, "iauSxpv", "v3", ref status);

        }

        static void t_taitt(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t a i t t
        **  - - - - - - - -
        **
        **  Test iauTaitt function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTaitt, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double t1, t2;
            int j;


            j = Sofa.iauTaitt(2453750.5, 0.892482639, out t1, out t2);

            vvd(t1, 2453750.5, 1e-6, "iauTaitt", "t1", ref status);
            vvd(t2, 0.892855139, 1e-12, "iauTaitt", "t2", ref status);
            viv(j, 0, "iauTaitt", "j", ref status);

        }

        static void t_taiut1(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ t a i u t 1
        **  - - - - - - - - -
        **
        **  Test iauTaiut1 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTaiut1, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double u1, u2;
            int j;


            j = Sofa.iauTaiut1(2453750.5, 0.892482639, -32.6659, out u1, out u2);

            vvd(u1, 2453750.5, 1e-6, "iauTaiut1", "u1", ref status);
            vvd(u2, 0.8921045614537037037, 1e-12, "iauTaiut1", "u2", ref status);
            viv(j, 0, "iauTaiut1", "j", ref status);

        }

        static void t_taiutc(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ t a i u t c
        **  - - - - - - - - -
        **
        **  Test iauTaiutc function.
        **
        **  Returned:
        **     status    LOGICAL     TRUE = success, FALSE = fail
        **
        **  Called:  iauTaiutc, vvd, viv
        **
        **  This revision:  2013 October 3
        */
        {
            double u1, u2;
            int j;


            j = Sofa.iauTaiutc(2453750.5, 0.892482639, out u1, out u2);

            vvd(u1, 2453750.5, 1e-6, "iauTaiutc", "u1", ref status);
            vvd(u2, 0.8921006945555555556, 1e-12, "iauTaiutc", "u2", ref status);
            viv(j, 0, "iauTaiutc", "j", ref status);

        }

        static void t_tcbtdb(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ t c b t d b
        **  - - - - - - - - -
        **
        **  Test iauTcbtdb function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTcbtdb, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double b1, b2;
            int j;


            j = Sofa.iauTcbtdb(2453750.5, 0.893019599, out b1, out b2);

            vvd(b1, 2453750.5, 1e-6, "iauTcbtdb", "b1", ref status);
            vvd(b2, 0.8928551362746343397, 1e-12, "iauTcbtdb", "b2", ref status);
            viv(j, 0, "iauTcbtdb", "j", ref status);

        }

        static void t_tcgtt(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t c g t t
        **  - - - - - - - -
        **
        **  Test iauTcgtt function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTcgtt, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double t1, t2;
            int j;


            j = Sofa.iauTcgtt(2453750.5, 0.892862531, out t1, out t2);

            vvd(t1, 2453750.5, 1e-6, "iauTcgtt", "t1", ref status);
            vvd(t2, 0.8928551387488816828, 1e-12, "iauTcgtt", "t2", ref status);
            viv(j, 0, "iauTcgtt", "j", ref status);

        }

        static void t_tdbtcb(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ t d b t c b
        **  - - - - - - - - -
        **
        **  Test iauTdbtcb function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTdbtcb, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double b1, b2;
            int j;


            j = Sofa.iauTdbtcb(2453750.5, 0.892855137, out b1, out b2);

            vvd(b1, 2453750.5, 1e-6, "iauTdbtcb", "b1", ref status);
            vvd(b2, 0.8930195997253656716, 1e-12, "iauTdbtcb", "b2", ref status);
            viv(j, 0, "iauTdbtcb", "j", ref status);

        }

        static void t_tdbtt(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t d b t t
        **  - - - - - - - -
        **
        **  Test iauTdbtt function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTdbtt, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double t1, t2;
            int j;


            j = Sofa.iauTdbtt(2453750.5, 0.892855137, -0.000201, out t1, out t2);

            vvd(t1, 2453750.5, 1e-6, "iauTdbtt", "t1", ref status);
            vvd(t2, 0.8928551393263888889, 1e-12, "iauTdbtt", "t2", ref status);
            viv(j, 0, "iauTdbtt", "j", ref status);

        }

        static void t_tf2a(ref int status)
        /*
        **  - - - - - - -
        **   t _ t f 2 a
        **  - - - - - - -
        **
        **  Test iauTf2a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTf2a, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double a;
            int j;


            j = Sofa.iauTf2a('+', 4, 58, 20.2, out a);

            vvd(a, 1.301739278189537429, 1e-12, "iauTf2a", "a", ref status);
            viv(j, 0, "iauTf2a", "j", ref status);

        }

        static void t_tf2d(ref int status)
        /*
        **  - - - - - - -
        **   t _ t f 2 d
        **  - - - - - - -
        **
        **  Test iauTf2d function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTf2d, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double d;
            int j;


            j = Sofa.iauTf2d(' ', 23, 55, 10.9, out d);

            vvd(d, 0.9966539351851851852, 1e-12, "iauTf2d", "d", ref status);
            viv(j, 0, "iauTf2d", "j", ref status);

        }

        static void t_tpors(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t p o r s
        **  - - - - - - - -
        **
        **  Test iauTpors function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTpors, vvd, viv
        **
        **  This revision:  2017 October 21
        */
        {
            double xi, eta, ra, dec, az1, bz1, az2, bz2;
            int n;


            xi = -0.03;
            eta = 0.07;
            ra = 1.3;
            dec = 1.5;

            n = Sofa.iauTpors(xi, eta, ra, dec, out az1, out bz1, out az2, out bz2);

            vvd(az1, 1.736621577783208748, 1e-13, "iauTpors", "az1", ref status);
            vvd(bz1, 1.436736561844090323, 1e-13, "iauTpors", "bz1", ref status);

            vvd(az2, 4.004971075806584490, 1e-13, "iauTpors", "az2", ref status);
            vvd(bz2, 1.565084088476417917, 1e-13, "iauTpors", "bz2", ref status);

            viv(n, 2, "iauTpors", "n", ref status);

        }

        static void t_tporv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t p o r v
        **  - - - - - - - -
        **
        **  Test iauTporv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTporv, iauS2c, vvd, viv
        **
        **  This revision:  2017 October 21
        */
        {
            double xi, eta, ra, dec;
            double[] v;
            double[] vz1, vz2;
            int n;


            xi = -0.03;
            eta = 0.07;
            ra = 1.3;
            dec = 1.5;
            Sofa.iauS2c(ra, dec, out v);

            n = Sofa.iauTporv(xi, eta, v, out vz1, out vz2);

            vvd(vz1[0], -0.02206252822366888610, 1e-15,
                "iauTporv", "x1", ref status);
            vvd(vz1[1], 0.1318251060359645016, 1e-14,
                "iauTporv", "y1", ref status);
            vvd(vz1[2], 0.9910274397144543895, 1e-14,
                "iauTporv", "z1", ref status);

            vvd(vz2[0], -0.003712211763801968173, 1e-16,
                "iauTporv", "x2", ref status);
            vvd(vz2[1], -0.004341519956299836813, 1e-16,
                "iauTporv", "y2", ref status);
            vvd(vz2[2], 0.9999836852110587012, 1e-14,
                "iauTporv", "z2", ref status);

            viv(n, 2, "iauTporv", "n", ref status);

        }

        static void t_tpsts(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t p s t s
        **  - - - - - - - -
        **
        **  Test iauTpsts function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTpsts, vvd
        **
        **  This revision:  2017 October 21
        */
        {
            double xi, eta, raz, decz, ra, dec;


            xi = -0.03;
            eta = 0.07;
            raz = 2.3;
            decz = 1.5;

            Sofa.iauTpsts(xi, eta, raz, decz, out ra, out dec);

            vvd(ra, 0.7596127167359629775, 1e-14, "iauTpsts", "ra", ref status);
            vvd(dec, 1.540864645109263028, 1e-13, "iauTpsts", "dec", ref status);

        }

        static void t_tpstv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t p s t v
        **  - - - - - - - -
        **
        **  Test iauTpstv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTpstv, iauS2c, vvd
        **
        **  This revision:  2017 October 21
        */
        {
            double xi, eta, raz, decz;
            double[] vz, v;


            xi = -0.03;
            eta = 0.07;
            raz = 2.3;
            decz = 1.5;
            Sofa.iauS2c(raz, decz, out vz);

            Sofa.iauTpstv(xi, eta, vz, out v);

            vvd(v[0], 0.02170030454907376677, 1e-15, "iauTpstv", "x", ref status);
            vvd(v[1], 0.02060909590535367447, 1e-15, "iauTpstv", "y", ref status);
            vvd(v[2], 0.9995520806583523804, 1e-14, "iauTpstv", "z", ref status);

        }

        static void t_tpxes(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t p x e s
        **  - - - - - - - -
        **
        **  Test iauTpxes function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTpxes, vvd, viv
        **
        **  This revision:  2017 October 21
        */
        {
            double ra, dec, raz, decz, xi, eta;
            int j;


            ra = 1.3;
            dec = 1.55;
            raz = 2.3;
            decz = 1.5;

            j = Sofa.iauTpxes(ra, dec, raz, decz, out xi, out eta);

            vvd(xi, -0.01753200983236980595, 1e-15, "iauTpxes", "xi", ref status);
            vvd(eta, 0.05962940005778712891, 1e-15, "iauTpxes", "eta", ref status);

            viv(j, 0, "iauTpxes", "j", ref status);

        }

        static void t_tpxev(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t p x e v
        **  - - - - - - - -
        **
        **  Test iauTpxev function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTpxev, iauS2c, vvd
        **
        **  This revision:  2017 October 21
        */
        {
            double ra, dec, raz, decz, xi, eta;
            double[] v, vz;
            int j;


            ra = 1.3;
            dec = 1.55;
            raz = 2.3;
            decz = 1.5;
            Sofa.iauS2c(ra, dec, out v);
            Sofa.iauS2c(raz, decz, out vz);

            j = Sofa.iauTpxev(v, vz, out xi, out eta);

            vvd(xi, -0.01753200983236980595, 1e-15, "iauTpxev", "xi", ref status);
            vvd(eta, 0.05962940005778712891, 1e-15, "iauTpxev", "eta", ref status);

            viv(j, 0, "iauTpxev", "j", ref status);

        }

        static void t_tr(ref int status)
        /*
        **  - - - - -
        **   t _ t r
        **  - - - - -
        **
        **  Test iauTr function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTr, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r = { new double[3], new double[3], new double[3] };
            double[][] rt = { new double[3], new double[3], new double[3] };


            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            Sofa.iauTr(r, ref rt);

            vvd(rt[0][0], 2.0, 0.0, "iauTr", "11", ref status);
            vvd(rt[0][1], 3.0, 0.0, "iauTr", "12", ref status);
            vvd(rt[0][2], 3.0, 0.0, "iauTr", "13", ref status);

            vvd(rt[1][0], 3.0, 0.0, "iauTr", "21", ref status);
            vvd(rt[1][1], 2.0, 0.0, "iauTr", "22", ref status);
            vvd(rt[1][2], 4.0, 0.0, "iauTr", "23", ref status);

            vvd(rt[2][0], 2.0, 0.0, "iauTr", "31", ref status);
            vvd(rt[2][1], 3.0, 0.0, "iauTr", "32", ref status);
            vvd(rt[2][2], 5.0, 0.0, "iauTr", "33", ref status);

        }

        static void t_trxp(ref int status)
        /*
        **  - - - - - - -
        **   t _ t r x p
        **  - - - - - - -
        **
        **  Test iauTrxp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTrxp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p = new double[3];
            double[][] r = { new double[3], new double[3], new double[3] };
            double[] trp = new double[3];


            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            p[0] = 0.2;
            p[1] = 1.5;
            p[2] = 0.1;

            Sofa.iauTrxp(r, p, ref trp);

            vvd(trp[0], 5.2, 1e-12, "iauTrxp", "1", ref status);
            vvd(trp[1], 4.0, 1e-12, "iauTrxp", "2", ref status);
            vvd(trp[2], 5.4, 1e-12, "iauTrxp", "3", ref status);

        }

        static void t_trxpv(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t r x p v
        **  - - - - - - - -
        **
        **  Test iauTrxpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTrxpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r = { new double[3], new double[3], new double[3] };
            double[][] pv = { new double[3], new double[3] };
            double[][] trpv = { new double[3], new double[3] };


            r[0][0] = 2.0;
            r[0][1] = 3.0;
            r[0][2] = 2.0;

            r[1][0] = 3.0;
            r[1][1] = 2.0;
            r[1][2] = 3.0;

            r[2][0] = 3.0;
            r[2][1] = 4.0;
            r[2][2] = 5.0;

            pv[0][0] = 0.2;
            pv[0][1] = 1.5;
            pv[0][2] = 0.1;

            pv[1][0] = 1.5;
            pv[1][1] = 0.2;
            pv[1][2] = 0.1;

            Sofa.iauTrxpv(r, pv, ref trpv);

            vvd(trpv[0][0], 5.2, 1e-12, "iauTrxpv", "p1", ref status);
            vvd(trpv[0][1], 4.0, 1e-12, "iauTrxpv", "p1", ref status);
            vvd(trpv[0][2], 5.4, 1e-12, "iauTrxpv", "p1", ref status);

            vvd(trpv[1][0], 3.9, 1e-12, "iauTrxpv", "v1", ref status);
            vvd(trpv[1][1], 5.3, 1e-12, "iauTrxpv", "v2", ref status);
            vvd(trpv[1][2], 4.1, 1e-12, "iauTrxpv", "v3", ref status);

        }

        static void t_tttai(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t t t a i
        **  - - - - - - - -
        **
        **  Test iauTttai function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTttai, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double a1, a2;
            int j;


            j = Sofa.iauTttai(2453750.5, 0.892482639, out a1, out a2);

            vvd(a1, 2453750.5, 1e-6, "iauTttai", "a1", ref status);
            vvd(a2, 0.892110139, 1e-12, "iauTttai", "a2", ref status);
            viv(j, 0, "iauTttai", "j", ref status);

        }

        static void t_tttcg(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t t t c g
        **  - - - - - - - -
        **
        **  Test iauTttcg function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTttcg, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double g1, g2;
            int j;


            j = Sofa.iauTttcg(2453750.5, 0.892482639, out g1, out g2);

            vvd(g1, 2453750.5, 1e-6, "iauTttcg", "g1", ref status);
            vvd(g2, 0.8924900312508587113, 1e-12, "iauTttcg", "g2", ref status);
            viv(j, 0, "iauTttcg", "j", ref status);

        }

        static void t_tttdb(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t t t d b
        **  - - - - - - - -
        **
        **  Test iauTttdb function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTttdb, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double b1, b2;
            int j;


            j = Sofa.iauTttdb(2453750.5, 0.892855139, -0.000201, out b1, out b2);

            vvd(b1, 2453750.5, 1e-6, "iauTttdb", "b1", ref status);
            vvd(b2, 0.8928551366736111111, 1e-12, "iauTttdb", "b2", ref status);
            viv(j, 0, "iauTttdb", "j", ref status);

        }

        static void t_ttut1(ref int status)
        /*
        **  - - - - - - - -
        **   t _ t t u t 1
        **  - - - - - - - -
        **
        **  Test iauTtut1 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauTtut1, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double u1, u2;
            int j;


            j = Sofa.iauTtut1(2453750.5, 0.892855139, 64.8499, out u1, out u2);

            vvd(u1, 2453750.5, 1e-6, "iauTtut1", "u1", ref status);
            vvd(u2, 0.8921045614537037037, 1e-12, "iauTtut1", "u2", ref status);
            viv(j, 0, "iauTtut1", "j", ref status);

        }

        static void t_ut1tai(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ u t 1 t a i
        **  - - - - - - - - -
        **
        **  Test iauUt1tai function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauUt1tai, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double a1, a2;
            int j;


            j = Sofa.iauUt1tai(2453750.5, 0.892104561, -32.6659, out a1, out a2);

            vvd(a1, 2453750.5, 1e-6, "iauUt1tai", "a1", ref status);
            vvd(a2, 0.8924826385462962963, 1e-12, "iauUt1tai", "a2", ref status);
            viv(j, 0, "iauUt1tai", "j", ref status);

        }

        static void t_ut1tt(ref int status)
        /*
        **  - - - - - - - -
        **   t _ u t 1 t t
        **  - - - - - - - -
        **
        **  Test iauUt1tt function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauUt1tt, vvd, viv
        **
        **  This revision:  2013 October 3
        */
        {
            double t1, t2;
            int j;


            j = Sofa.iauUt1tt(2453750.5, 0.892104561, 64.8499, out t1, out t2);

            vvd(t1, 2453750.5, 1e-6, "iauUt1tt", "t1", ref status);
            vvd(t2, 0.8928551385462962963, 1e-12, "iauUt1tt", "t2", ref status);
            viv(j, 0, "iauUt1tt", "j", ref status);

        }

        static void t_ut1utc(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ u t 1 u t c
        **  - - - - - - - - -
        **
        **  Test iauUt1utc function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauUt1utc, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double u1, u2;
            int j;


            j = Sofa.iauUt1utc(2453750.5, 0.892104561, 0.3341, out u1, out u2);

            vvd(u1, 2453750.5, 1e-6, "iauUt1utc", "u1", ref status);
            vvd(u2, 0.8921006941018518519, 1e-12, "iauUt1utc", "u2", ref status);
            viv(j, 0, "iauUt1utc", "j", ref status);

        }

        static void t_utctai(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ u t c t a i
        **  - - - - - - - - -
        **
        **  Test iauUtctai function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauUtctai, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double u1, u2;
            int j;


            j = Sofa.iauUtctai(2453750.5, 0.892100694, out u1, out u2);

            vvd(u1, 2453750.5, 1e-6, "iauUtctai", "u1", ref status);
            vvd(u2, 0.8924826384444444444, 1e-12, "iauUtctai", "u2", ref status);
            viv(j, 0, "iauUtctai", "j", ref status);

        }

        static void t_utcut1(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ u t c u t 1
        **  - - - - - - - - -
        **
        **  Test iauUtcut1 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauUtcut1, vvd, viv
        **
        **  This revision:  2013 August 7
        */
        {
            double u1, u2;
            int j;


            j = Sofa.iauUtcut1(2453750.5, 0.892100694, 0.3341, out u1, out u2);

            vvd(u1, 2453750.5, 1e-6, "iauUtcut1", "u1", ref status);
            vvd(u2, 0.8921045608981481481, 1e-12, "iauUtcut1", "u2", ref status);
            viv(j, 0, "iauUtcut1", "j", ref status);

        }

        static void t_xy06(ref int status)
        /*
        **  - - - - - - -
        **   t _ x y 0 6
        **  - - - - - - -
        **
        **  Test iauXy06 function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauXy06, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y;


            Sofa.iauXy06(2400000.5, 53736.0, out x, out y);

            vvd(x, 0.5791308486706010975e-3, 1e-15, "iauXy06", "x", ref status);
            vvd(y, 0.4020579816732958141e-4, 1e-16, "iauXy06", "y", ref status);

        }

        static void t_xys00a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ x y s 0 0 a
        **  - - - - - - - - -
        **
        **  Test iauXys00a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauXys00a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y, s;


            Sofa.iauXys00a(2400000.5, 53736.0, out x, out y, out s);

            vvd(x, 0.5791308472168152904e-3, 1e-14, "iauXys00a", "x", ref status);
            vvd(y, 0.4020595661591500259e-4, 1e-15, "iauXys00a", "y", ref status);
            vvd(s, -0.1220040848471549623e-7, 1e-18, "iauXys00a", "s", ref status);

        }

        static void t_xys00b(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ x y s 0 0 b
        **  - - - - - - - - -
        **
        **  Test iauXys00b function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauXys00b, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y, s;


            Sofa.iauXys00b(2400000.5, 53736.0, out x, out y, out s);

            vvd(x, 0.5791301929950208873e-3, 1e-14, "iauXys00b", "x", ref status);
            vvd(y, 0.4020553681373720832e-4, 1e-15, "iauXys00b", "y", ref status);
            vvd(s, -0.1220027377285083189e-7, 1e-18, "iauXys00b", "s", ref status);

        }

        static void t_xys06a(ref int status)
        /*
        **  - - - - - - - - -
        **   t _ x y s 0 6 a
        **  - - - - - - - - -
        **
        **  Test iauXys06a function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauXys06a, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double x, y, s;


            Sofa.iauXys06a(2400000.5, 53736.0, out x, out y, out s);

            vvd(x, 0.5791308482835292617e-3, 1e-14, "iauXys06a", "x", ref status);
            vvd(y, 0.4020580099454020310e-4, 1e-15, "iauXys06a", "y", ref status);
            vvd(s, -0.1220032294164579896e-7, 1e-18, "iauXys06a", "s", ref status);

        }

        static void t_zp(ref int status)
        /*
        **  - - - - -
        **   t _ z p
        **  - - - - -
        **
        **  Test iauZp function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauZp, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[] p = new double[3];


            p[0] = 0.3;
            p[1] = 1.2;
            p[2] = -2.5;

            Sofa.iauZp(out p);

            vvd(p[0], 0.0, 0.0, "iauZp", "1", ref status);
            vvd(p[1], 0.0, 0.0, "iauZp", "2", ref status);
            vvd(p[2], 0.0, 0.0, "iauZp", "3", ref status);

        }

        static void t_zpv(ref int status)
        /*
        **  - - - - - -
        **   t _ z p v
        **  - - - - - -
        **
        **  Test iauZpv function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauZpv, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] pv = { new double[3], new double[3] };


            pv[0][0] = 0.3;
            pv[0][1] = 1.2;
            pv[0][2] = -2.5;

            pv[1][0] = -0.5;
            pv[1][1] = 3.1;
            pv[1][2] = 0.9;

            Sofa.iauZpv(out pv);

            vvd(pv[0][0], 0.0, 0.0, "iauZpv", "p1", ref status);
            vvd(pv[0][1], 0.0, 0.0, "iauZpv", "p2", ref status);
            vvd(pv[0][2], 0.0, 0.0, "iauZpv", "p3", ref status);

            vvd(pv[1][0], 0.0, 0.0, "iauZpv", "v1", ref status);
            vvd(pv[1][1], 0.0, 0.0, "iauZpv", "v2", ref status);
            vvd(pv[1][2], 0.0, 0.0, "iauZpv", "v3", ref status);

        }

        static void t_zr(ref int status)
        /*
        **  - - - - -
        **   t _ z r
        **  - - - - -
        **
        **  Test iauZr function.
        **
        **  Returned:
        **     status    int         FALSE = success, TRUE = fail
        **
        **  Called:  iauZr, vvd
        **
        **  This revision:  2013 August 7
        */
        {
            double[][] r = { new double[3], new double[3], new double[3] };


            r[0][0] = 2.0;
            r[1][0] = 3.0;
            r[2][0] = 2.0;

            r[0][1] = 3.0;
            r[1][1] = 2.0;
            r[2][1] = 3.0;

            r[0][2] = 3.0;
            r[1][2] = 4.0;
            r[2][2] = 5.0;

            Sofa.iauZr(out r);

            vvd(r[0][0], 0.0, 0.0, "iauZr", "00", ref status);
            vvd(r[1][0], 0.0, 0.0, "iauZr", "01", ref status);
            vvd(r[2][0], 0.0, 0.0, "iauZr", "02", ref status);

            vvd(r[0][1], 0.0, 0.0, "iauZr", "10", ref status);
            vvd(r[1][1], 0.0, 0.0, "iauZr", "11", ref status);
            vvd(r[2][1], 0.0, 0.0, "iauZr", "12", ref status);

            vvd(r[0][2], 0.0, 0.0, "iauZr", "20", ref status);
            vvd(r[1][2], 0.0, 0.0, "iauZr", "21", ref status);
            vvd(r[2][2], 0.0, 0.0, "iauZr", "22", ref status);

        }









        static void Main(string[] args)
        /*
        **  - - - - -
        **   m a i n
        **  - - - - -
        **
        **  This revision:  2021 April 18
        */
        {
            int status;


            /* If any command-line argument, switch to verbose reporting. */
            if (args.Length > 0)
            {
                if (args[0].Length > 0)
                {
                    verbose = 1;
                }
            }

            /* Preset the ref status to FALSE = success. */
            status = 0;

            /* Test all of the SOFA functions. */
            t_a2af(ref status);
            t_a2tf(ref status);
            t_ab(ref status);
            t_ae2hd(ref status);
            t_af2a(ref status);
            t_anp(ref status);
            t_anpm(ref status);
            t_apcg(ref status);
            t_apcg13(ref status);
            t_apci(ref status);
            t_apci13(ref status);
            t_apco(ref status);
            t_apco13(ref status);
            t_apcs(ref status);
            t_apcs13(ref status);
            t_aper(ref status);
            t_aper13(ref status);
            t_apio(ref status);
            t_apio13(ref status);
            t_atcc13(ref status);
            t_atccq(ref status);
            t_atci13(ref status);
            t_atciq(ref status);
            t_atciqn(ref status);
            t_atciqz(ref status);
            t_atco13(ref status);
            t_atic13(ref status);
            t_aticq(ref status);
            t_aticqn(ref status);
            t_atio13(ref status);
            t_atioq(ref status);
            t_atoc13(ref status);
            t_atoi13(ref status);
            t_atoiq(ref status);
            t_bi00(ref status);
            t_bp00(ref status);
            t_bp06(ref status);
            t_bpn2xy(ref status);
            t_c2i00a(ref status);
            t_c2i00b(ref status);
            t_c2i06a(ref status);
            t_c2ibpn(ref status);
            t_c2ixy(ref status);
            t_c2ixys(ref status);
            t_c2s(ref status);
            t_c2t00a(ref status);
            t_c2t00b(ref status);
            t_c2t06a(ref status);
            t_c2tcio(ref status);
            t_c2teqx(ref status);
            t_c2tpe(ref status);
            t_c2txy(ref status);
            t_cal2jd(ref status);
            t_cp(ref status);
            t_cpv(ref status);
            t_cr(ref status);
            t_d2dtf(ref status);
            t_d2tf(ref status);
            t_dat(ref status);
            t_dtdb(ref status);
            t_dtf2d(ref status);
            t_eceq06(ref status);
            t_ecm06(ref status);
            t_ee00(ref status);
            t_ee00a(ref status);
            t_ee00b(ref status);
            t_ee06a(ref status);
            t_eect00(ref status);
            t_eform(ref status);
            t_eo06a(ref status);
            t_eors(ref status);
            t_epb(ref status);
            t_epb2jd(ref status);
            t_epj(ref status);
            t_epj2jd(ref status);
            t_epv00(ref status);
            t_eqec06(ref status);
            t_eqeq94(ref status);
            t_era00(ref status);
            t_fad03(ref status);
            t_fae03(ref status);
            t_faf03(ref status);
            t_faju03(ref status);
            t_fal03(ref status);
            t_falp03(ref status);
            t_fama03(ref status);
            t_fame03(ref status);
            t_fane03(ref status);
            t_faom03(ref status);
            t_fapa03(ref status);
            t_fasa03(ref status);
            t_faur03(ref status);
            t_fave03(ref status);
            t_fk425(ref status);
            t_fk45z(ref status);
            t_fk524(ref status);
            t_fk52h(ref status);
            t_fk54z(ref status);
            t_fk5hip(ref status);
            t_fk5hz(ref status);
            t_fw2m(ref status);
            t_fw2xy(ref status);
            t_g2icrs(ref status);
            t_gc2gd(ref status);
            t_gc2gde(ref status);
            t_gd2gc(ref status);
            t_gd2gce(ref status);
            t_gmst00(ref status);
            t_gmst06(ref status);
            t_gmst82(ref status);
            t_gst00a(ref status);
            t_gst00b(ref status);
            t_gst06(ref status);
            t_gst06a(ref status);
            t_gst94(ref status);
            t_h2fk5(ref status);
            t_hd2ae(ref status);
            t_hd2pa(ref status);
            t_hfk5z(ref status);
            t_icrs2g(ref status);
            t_ir(ref status);
            t_jd2cal(ref status);
            t_jdcalf(ref status);
            t_ld(ref status);
            t_ldn(ref status);
            t_ldsun(ref status);
            t_lteceq(ref status);
            t_ltecm(ref status);
            t_lteqec(ref status);
            t_ltp(ref status);
            t_ltpb(ref status);
            t_ltpecl(ref status);
            t_ltpequ(ref status);
            t_moon98(ref status);
            t_num00a(ref status);
            t_num00b(ref status);
            t_num06a(ref status);
            t_numat(ref status);
            t_nut00a(ref status);
            t_nut00b(ref status);
            t_nut06a(ref status);
            t_nut80(ref status);
            t_nutm80(ref status);
            t_obl06(ref status);
            t_obl80(ref status);
            t_p06e(ref status);
            t_p2pv(ref status);
            t_p2s(ref status);
            t_pap(ref status);
            t_pas(ref status);
            t_pb06(ref status);
            t_pdp(ref status);
            t_pfw06(ref status);
            t_plan94(ref status);
            t_pmat00(ref status);
            t_pmat06(ref status);
            t_pmat76(ref status);
            t_pm(ref status);
            t_pmp(ref status);
            t_pmpx(ref status);
            t_pmsafe(ref status);
            t_pn(ref status);
            t_pn00(ref status);
            t_pn00a(ref status);
            t_pn00b(ref status);
            t_pn06a(ref status);
            t_pn06(ref status);
            t_pnm00a(ref status);
            t_pnm00b(ref status);
            t_pnm06a(ref status);
            t_pnm80(ref status);
            t_pom00(ref status);
            t_ppp(ref status);
            t_ppsp(ref status);
            t_pr00(ref status);
            t_prec76(ref status);
            t_pv2p(ref status);
            t_pv2s(ref status);
            t_pvdpv(ref status);
            t_pvm(ref status);
            t_pvmpv(ref status);
            t_pvppv(ref status);
            t_pvstar(ref status);
            t_pvtob(ref status);
            t_pvu(ref status);
            t_pvup(ref status);
            t_pvxpv(ref status);
            t_pxp(ref status);
            t_refco(ref status);
            t_rm2v(ref status);
            t_rv2m(ref status);
            t_rx(ref status);
            t_rxp(ref status);
            t_rxpv(ref status);
            t_rxr(ref status);
            t_ry(ref status);
            t_rz(ref status);
            t_s00a(ref status);
            t_s00b(ref status);
            t_s00(ref status);
            t_s06a(ref status);
            t_s06(ref status);
            t_s2c(ref status);
            t_s2p(ref status);
            t_s2pv(ref status);
            t_s2xpv(ref status);
            t_sepp(ref status);
            t_seps(ref status);
            t_sp00(ref status);
            t_starpm(ref status);
            t_starpv(ref status);
            t_sxp(ref status);
            t_sxpv(ref status);
            t_taitt(ref status);
            t_taiut1(ref status);
            t_taiutc(ref status);
            t_tcbtdb(ref status);
            t_tcgtt(ref status);
            t_tdbtcb(ref status);
            t_tdbtt(ref status);
            t_tf2a(ref status);
            t_tf2d(ref status);
            t_tpors(ref status);
            t_tporv(ref status);
            t_tpsts(ref status);
            t_tpstv(ref status);
            t_tpxes(ref status);
            t_tpxev(ref status);
            t_tr(ref status);
            t_trxp(ref status);
            t_trxpv(ref status);
            t_tttai(ref status);
            t_tttcg(ref status);
            t_tttdb(ref status);
            t_ttut1(ref status);
            t_ut1tai(ref status);
            t_ut1tt(ref status);
            t_ut1utc(ref status);
            t_utctai(ref status);
            t_utcut1(ref status);
            t_xy06(ref status);
            t_xys00a(ref status);
            t_xys00b(ref status);
            t_xys06a(ref status);
            t_zp(ref status);
            t_zpv(ref status);
            t_zr(ref status);

            /* Report, set up an appropriate exit status, and finish. */
            if (status != 0)
            {
                Console.Write("t_sofa_c validation failed!\n");
            }
            else
            {
                Console.Write("t_sofa_c validation successful\n");
            }

            Console.ReadKey();
        }
    }
}
