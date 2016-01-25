package com.example.android.aerotoolbox;

/**
 * Created by Austin on 12/20/2015.
 * Collection of various functions, mostly Newton's method, for isentropic flow, normal shocks and
 * p-m flow. Not enough to warrant separate classes like ObliqueShock, Reynolds, or Atmosphere.
 */
public class AeroCalc {

    public static double Isentropic(double M, double A, double tol, int max_it){
        double f, df, x1, t;
        for (int i=0; i<max_it;i++) {
            f = Math.pow(M,6)+15*Math.pow(M,4)+75*Math.pow(M,2)-216*M*A+125;
            df = 6*Math.pow(M,5)+60*Math.pow(M,3)+150*M-216*A;
            x1 = M-f/df;
            t = Math.abs(M-x1);
            M=x1;
            if (t<tol) {
                break;
            }
        }
        return M;
    }

    public static double Normal_pt21(double M1, double pt21, double tol, int max_it){
        double f, df, x1, t, M1sq;
        M1sq=M1*M1;
        for (int i=0; i<max_it;i++) {
            f = Math.pow((6*M1sq)/(M1sq+5),3.5)*Math.pow(6/(7*M1sq-1),2.5)-pt21;
            //PROBLEM HERE WITH  f or df;
            df = 7/2*(6/(M1sq+5)-(6*M1sq)/Math.pow(M1sq+5,2))*Math.pow(6/(7*M1sq - 1),5/2)*Math.pow((6*M1sq)/(M1sq + 5),5/2)
                    -105*Math.pow(6/(7*M1sq - 1),1.5)*Math.pow(6*M1sq/(M1sq + 5),3.5)/Math.pow(7*M1sq - 1,2);
            x1 = M1sq-f/df;
            t = Math.abs(M1sq-x1);
            M1sq=x1;
            if (t<tol) {
                break;
            }
        }
        return Math.sqrt(M1sq);
    }

    public static double Normal_pt2p1(double M1, double pt2p1, double tol, int max_it){
        double f, df, x1, t, M1sq;
        //M1sq=10;
        M1sq=M1*M1;
        for (int i=0; i<max_it;i++) {
            f = Math.pow((6*M1sq)/5,3.5)*Math.pow(6/(7*M1sq-1),2.5)-pt2p1;
            df = 21/5*Math.pow(6*M1sq/5,2.5)*Math.pow(6/(7*M1sq - 1),2.5)
                    -105*Math.pow(6*M1sq/5,3.5)*Math.pow(6/(7*M1sq - 1),1.5)/Math.pow(7*M1sq-1,2);
            x1 = M1sq-f/df;
            t = Math.abs(M1sq-x1);
            M1sq=x1;
            if (t<tol) {
                break;
            }
        }
        return Math.sqrt(M1sq);
    }

    public static double Expansion_nu_to_M(double M_guess, double nu, double tol, int max_it){
        double f, df, x1, t, Ms, c, sMs;
        nu *= Math.PI/180;
        Ms=M_guess*M_guess-1;
        c=Math.sqrt(6);
        for (int i=0; i<max_it;i++) {
            sMs=Math.sqrt(Ms);
            f = c*Math.atan(sMs/c)-Math.atan(sMs)-nu; //harcoded gamma=1.4
            df = c/(12*(Ms/6+1)*(sMs/c))-1/(2*sMs*(Ms+1));
            x1 = Ms-f/df;
            t = Math.abs(Ms-x1);
            Ms=x1;
            if (t<tol) {
                break;
            }
        }
        return Math.sqrt(Ms+1);
    }

    public static double M_to_nu(double M_in) {
        double nu, sMs, c;
        sMs=Math.sqrt(M_in*M_in-1);
        c=Math.sqrt(6);
        nu=180/Math.PI*(c*Math.atan(sMs/c)-Math.atan(sMs));
        return nu;
    }

}
