package com.example.android.aerotoolbox;

/**
 * Created by Austin on 12/24/2015.
 */
public class Atmosphere {

    private double h, rho_ratio, T_ratio, P_ratio, rho, T, P, mu, a;
    private double T_sl = 518.67;                                 //Sea level Temperature
    private double P_sl = 2116.2;                                 //Sea level Pressure
    private double rho_sl = .0023769;                              //Sea level Density
    private double mu_sl = 3.737e-7;                              //Sea level Coefficient of Viscosity
    private double gamma = 1.4;                                  //The ratio of specific heats
    private double R = 1716;                                     //Univeral Gass Constant (EE units)
    double [] h_cut={36089.23885,82020.99738};
    double [] p_cut={472.7/P_sl,51.93/P_sl};
    double [] rho_cut={0.297116,0.032675884511227};
    // constructor
    public Atmosphere(){
    };

    //methods
    public void CalculateFromHeight(double h_in) {
        h=h_in;
        if (h < h_cut[0]) {
            T_ratio = 1-6.87558563248308e-6 * h;
            P_ratio = Math.pow(T_ratio, 5.2561);
        }
        else if(h > h_cut[1]){
            T_ratio=0.0000031626474339*h + 0.492571670940253;
            P_ratio=0.0245546508758944/ Math.pow(0.000004207534619643328 * h + 0.65507848568790397, 2847./250.);
        }
        else {
            T_ratio = .75189;
            P_ratio = 0.223358536281710 * Math.exp(4.80637968933164e-05 * (h_cut[0] - h));
        }
        rho_ratio = P_ratio / T_ratio;
        T = T_sl * T_ratio;
        P = P_sl * P_ratio;
        rho = rho_sl * rho_ratio;
        mu = mu_sl*Math.pow(T/T_sl,1.5)*(T_sl+198.72)/(T+198.72);
        a = Math.sqrt(gamma * R * T);
    }

    public void T_to_h(double T_in){
        T_ratio=T_in/T_sl;
        T_ratio_to_h(T_ratio);
    }

    private void T_ratio_to_h(double T_ratio){
        h=(1-T_ratio)/6.87558563248308e-6;
        CalculateFromHeight(h);
    }

    public void P_to_h(double P_in){
        P_ratio=P_in/P_sl;
        P_ratio_to_h(P_ratio);
    }

    private void P_ratio_to_h(double P_ratio_in){
        if (P_ratio_in>p_cut[0]) {
            h = (1 - Math.pow(P_ratio_in, 1 / 5.2561))/6.87558563248308e-6;
        }
        else if (P_ratio_in<p_cut[1]) {
            h=237668.8703478261/Math.pow(40.72548231511254*P_ratio_in,250./2847.) - 155691.7636826087;
        }
        else {
           h=h_cut[0]-20805.68046298183381*Math.log(4.4771067031830572759*P_ratio_in);
        }
        CalculateFromHeight(h);
    }

    public void rho_to_h(double rho_in){
        rho_ratio=rho_in/rho_sl;
        rho_ratio_to_h(rho_ratio);
    }

    private void rho_ratio_to_h(double rho_ratio_in){
        if (rho_ratio_in>rho_cut[0]) { //<11k meters
            T_ratio = Math.pow(rho_ratio_in, 1.0 / 4.2561);
            T_ratio_to_h(T_ratio);
        }
        else if (rho_ratio_in<rho_cut[1]) { // over 25k meters
            double h_guess=85000;
            double tol=1e-10;
            int max_it=100;
            h=rho_ratio_to_h(h_guess, rho_ratio_in, tol, max_it);
            CalculateFromHeight(h);
        }
        else { // 11-25 k meters
            T_ratio = .75189;
            P_ratio=rho_ratio_in*T_ratio;
            P_ratio_to_h(P_ratio);
        }
    }

    public double [] values_english(){
        return new double [] {h, P, T, rho, mu, a};
    }

    public double [] values_metric(){
        return new double [] {h/3.2808399, P*101325/2116.2, T/1.8, rho*1.225/0.0023769,
                              mu*47.880106020829835, a/3.2808399};
    }

    public double [] values_ratio(){
        return new double [] {P_ratio, T_ratio, rho_ratio};
    }

    public void mu_to_h(double mu_in){
        double T_guess=400;
        double tol=1e-15;
        int max_it=100;
        mu=mu_in;
        T=mu_to_T(mu,T_guess,tol, max_it);
        T_to_h(T);
    }


    private static double rho_ratio_to_h(double h, double rho_ratio, double tol, int max_it){
        //For high altitude only (above 82000 feet)
        double f, df, x1, t;
        for (int i=0; i<max_it;i++) {

            f = (0.2235381199111769*Math.exp(1.73 - 0.0000478535999272625 * h))/(0.0000031626474339*h + 0.492571670940253)-rho_ratio;
            df = - (0.0000106971037587218*Math.exp(1.73 - 0.0000478535999272625*h))/(0.0000031626474339*h + 0.492571670940253)
                    - (0.00000070697226131591415*Math.exp(1.73 - 0.0000478535999272625*h))/Math.pow(0.0000031626474339*h + 0.492571670940253,2);
            x1 = h-f/df;
            t = Math.abs(h-x1);
            h=x1;
            if (t<tol) {
                break;
            }
        }
        return h;
    }

    private double mu_to_T(double mu_in, double T, double tol, int max_it){
        double f, df, x1, t;
        double C=198.72;
        for (int i=0; i<max_it;i++) {
            f=mu_sl*Math.pow(T/T_sl,1.5)*(C + T_sl)/(C + T)-mu_in;
            df=3*mu_sl*Math.sqrt(T/T_sl)*(C + T_sl)/(2*T_sl*(C + T))
                    - (mu_sl*Math.pow(T/T_sl,1.5)*(C + T_sl))/Math.pow(C + T,2.);
            x1 = T-f/df;
            t = Math.abs(T-x1);
            T=x1;
            if (t<tol) {
                break;
            }
        }
        return T;
    }
}