package com.example.android.aerotoolbox;

/**
 * Created by Austin on 12/29/2015.
 */
public class ObliqueShock {

    private double M1, M2, delta, theta, delta_max, theta_max, P_ratio, T_ratio, rho_ratio;
    //theta_max is not a true max theta, it is the theta for max delta

    // constructor
    public ObliqueShock(){};

    public void CalcMaxTheta(double M1)
    {
        double M1s=M1*M1;
        theta_max=180/Math.PI*Math.asin(Math.sqrt((3*M1s-5+Math.sqrt(3*(3*M1s*M1s+4*M1s+20)))/(7*M1s))); //eqn 168, gamma hardcoded
    }

    public double CalcDelta(double M1_in, double theta_in){

        double M1sq=M1_in*M1_in;
        double sinsq=Math.pow(Math.sin(theta_in*Math.PI/180),2);
        double M1st=M1sq*sinsq;

        return 180/Math.PI*Math.atan((5*(M1st-1))/(Math.tan(Math.PI/180*theta_in)*(5+M1sq*(6-5*sinsq))));
    }

    public void CalcFromM1theta(double M1_in, double theta_in){
        M1=M1_in;
        theta=theta_in;
        double M1sq=M1*M1;
        double sinsq=Math.pow(Math.sin(theta*Math.PI/180),2);
        double M1st=M1sq*sinsq;

        CalcMaxTheta(M1); //Calculate theta_max
        delta_max=CalcDelta(M1_in,theta_max);
        delta=CalcDelta(M1_in,theta_in);

        //delta=180/Math.PI*Math.atan((5*(M1st-1))/(Math.tan(Math.PI/180*theta)*(5+M1sq*(6-5*sinsq))));

        P_ratio=(7*M1st-1)/6;
        rho_ratio=6*M1st/(M1st+5);
        T_ratio=P_ratio/rho_ratio;
        M2=Math.sqrt((36*M1sq*M1st-5*(M1st-1)*(7*M1st+5))/((7*M1st-1)*(M1st+5))); //has gama=1.4

    }

    public void CalcFromdeltatheta(double delta_in, double theta_in){
        theta=theta_in;
        delta=delta_in;

        double cott=1/Math.tan(Math.PI/180.*theta);
        double tand=Math.tan(Math.PI / 180. * delta);
        double sin2t=Math.sin(2 * Math.PI / 180. * theta);
        double cos2t=Math.cos(2 * Math.PI / 180. * theta);
        //double temp1 = (5*sin2t-tand*(7+5*cos2t));
        //double temp2 = 10*(cott+tand);
        double temp=(10*(cott+tand))/(5*sin2t-tand*(7+5*cos2t));
        if (temp<0){
            delta=-371;
            return;
        }
        M1=Math.sqrt(temp); //eqn 148b
        CalcFromM1theta(M1,theta_in);
    }

    public double get_delta(){
        return delta;
    }


    public double get_delta_max(){
        return delta_max;
    }

    public double get_theta_max(){
        return theta_max;
    }

    public double get_theta(){
        return theta;
    }

    public void CalcFromM1delta(double M1_in, double delta_in, boolean strong){
        M1=M1_in;
        delta=delta_in;

        CalcMaxTheta(M1); //Calculate theta_max
        delta_max=CalcDelta(M1_in,theta_max);
        if (delta>delta_max){
            delta=-370;
            return;
        }

        double M1sq=M1*M1;
        double sd=Math.pow(Math.sin(delta * Math.PI / 180.), 2);
        double y=1.4;
        double b=-(M1sq+2)/M1sq-y*sd;

        if(strong){ //
            theta=Newton(M1, delta, 1);
        }
        else{ //weak, lower theta value make initial guess -b/3
            theta=Newton(M1, delta, -b/3); //average of roots of derivative
        }

        CalcFromM1theta(M1, theta);
    }

    private double Newton(double M1, double delta, double x){
        double f, df, x1, t, M1s, M1q, sd, b, c, d, theta;
        double tol=1e-10;
        double y=1.4;
        int max_it=100;

        M1s=M1*M1;
        M1q=M1s*M1s;
        sd=Math.pow(Math.sin(delta*Math.PI/180.),2);

        b=-(M1s+2)/M1s-y*sd;
        c=(2*M1s+1)/M1q+(Math.pow(y+1,2)/4+(y-1)/M1s)*sd;
        d=-(1-sd)/M1q;

        for (int i=0; i<max_it;i++) {
            f=x*x*x+b*x*x+c*x+d;
            df=3*x*x+2*b*x+c;
            x1 = x-f/df;
            t = Math.abs(x-x1);
            x=x1;
            if (t<tol) {
                break;
            }
        }
        //theta=Math.asin(Math.sqrt(delta) * 180. / Math.PI);
        theta=Math.asin(Math.sqrt(x)) * 180. / Math.PI;
        return theta;
    }

    public double [] values(){
        double values [] = {M1, M2, delta, theta, delta_max, P_ratio, T_ratio, rho_ratio};
        return values;
    }

}
