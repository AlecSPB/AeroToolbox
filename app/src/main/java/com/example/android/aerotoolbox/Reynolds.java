package com.example.android.aerotoolbox;

/**
 * Created by Austin on 1/24/2016.
 */
public class Reynolds {

    public double Re, delta_s;
    private double cf, twall, ufric;

    private double [] density_units   = {1, 1.225/0.0023769, 1.225/0.0765};//0.06244897959,0.00194032653};
    private double [] length_units    = {1, 1/3.280839895, 1/(3.280839895*12)};
    private double [] velocity_units  = length_units;
    private double [] viscosity_units = {1, 47.880106020829835};

    // constructor
    public Reynolds(){};

    public void Calculate(double density, double velocity, double length, double viscosity,
                        int density_index, int velocity_index, int length_index, int viscosity_index){

        if (viscosity==0) {
            Re=Double.POSITIVE_INFINITY;
        }
        else{
            density*= density_units[density_index];
            velocity*= velocity_units[velocity_index];
            length*= length_units[length_index];
            viscosity/= viscosity_units[viscosity_index];

            Re = density * velocity * length / viscosity;

        }

        if (density==0) {
            delta_s = Double.POSITIVE_INFINITY;
        }
        else {
            cf = .026 / Math.pow(Re, 1 / 7.);
            twall = cf * density * velocity * velocity / 2;
            ufric = Math.sqrt(twall / density);
            delta_s = 1 * viscosity / (ufric * density); //in meters
            delta_s /= length_units[length_index];
        }


    }
}
