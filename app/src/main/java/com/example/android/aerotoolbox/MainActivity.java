/*
 * Copyright (C) 2012 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.example.android.aerotoolbox;

import android.content.Context;
import android.os.Bundle;
import android.support.v4.app.FragmentActivity;
import android.support.v4.app.FragmentTransaction;
import android.view.View;
import android.widget.EditText;
import android.widget.RadioButton;
import android.widget.Spinner;
import android.widget.TextView;
import android.widget.Toast;

public class MainActivity extends FragmentActivity
        implements MenuFragment.OnHeadlineSelectedListener {

    /** Called when the activity is first created. */
    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.news_articles);

        // Check whether the activity is using the layout version with
        // the fragment_container FrameLayout. If so, we must add the first fragment
        if (findViewById(R.id.fragment_container) != null) {

            // However, if we're being restored from a previous state,
            // then we don't need to do anything and should return or else
            // we could end up with overlapping fragments.
            if (savedInstanceState != null) {
                return;
            }

            // Create an instance of ExampleFragment
            MenuFragment firstFragment = new MenuFragment();

            // In case this activity was started with special instructions from an Intent,
            // pass the Intent's extras to the fragment as arguments
            firstFragment.setArguments(getIntent().getExtras());

            // Add the fragment to the 'fragment_container' FrameLayout
            getSupportFragmentManager().beginTransaction()
                    .add(R.id.fragment_container, firstFragment).commit();
        }
    }
/*
    Spinner spinner = (Spinner) findViewById(R.id.spinner);
    // Create an ArrayAdapter using the string array and a default spinner layout
    ArrayAdapter<CharSequence> adapter = ArrayAdapter.createFromResource(this,
            R.array.reynold_density, android.R.layout.simple_spinner_item);
    // Specify the layout to use when the list of choices appears
    adapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);

    // Apply the adapter to the spinner
    spinner.setAdapter(adapter);
*/
    public void onArticleSelected(int position) {
        // The user selected the headline of an article from the HeadlinesFragment

        // Capture the article fragment from the activity layout
        com.example.android.aerotoolbox.ContentFragment articleFrag = (com.example.android.aerotoolbox.ContentFragment)
                getSupportFragmentManager().findFragmentById(R.id.article_fragment);

        if (articleFrag != null) {
            // If article frag is available, we're in two-pane layout...

            // Call a method in the ArticleFragment to update its content
            articleFrag.updateArticleView(position);

        } else {
            // If the frag is not available, we're in the one-pane layout and must swap frags...

            // Create fragment and give it an argument for the selected article
            com.example.android.aerotoolbox.ContentFragment newFragment = new com.example.android.aerotoolbox.ContentFragment();
            Bundle args = new Bundle();
            args.putInt(com.example.android.aerotoolbox.ContentFragment.ARG_POSITION, position);
            newFragment.setArguments(args);
            FragmentTransaction transaction = getSupportFragmentManager().beginTransaction();

            // Replace whatever is in the fragment_container view with this fragment,
            // and add the transaction to the back stack so the user can navigate back
            transaction.replace(R.id.fragment_container, newFragment);
            transaction.addToBackStack(null);

            // Commit the transaction
            transaction.commit();
        }
    }

    private double NewtownIsentropic(double M, double A, double tol, int max_it){
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

    private void IsentropicSetMessage(int[] prefixes, int[] labels, double[] values, int title){
        String message;
        double M = values[0];
        TextView IsentropicTextView;
        double tol=1e-9;
        if (M>1+tol){
            message = "Supersonic Solution";
        }
        else if (M<1-tol){
            message = "Subsonic Solution";
        }
        else{
            message = "Sonic Solution";
        }

        IsentropicTextView = (TextView) findViewById(title);
        IsentropicTextView.setText(message);

        for(int i=0;i<5;i++) {
            message = getString(prefixes[i]) + ": " + String.format("%1$,.9f", values[i]);
            IsentropicTextView = (TextView) findViewById(labels[i]);
            IsentropicTextView.setText(message);
        }
    }

    private void IsentropicClearMessage(int [] labels, int title){
        String message="";
        TextView IsentropicTextView;

        IsentropicTextView = (TextView) findViewById(title);
        IsentropicTextView.setText(message);

        for(int i=0;i<5;i++) {
            IsentropicTextView = (TextView) findViewById(labels[i]);
            IsentropicTextView.setText(message);
        }
    }

    public void IsentropicCalculate(View view) {
        int labels [] = {R.id.isentropic_mach_output, R.id.isentropic_pressure_output,
                R.id.isentropic_temperature_output, R.id.isentropic_density_output, R.id.isentropic_area_output};
        int labels2 [] = {R.id.isentropic_mach_output2, R.id.isentropic_pressure_output2,
                R.id.isentropic_temperature_output2, R.id.isentropic_density_output2, R.id.isentropic_area_output2};
        int prefixes [] = {R.string.Mach, R.string.p_ratio, R.string.t_ratio, R.string.rho_ratio, R.string.a_ratio};
        int title = R.id.isentropic_label;
        int title2 = R.id.isentropic_label2;
        String message;
        TextView IsentropicTextView;
        double M, p_ratio,t_ratio,rho_ratio,a_ratio;
        RadioButton radio_mach, radio_pressure, radio_temperature, radio_density, radio_area;
        EditText text;

        try {

            radio_mach = (RadioButton) findViewById(R.id.radio_mach);
            radio_pressure = (RadioButton) findViewById(R.id.radio_pressure);
            radio_temperature = (RadioButton) findViewById(R.id.radio_temperature);
            radio_density = (RadioButton) findViewById(R.id.radio_density);
            radio_area = (RadioButton) findViewById(R.id.radio_area);
            if (radio_mach.isChecked()) {
                text = (EditText) findViewById(R.id.isentropic_mach_input);
                M = Double.parseDouble(text.getText().toString().trim());
                t_ratio = 1 / (1 + M * M / 5);
                p_ratio = Math.pow(t_ratio, 3.5);
                rho_ratio = p_ratio / t_ratio;
                a_ratio = 1/(M*Math.pow(1.2*t_ratio,3));
            }
            else if (radio_pressure.isChecked()) {
                text = (EditText) findViewById(R.id.isentropic_pressure_input);
                p_ratio = Double.parseDouble(text.getText().toString().trim());
                t_ratio = Math.pow(p_ratio, 1 / 3.5);
                M = Math.sqrt(5/t_ratio-5);
                rho_ratio = p_ratio / t_ratio;
                a_ratio = 1/(M*Math.pow(1.2*t_ratio,3));
            }
            else if (radio_temperature.isChecked()) {
                text = (EditText) findViewById(R.id.isentropic_temperature_input);
                t_ratio = Double.parseDouble(text.getText().toString().trim());
                M = Math.sqrt(5 / t_ratio - 5);
                p_ratio = Math.pow(t_ratio, 3.5);
                rho_ratio = p_ratio / t_ratio;
                a_ratio = 1/(M*Math.pow(1.2*t_ratio,3));
            }
            else if (radio_density.isChecked()) {
                text = (EditText) findViewById(R.id.isentropic_density_input);
                rho_ratio = Double.parseDouble(text.getText().toString().trim());
                t_ratio = Math.pow(rho_ratio, 0.4);
                p_ratio = t_ratio * rho_ratio;
                M = Math.sqrt(5 / t_ratio - 5);
                a_ratio = 1/(M*Math.pow(1.2*t_ratio,3));
            }
            else if (radio_area.isChecked()) {
                text = (EditText) findViewById(R.id.isentropic_area_input);
                a_ratio = Double.parseDouble(text.getText().toString().trim());
                if (a_ratio==1){
                    M=1;
                    t_ratio = 1 / (1 + M * M / 5);
                    p_ratio = Math.pow(t_ratio, 3.5);
                    rho_ratio = p_ratio / t_ratio;
                    IsentropicClearMessage(labels2, title2);
                }else{
                    double tol=1e-10; int max_it=100;
                    M=0;
                    M=NewtownIsentropic(M, a_ratio, tol, max_it);
                    t_ratio = 1 / (1 + M * M / 5);
                    p_ratio = Math.pow(t_ratio, 3.5);
                    rho_ratio = p_ratio / t_ratio;

                    double M2=10, t_ratio2, p_ratio2, rho_ratio2;
                    M2=NewtownIsentropic(M2, a_ratio, tol, max_it);
                    t_ratio2 = 1 / (1 + M2 * M2 / 5);
                    p_ratio2 = Math.pow(t_ratio2, 3.5);
                    rho_ratio2 = p_ratio2 / t_ratio2;
                    double values2 [] ={M2,p_ratio2,t_ratio2,rho_ratio2,a_ratio};
                    IsentropicSetMessage(prefixes, labels2, values2, title2);
                }

            }
            else{
                return;
            }
        }
        catch (NumberFormatException e){
            return;
        }
        double values [] ={M,p_ratio,t_ratio,rho_ratio,a_ratio};
        IsentropicSetMessage(prefixes, labels, values, title);

        if (!radio_area.isChecked()){
            IsentropicClearMessage(labels2, title2);
        }
    }

    public void ReynoldsClear(View view) {
        EditText text;
        text = (EditText) findViewById(R.id.reynold_density);
        text.setText("");
        text = (EditText) findViewById(R.id.reynold_velocity);
        text.setText("");
        text = (EditText) findViewById(R.id.reynold_length);
        text.setText("");
        text = (EditText) findViewById(R.id.reynold_viscosity);
        text.setText("");
    }

    public void ReynoldsReset(View view) {
        EditText text;
        Spinner spinner;
        int density_index, viscosity_index;

        spinner = (Spinner) findViewById(R.id.spinner_density);
        density_index = spinner.getSelectedItemPosition();
        spinner = (Spinner) findViewById(R.id.spinner_viscosity);
        viscosity_index = spinner.getSelectedItemPosition();

        String [] density_vals   = {"1.225", "0.0023769", "0.0765"};
        String [] viscosity_vals = {".00001789", ".0000003737"};

        text = (EditText) findViewById(R.id.reynold_density);
        text.setText(density_vals[density_index]);
        text = (EditText) findViewById(R.id.reynold_velocity);
        text.setText("100");
        text = (EditText) findViewById(R.id.reynold_length);
        text.setText("1");
        text = (EditText) findViewById(R.id.reynold_viscosity);
        text.setText(viscosity_vals[viscosity_index]);
    }

    public void ReynoldsCalculate(View view) {
        EditText text;
        double density, velocity, length, viscosity, Re;
        Spinner spinner;
        int density_index, velocity_index, length_index, viscosity_index;
        try {
            text = (EditText) findViewById(R.id.reynold_density);
            density = Double.parseDouble(text.getText().toString().trim());
            text = (EditText) findViewById(R.id.reynold_velocity);
            velocity = Double.parseDouble(text.getText().toString().trim());
            text = (EditText) findViewById(R.id.reynold_length);
            length = Double.parseDouble(text.getText().toString().trim());
            text = (EditText) findViewById(R.id.reynold_viscosity);
            viscosity = Double.parseDouble(text.getText().toString().trim());
        }
        catch (NumberFormatException e){
            return;
        }

        spinner = (Spinner) findViewById(R.id.spinner_density);
        density_index = spinner.getSelectedItemPosition();
        spinner = (Spinner) findViewById(R.id.spinner_velocity);
        velocity_index = spinner.getSelectedItemPosition();
        spinner = (Spinner) findViewById(R.id.spinner_length);
        length_index = spinner.getSelectedItemPosition();
        spinner = (Spinner) findViewById(R.id.spinner_viscosity);
        viscosity_index = spinner.getSelectedItemPosition();

        double [] density_units   = {1, 1.225/0.0023769, 1.225/0.0765};//0.06244897959,0.00194032653};
        double [] length_units    = {1, 1/3.280839895, 1/(3.280839895*12)};
        double [] velocity_units  = length_units;
        double [] viscosity_units = {1, 47.880106020829835};

        if (viscosity==0) {
            Re=Double.POSITIVE_INFINITY;
        }
        else{
            Re = density * velocity * length / viscosity
            * density_units[density_index]
            * velocity_units[velocity_index]
            * length_units[length_index]
            / viscosity_units[viscosity_index];
        }

        String message = "Reynold's Number: " + String.format("%1$,.2f", Re);

        TextView ReynoldTextView = (TextView) findViewById(R.id.reynold_number_value);
        ReynoldTextView.setText(message);
    }

    private double NewtownNormal_pt21(double M1, double pt21, double tol, int max_it){
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
    private double NewtownNormal_pt2p1(double M1, double pt2p1, double tol, int max_it){
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

    public class NormalValues {
        public double M1, M2, p21, pt21, pt2p1, t21, rho21;

        // constructor
        public NormalValues(double M1_in) {
            double M1sq;
            M1=M1_in;
            M1sq=M1_in*M1_in;

            M2=Math.sqrt((M1sq+5)/(7*M1sq-1));
            p21=(7*M1sq-1)/6;
            pt21=Math.pow((6*M1sq)/(M1sq+5),3.5)*Math.pow(6/(7*M1sq-1),2.5);
            pt2p1=Math.pow((6*M1sq)/(5),3.5)*Math.pow(6/(7*M1sq-1),2.5);
            t21=(7*M1sq-1)*(M1sq+5)/(36*M1sq);
            rho21=(6*M1sq)/(M1sq+5);
        }

        //methods
        public double [] values(){
            double [] normal_values={M1,M2,p21,pt21,pt2p1,t21,rho21};
            return normal_values;
        }
    }

    public void ShowToast(CharSequence text){
        Context context = getApplicationContext();
        int duration = Toast.LENGTH_SHORT;
        Toast toast = Toast.makeText(context, text, duration);
        toast.show();
    }

    public void NormalCalculate(View view) {
        double values[];
        double M1, M2, p21, pt21, pt2p1, t21, rho21;
        NormalValues OutputValues;
        RadioButton radio_m1, radio_m2, radio_p21, radio_pt21, radio_pt2p1, radio_t21, radio_rho21;
        EditText text;

        int labels[] = {R.id.normal_m1_output, R.id.normal_m2_output, R.id.normal_p21_output,
                R.id.normal_pt21_output, R.id.normal_pt2p1_output, R.id.normal_t21_output, R.id.normal_rho21_output};
        String prefixes[] = getResources().getStringArray(R.array.normal_array);
        int title = R.id.normal_label;

        try {
            radio_m1 = (RadioButton) findViewById(R.id.normal_radio_m1);
            radio_m2 = (RadioButton) findViewById(R.id.normal_radio_m2);
            radio_p21 = (RadioButton) findViewById(R.id.normal_radio_p21);
            radio_pt21 = (RadioButton) findViewById(R.id.normal_radio_pt21);
            radio_pt2p1 = (RadioButton) findViewById(R.id.normal_radio_pt2p1);
            radio_t21 = (RadioButton) findViewById(R.id.normal_radio_t21);
            radio_rho21 = (RadioButton) findViewById(R.id.normal_radio_rho21);

            if (radio_m1.isChecked()) {
                text = (EditText) findViewById(R.id.normal_m1_input);
                M1 = Double.parseDouble(text.getText().toString().trim());
                if (M1<1){
                    ShowToast("Upstream Mach number must be greater than 1");
                    return;
                }
            }
            else if (radio_m2.isChecked()) {
                text = (EditText) findViewById(R.id.normal_m2_input);
                M2 = Double.parseDouble(text.getText().toString().trim());
                if (M2>1){
                    ShowToast("Downstream Mach number must be less than 1");
                    return;
                }
                M1=Math.sqrt((M2*M2 + 5)/(7*M2*M2 - 1));
            }
            else if (radio_p21.isChecked()) {
                text = (EditText) findViewById(R.id.normal_p21_input);
                p21 = Double.parseDouble(text.getText().toString().trim());
                if (p21<1){
                    ShowToast("Static pressure ratio must be greater than 1");
                    return;
                }
                M1=Math.sqrt((6*p21+1)/7);
            }
            else if (radio_pt21.isChecked()) {
                text = (EditText) findViewById(R.id.normal_pt21_input);
                pt21 = Double.parseDouble(text.getText().toString().trim());
                if (pt21>1){
                    ShowToast("Total pressure ratio must be less than 1");
                    return;
                }
                M1=NewtownNormal_pt21(1.01,pt21,1e-10,100);
            }
            else if (radio_pt2p1.isChecked()) {
                text = (EditText) findViewById(R.id.normal_pt2p1_input);
                pt2p1 = Double.parseDouble(text.getText().toString().trim());
                if (pt2p1<1.892929159){
                    ShowToast("Value must be greater than Pt/P for sonic flow, ~1.8929");
                    return;
                }
                M1=NewtownNormal_pt2p1(1.01,pt2p1,1e-10,40);
            }
            else if (radio_t21.isChecked()) {
                text = (EditText) findViewById(R.id.normal_t21_input);
                t21 = Double.parseDouble(text.getText().toString().trim());
                if (t21<1){
                    ShowToast("Static temperature ratio must be greater than 1");
                    return;
                }
                M1=Math.sqrt((18*t21+6*Math.sqrt(9*t21*t21-17*t21 + 9)-17)/7);
            }
            else if (radio_rho21.isChecked()) {
                text = (EditText) findViewById(R.id.normal_rho21_input);
                rho21 = Double.parseDouble(text.getText().toString().trim());
                if (rho21<1){
                    ShowToast("density ratio must be greater than 1");
                    return;
                }
                if (rho21>=6){
                    ShowToast("Theoretical limit of density ratio is 6");
                    return;
                }
                M1=Math.sqrt(-(5*rho21)/(rho21 - 6));
            }
            else{
                return; // Nothing checked
            }

            OutputValues = new NormalValues(M1);
            values=OutputValues.values();
            NormalSetMessage(prefixes, labels, values, title);
        }

        catch (NumberFormatException e){
            return;
        }

    }
    private void NormalSetMessage(String[] prefixes, int[] labels, double[] values, int title){
        String message;
        TextView CurrentTextView;

        for(int i=0;i<7;i++) {
            message = prefixes[i] + ": " + String.format("%1$,.9f", values[i]);
            CurrentTextView = (TextView) findViewById(labels[i]);
            CurrentTextView.setText(message);
        }

    }

}
