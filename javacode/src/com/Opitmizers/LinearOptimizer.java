package com.Opitmizers;

import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import java.util.Arrays;

/**
 * Created by mkleyman on 7/19/2016.
 */
public class LinearOptimizer extends Optimizer {
    public static String name = "Linear";
    private char mode ;

    public LinearOptimizer(){
        this.mode = 'p';
    }

    public LinearOptimizer(char mode){
        this.mode = mode;
    }

    public void setMode(char mode){
        this.mode = mode;
    }

    public String getName() {
        return name;
    }

    public double[] optimizePearson(Aligner aligner, SplineDictionary splineDict,
                                               double threshold){
        double[] best = {0.0, 0.0};
        double[] attempt;
        double opt_val;
        double max = Double.NEGATIVE_INFINITY;
        double[] params = {-5000.0,3000.0,0.2,250.0};
        double first_split = (params[1]-params[0])/500.0;
        double second_split = (params[3]-params[2])/500.0;
        for(double b = params[0]; b<params[1]; b+=first_split){
            for(double m = params[2]; m<params[3]; m+= (second_split)){
                attempt = new double[]{b,m};
                PolynomialFunction poly = new PolynomialFunction(attempt);
                opt_val = this.optVal(aligner,splineDict,poly,threshold,mode);
                if (opt_val>max){
                    max = opt_val;
                    best = attempt;
                }

            }
        }
        //System.out.println(Arrays.toString(best));
        //System.out.println(max);

        return Doubles.concat(new double[]{max,threshold},best);
    }
}
