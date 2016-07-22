package com.Opitmizers;

import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import java.util.Arrays;

/**
 * Created by mkleyman on 7/19/2016.
 */
public class LinearOptimizer implements Optimizer {
    public static String name = "Linear";

    public String getName() {
        return name;
    }

    public double[] optimizePearson(Aligner aligner, SplineDictionary splineDict,
                                               double threshold){
        double[] best = {0.0, 0.0};
        double[] attempt;
        double opt_val;
        double max = Double.NEGATIVE_INFINITY;
        double[] params = {-3000.0,2000.0,0.2,60.0};
        for(double b = params[0]; b<params[1]; b+=((params[0]+params[1])/50.0)){
            for(double m = params[2]; m<params[3]; m+= ((params[2]+params[3])/50.0)){
                attempt = new double[]{b,m};
                PolynomialFunction poly = new PolynomialFunction(attempt);
                opt_val = aligner.align_polynomial_pearson(splineDict, poly, threshold);
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
