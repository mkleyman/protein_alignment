package com.Opitmizers;

import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import java.util.Arrays;
/**
 * Created by mkleyman on 7/20/2016.
 */
public class CubicOptmizer implements Optimizer {
    public static String name = "Cubic";

    public double[] optimizePearson(Aligner aligner, SplineDictionary splineDict,
                                    double threshold){
        double[] best = {0.0, 0.0, 0.0};
        double[] attempt;
        double opt_val;
        double max = Double.NEGATIVE_INFINITY;
        double[] params = {-3000,2000,-180,180,-30,30,-7.0,7.0};
        for(double b = params[0]; b<params[1]; b+=((params[0]+params[1])/10.0)){
            for(double m = params[2]; m<params[3]; m+= ((params[2]+params[3])/10.0)){
                for(double a=m+params[4];a<params[5]; a+= ((params[4]+params[5])/10.0)) {
                    for(double c=params[6];a<params[7]; a+= ((params[6]+params[7])/10.0)) {
                        attempt = new double[]{b, m, a, c};
                        PolynomialFunction poly = new PolynomialFunction(attempt);
                        opt_val = aligner.align_polynomial_pearson(splineDict, poly, threshold);
                        if (opt_val > max) {
                            max = opt_val;
                            best = attempt;
                        }
                    }
                }

            }
        }
        //System.out.println(Arrays.toString(best));
        //System.out.println(max);

        return Doubles.concat(new double[]{max, threshold}, best);
    }


    public String getName() {
        return name;
    }
}
