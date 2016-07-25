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
        double first_split = (params[1]-params[0])/100;
        double second_split = (params[3]-params[2])/100.0;
        double third_split =  (params[5]-params[4])/100.0;
        double fourth_split = (params[7]-params[6])/100.0;
        for(double b = params[0]; b<params[1]; b+=first_split){
            for(double m = params[2]; m<params[3]; m+= (second_split)){
                for(double a=params[4];a<params[5]; a+= third_split) {
                    for(double c=params[6];a<params[7]; a+= fourth_split) {
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
