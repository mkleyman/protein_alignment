package com.Opitmizers;

import com.TransformationFunctions.LogarithmTransform;
import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.primitives.Doubles;

import java.util.Arrays;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class LogarithmOptimizer implements Optimizer {
    public static String name = "Log";

    public String getName() {
        return name;
    }
    public double[] optimizePearson(Aligner aligner, SplineDictionary splineDict,
                                               double threshold){
        double[] best = {0.0, 0.0};
        double opt_val;
        double max = Double.NEGATIVE_INFINITY;
        double[] params = {-3000,2000,1.0,3500};
        for(double b = params[0]; b<params[1]; b+=((params[0]+params[1])/50.0)){
            for(double m = params[2]; m<params[3]; m+= ((params[2]+params[3])/50.0)){
                LogarithmTransform poly = new LogarithmTransform(b,m);
                opt_val = aligner.align_polynomial_pearson(splineDict, poly, threshold);
                if (opt_val>max){
                    max = opt_val;
                    best = new double[]{b,m};
                }

            }
        }
        //System.out.println(Arrays.toString(best));
        //System.out.println(max);

        return Doubles.concat(new double[]{max, threshold}, best);
    }
}
