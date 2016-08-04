package com.Opitmizers;

import com.TransformationFunctions.LogarithmTransform;
import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.primitives.Doubles;

import java.util.Arrays;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class LogarithmOptimizer extends Optimizer {
    public static String name = "Log";
    private char mode ;

    public LogarithmOptimizer(){
        this.mode = 'p';
    }

    public LogarithmOptimizer(char mode){
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
        double opt_val;
        double max = Double.NEGATIVE_INFINITY;
        double[] params = {-8000,3000,1.0,5000};
        double first_split = (params[1]-params[0])/500.0;
        double second_split = (params[3]-params[2])/500.0;
        for(double b = params[0]; b<params[1]; b+=first_split){
            for(double m = params[2]; m<params[3]; m+= (second_split)){
                LogarithmTransform poly = new LogarithmTransform(b,m);
                opt_val = this.optVal(aligner,splineDict,poly,threshold,mode);
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
