package com.Opitmizers;

import com.TransformationFunctions.LogitTransform;
import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.primitives.Doubles;

import java.util.Arrays;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class LogitOptimizer extends Optimizer {
    public static String name = "Logit";
    private char mode ;

    public LogitOptimizer(){
        this.mode = 'p';
    }

    public LogitOptimizer(char mode){
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
        double[] best = {0.0, 0.0, 0.0};
        double opt_val;
        double max = Double.NEGATIVE_INFINITY;
        double[] params = {-3000,2000,1,16,47,100};
        double first_split = (params[1]-params[0])/100.0;
        double second_split = (params[3]-params[2])/100.0;
        double third_split = (params[5]-params[4])/100.0;
        for(double b = params[0]; b<params[1]; b+=first_split){
            for(double m = params[2]; m<params[3]; m+= (second_split)){
                for(double a=m+params[4];a<params[5]; a+= third_split) {
                    LogitTransform poly = new LogitTransform (b,m,a);
                    opt_val = this.optVal(aligner,splineDict,poly,threshold,mode);
                    if (opt_val > max) {
                        max = opt_val;
                        best = new double[]{b, m, a};;
                    }
                }

            }
        }
        //System.out.println(Arrays.toString(best));
        //System.out.println(max);

        return Doubles.concat(new double[]{max, threshold}, best);
    }
}
