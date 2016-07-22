package com.TransformationFunctions;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Logit;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class LogitTransform  implements UnivariateFunction{
    private Logit logit;
    private double shift;

    public LogitTransform(double shift, double min, double max){
        this.logit = new Logit(min,max);
        this.shift = shift;
    }

    public double value(double x){
        return (logit.value(x)+shift);
    }


}
