package com.TransformationFunctions;

import org.apache.commons.math3.analysis.UnivariateFunction;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class ExpTransform implements UnivariateFunction {
    double shift;
    double scale;

    public ExpTransform(double shift, double scale){
        this.scale = scale;
        this.shift = shift;
    }

    public double value(double x){
        return (Math.exp(x*scale)+shift);
    }
}
