package com.TransformationFunctions;

import org.apache.commons.math3.analysis.UnivariateFunction;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class SqrtTransform implements UnivariateFunction{
    double shift;
    double scale;

    public SqrtTransform(double shift, double scale){
        this.scale = scale;
        this.shift = shift;
    }

    public double value(double x){
        return (scale*Math.sqrt(x)+shift);
    }
}
