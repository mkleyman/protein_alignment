package com.TransformationFunctions;

import org.apache.commons.math3.analysis.UnivariateFunction;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class LogarithmTransform implements UnivariateFunction{
    double shift;
    double scale;

    public LogarithmTransform(double shift, double scale){
        this.scale = scale;
        this.shift = shift;
    }

    public double value(double x){
        return (scale*Math.log(x)+shift);
    }
}
