package com.company;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import umontreal.ssj.functions.MathFunction;

/**
 * Created by mkleyman on 8/5/2016.
 */
public class SplineFun implements MathFunction {

    UnivariateFunction  poly;
    public SplineFun(UnivariateFunction pol){
        this.poly = pol;
    }



    public double evaluate(double x){
        return poly.value(x);

    }
}
