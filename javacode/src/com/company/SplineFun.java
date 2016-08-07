package com.company;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import umontreal.ssj.functions.MathFunction;

/**
 * Created by mkleyman on 8/5/2016.
 */
public class SplineFun implements MathFunction {

    PolynomialSplineFunction poly;
    public SplineFun(PolynomialSplineFunction pol){
        this.poly = pol;
    }

    public double evaluate(double x){
        return poly.value(x);
    }
}
