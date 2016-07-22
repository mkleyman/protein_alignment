package com.company;

import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.*;

/**
 * Created by mkleyman on 7/18/2016.
 */
public class SplineDictionary {


    private Map<String, PolynomialSplineFunction> splineDict;

    public SplineDictionary(TreeBasedTable<String,Double,Double> expressionTable, AkimaSplineInterpolator splineMaker,
                            Set<String> proteinSet ){
        this.splineDict = new HashMap<>();
        double[] times = Doubles.toArray(expressionTable.columnKeySet());
        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                double[] expressionVals = Doubles.toArray(expressionTable.row(row).values());
                splineDict.put(row, splineMaker.interpolate(times, expressionVals));
            }
        }
    }

    public SplineDictionary(TreeBasedTable<String,Double,Double> expressionTable,
                            AkimaSplineInterpolator splineMaker, Set<String> proteinSet, long rand){
        this.splineDict = new HashMap<>();
        double[] times = Doubles.toArray(expressionTable.columnKeySet());
        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                List<Double> randomExpressions = new ArrayList<Double>(expressionTable.row(row).size());
                randomExpressions.addAll(expressionTable.row(row).values());
                Collections.shuffle(randomExpressions, new Random(rand));
                double[] expressionVals = Doubles.toArray(randomExpressions);
                this.splineDict.put(row, splineMaker.interpolate(times, expressionVals));
            }
        }
    }

    public PolynomialSplineFunction getSpline(String homolog){
        if(this.splineDict.containsKey(homolog)){
            return this.splineDict.get(homolog);
        }
        else return null;
    }
}
