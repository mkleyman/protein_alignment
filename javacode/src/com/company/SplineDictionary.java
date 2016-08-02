package com.company;

import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import umontreal.ssj.functionfit.SmoothingCubicSpline;
import umontreal.ssj.functions.MathFunction;


import javax.script.*;

import java.util.*;

/**
 * Created by mkleyman on 7/18/2016.
 */
public class SplineDictionary {


    private Map<String,MathFunction> splineDict;

    public SplineDictionary(TreeBasedTable<String,Double,Double> expressionTable, AkimaSplineInterpolator splineMaker,
                            Set<String> proteinSet ){
        this.splineDict = new HashMap<>();
        double[] times = Doubles.toArray(expressionTable.columnKeySet());
        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                double[] expressionVals = Doubles.toArray(expressionTable.row(row).values());
                //SmoothingCubicSpline scs = new SmoothingCubicSpline(times,expressionVals,20);
                splineDict.put(row, new SmoothingCubicSpline(times,expressionVals, 0.1));

                //System.out.println(splineDict.get(row).getKnots().length);
            }
        }
    }

    public SplineDictionary(TreeBasedTable<String,Double,Double> expressionTable,
                            AkimaSplineInterpolator splineMaker, Set<String> proteinSet, long rand){
        this.splineDict = new HashMap<>();
        double[] times = Doubles.toArray(expressionTable.columnKeySet());




        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                List<Double> randomExpressions = new ArrayList<>(expressionTable.row(row).size());
                randomExpressions.addAll(expressionTable.row(row).values());
                Collections.shuffle(randomExpressions, new Random(rand));
                double[] expressionVals = Doubles.toArray(randomExpressions);
                this.splineDict.put(row, new SmoothingCubicSpline(times,expressionVals,0.1));
            }
        }
    }

    public MathFunction getSpline(String homolog){
        if(this.splineDict.containsKey(homolog)){
            return this.splineDict.get(homolog);
        }
        else return null;
    }


}
