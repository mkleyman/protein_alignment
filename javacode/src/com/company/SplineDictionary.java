package com.company;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.interpolation.*;
import umontreal.ssj.functionfit.SmoothingCubicSpline;
import umontreal.ssj.functions.MathFunction;


import javax.script.*;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by mkleyman on 7/18/2016.
 */
public class SplineDictionary {


    private Map<String,MathFunction> splineDict;

    public SplineDictionary(){
        splineDict = new HashMap<>();
    }

    public SplineDictionary(TreeBasedTable<String,Double,Double> expressionTable, AkimaSplineInterpolator splineMaker,
                            Set<String> proteinSet ){
        this.splineDict = new HashMap<>();
        double[] times = Doubles.toArray(expressionTable.columnKeySet());
        UnivariateInterpolator interpolator = new LoessInterpolator(0.6,2);
        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                double[] expressionVals = Doubles.toArray(expressionTable.row(row).values());
                //SmoothingCubicSpline scs = new SmoothingCubicSpline(times,expressionVals,20);
                //splineDict.put(row, new SmoothingCubicSpline(times,expressionVals, 0.3));
                //splineDict.put(row, BSpline.createApproxBSpline(times,expressionVals,6,3));
                splineDict.put(row,new SplineFun(interpolator.interpolate(times,expressionVals)));
                //splineDict.put(row, new SplineFun(new CubicInterpolation(times,expressionVals,0)));
                //System.out.println(splineDict.get(row).getKnots().length);
            }
        }
    }

    public void retainAll(Collection<String> homologs){
        splineDict.keySet().retainAll(homologs);
    }


    public SplineDictionary(TreeBasedTable<String,Double,Double> expressionTable,
                            AkimaSplineInterpolator splineMaker, Set<String> proteinSet, long rand){
        this.splineDict = new HashMap<>();
        double[] times = Doubles.toArray(expressionTable.columnKeySet());

        UnivariateInterpolator interpolator = new LoessInterpolator(0.6,2);
        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                List<Double> randomExpressions = new ArrayList<>(expressionTable.row(row).size());
                randomExpressions.addAll(expressionTable.row(row).values());
                Collections.shuffle(randomExpressions, new Random(rand));
                double[] expressionVals = Doubles.toArray(randomExpressions);
                //this.splineDict.put(row, new SmoothingCubicSpline(times,expressionVals,0.1));
                splineDict.put(row,new SplineFun(interpolator.interpolate(times,expressionVals)));

            }
        }
    }

    public void addTable(TreeBasedTable<String,Double,Double> expressionTable,
                            Set<String> proteinSet ){
        double[] times = Doubles.toArray(expressionTable.columnKeySet());
        UnivariateInterpolator interpolator = new LoessInterpolator(0.6,2);
        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                double[] expressionVals = Doubles.toArray(expressionTable.row(row).values());
                //SmoothingCubicSpline scs = new SmoothingCubicSpline(times,expressionVals,20);
                //splineDict.put(row, new SmoothingCubicSpline(times,expressionVals, 0.3));
               // splineDict.put(row, BSpline.createApproxBSpline(times, expressionVals, 6,3));
                //splineDict.put(row, new SplineFun(new CubicInterpolation(times,expressionVals,0)));

                splineDict.put(row,new SplineFun(interpolator.interpolate(times,expressionVals)));

                //System.out.println(splineDict.get(row).getKnots().length);
            }
        }
    }


    public void addTable(TreeBasedTable<String,Double,Double> expressionTable,
                            Set<String> proteinSet, long rand){

        double[] times = Doubles.toArray(expressionTable.columnKeySet());
        UnivariateInterpolator interpolator = new LoessInterpolator(0.6,2);
        for(String row: expressionTable.rowKeySet()){
            if(proteinSet.contains(row)) {
                List<Double> randomExpressions = new ArrayList<>(expressionTable.row(row).size());
                randomExpressions.addAll(expressionTable.row(row).values());
                Collections.shuffle(randomExpressions, new Random(rand));
                double[] expressionVals = Doubles.toArray(randomExpressions);
                //this.splineDict.put(row, new SmoothingCubicSpline(times,expressionVals,0.05));
                splineDict.put(row,new SplineFun(interpolator.interpolate(times,expressionVals)));

            }
        }
    }

    public void parseInFile(String filename) throws IOException{
        File file = new File(filename);
        Scanner scan = new Scanner(file);
        LoessInterpolator interpolator = new LoessInterpolator();
        String[] header = scan.nextLine().split(",");

        double[] times = new double[header.length-1];

        //read in times for column names
        for(int i=1;i<header.length;i++){
            times[i-1] = Double.parseDouble(header[i]);
        }
        //read in protein expression
        double[] expression = new double[times.length];
        while(scan.hasNextLine()){
            String[] row = scan.nextLine().split(",");
            String protein = row[0];
            for(int i=1;i<row.length;i++){
                expression[i-1] = Double.parseDouble(row[i]);
            }
            System.out.println(Arrays.toString(times));
            ///splineDict.put(protein, new SplineFun(new KrigingInterpolation1D(times,expression)));

        }
        scan.close();
        System.out.println("splines parsed in");
    }



    public MathFunction getSpline(String homolog){
        if(this.splineDict.containsKey(homolog)){
            return this.splineDict.get(homolog);
        }
        else return null;
    }


}
