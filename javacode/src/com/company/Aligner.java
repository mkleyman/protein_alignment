package com.company;

import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;


import java.util.*;

/**
 * Created by mkleyman on 7/18/2016.
 */
public class Aligner {
    private final double[] compTimes;
    private final Set<String> proteinList;
    private final Map<String,String> homologDict;
    private final double[] times;
    private final double[] testTimes;
    private final Table<String,Double,Double> refTable;
    private static PearsonsCorrelation pearson = new PearsonsCorrelation();
    private static SpearmansCorrelation spearman = new SpearmansCorrelation();




    public Aligner(double[] compTimes, Set<String> proteinList, Map<String,String> homologDict,
                   double[] times, Table<String, Double, Double> refTable, double[] testTimes){
        this.compTimes= compTimes;
        this.proteinList = proteinList;
        this.homologDict = homologDict;
        this.times = times;
        this.refTable = refTable;
        this.testTimes = testTimes;

    }
    public static double pearson_threshold(double[] a, double[] b, double threshold){
        if (pearson.correlation(a,b)>threshold) return 1.0;
        else return 0.0;
    }


    public static double spearman_threshold(double[] a, double[] b, double threshold){
        if (spearman.correlation(a,b)>threshold) return 1.0;
        else return 0.0;
    }

    public double align_polynomial_spearman(SplineDictionary splineDict, double[] parameters,
                                           double threshold){
        PolynomialFunction poly = new PolynomialFunction(parameters);
        double total = 0.0;
        double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        double[] bounds = {transformedTimes[0],transformedTimes[transformedTimes.length-1],
                (transformedTimes[0]+transformedTimes[transformedTimes.length-1])/2.0};
        if(!checkBounds(bounds)) return Double.NEGATIVE_INFINITY;
        for(String protein:this.proteinList){
            double[] refVals = Doubles.toArray(this.refTable.row(protein).values());
            double[] compVals = Arrays.stream(transformedTimes).map(e->
                    splineDict.getSpline(this.homologDict.get(protein)).value(e)).toArray();
            total+= spearman_threshold(refVals,compVals,threshold);

        }
        return total;
    }

    public double align_polynomial_pearson(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        //PolynomialFunction poly = new PolynomialFunction(parameters);
        double total = 0.0;
        double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //create time surface
        double[] sTimes = new double[50];
        int index =0;
        //make a static variable
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return Double.NEGATIVE_INFINITY;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return Double.NEGATIVE_INFINITY;
        for(String protein:this.proteinList){
            double[] refVals = Doubles.toArray(this.refTable.row(protein).values());
            int[] indices = getBoundIndices(transformedTimes);
            //System.out.println(Arrays.toString(indices));
            refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            transformedTimes = Arrays.copyOfRange(transformedTimes,indices[0],indices[1]);
            //try {
                PolynomialSplineFunction spline = splineDict.getSpline(this.homologDict.get(protein));
                double[] compVals = Arrays.stream(transformedTimes).map(e ->
                        spline.value(e)).toArray();
                total+= pearson_threshold(refVals, compVals, threshold);
            /*}catch(Exception e){
                //System.out.println(Arrays.toString(indices));

                return Double.NEGATIVE_INFINITY;
            }*/



        }
        return total;
    }

    public List<String> align_polynomial_pearson_list(SplineDictionary splineDict, double[] parameters
            , double threshold){
        LinkedList<String>  passList = new LinkedList();
        PolynomialFunction poly = new PolynomialFunction(parameters);
        double total = 0.0;
        double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //create time surface
        double[] sTimes = new double[50];
        int index =0;
        //make a static variable
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return null;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return null;
        for(String protein:this.proteinList){
            double[] refVals = Doubles.toArray(this.refTable.row(protein).values());
            int[] indices = getBoundIndices(transformedTimes);
            //System.out.println(Arrays.toString(indices));
            refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            transformedTimes = Arrays.copyOfRange(transformedTimes,indices[0],indices[1]);
            try {
                PolynomialSplineFunction spline = splineDict.getSpline(this.homologDict.get(protein));
                double[] compVals = Arrays.stream(transformedTimes).map(e ->
                        spline.value(e)).toArray();
                if (pearson_threshold(refVals, compVals, threshold) >0){
                    passList.add(protein);
                }
            }catch(Exception e){
                //System.out.println(Arrays.toString(indices));

                return null;
            }



        }
        return passList;
    }

    public static double extrapolateValue(PolynomialSplineFunction spline, double v){
        double[] knots = spline.getKnots();
        PolynomialFunction[] splines = spline.getPolynomials();
        if(knots[0]>v){
            return splines[0].value(v);
        }
        else if(knots[knots.length-1]<v){
            return splines[splines.length-1].value(v);
        }
        else return spline.value(v);
    }

    public int[] getBoundIndices(double[] transformedTimes){
        int[] bounds = new int[2];
        bounds[0] = Arrays.binarySearch(transformedTimes, this.compTimes[0]);
        if(bounds[0]<0) bounds[0] = -1*(bounds[0]+1);
        //System.out.println(Arrays.binarySearch(transformedTimes, this.compTimes[1]));
        bounds[1] = Math.min(Arrays.binarySearch(transformedTimes, this.compTimes[compTimes.length-1]), compTimes.length);
        if(bounds[1]<0) bounds[1] = -1*(bounds[1]+1);
        return bounds;
    }

    public boolean checkBounds(double[] refTimes){
        if (!(Math.max(0.0, refTimes[0]-compTimes[0])+Math.max(0.0, compTimes[compTimes.length-1]
                -refTimes[refTimes.length-1]) >= (compTimes[0]+compTimes[compTimes.length-1])/2.0)){
                return false;
        }
        int tot = 0;
        for(double time: refTimes){
            if(time>= compTimes[0] &&  time<= compTimes[compTimes.length-1]) tot++;
        }
        return (tot >= (int)((double)refTimes.length)/2.0);

    }
}
