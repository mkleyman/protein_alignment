package com.company;

import JavaMI.MutualInformation;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.geometry.euclidean.oned.Vector1D;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umontreal.ssj.functions.MathFunction;


import java.util.*;

/**
 * Created by mkleyman on 7/18/2016.
 */
public class Aligner {
    private final double[] compTimes;
    public final Set<String> proteinList;
    public final Map<String,String> homologDict;
    private final double[] times;
    private final double[] testTimes;
    public final Table<String,Double,Double> refTable;
    private static PearsonsCorrelation pearson = new PearsonsCorrelation();
    private static SpearmansCorrelation spearman = new SpearmansCorrelation();
    private SplineDictionary refDict;
    private Map<String, Double> infoMapRef;
    private Map<String, Double> infoMapComp;




    public Aligner(double[] compTimes, Set<String> proteinList, Map<String,String> homologDict,
                   double[] times, Table<String, Double, Double> refTable, double[] testTimes,
                   SplineDictionary refDict){
        this.compTimes= compTimes;
        this.proteinList = proteinList;
        this.homologDict = homologDict;
        this.times = times;
        this.refTable = refTable;
        this.testTimes = testTimes;
        this.refDict = refDict;

    }

    public Aligner(double[] compTimes, Set<String> proteinList, Map<String,String> homologDict,
                   double[] times, Table<String, Double, Double> refTable, double[] testTimes,
                   SplineDictionary refDict, Map<String,Double> informationMapRef,
                   Map<String,Double> informationMapComp ){
        this.compTimes= compTimes;
        this.proteinList = proteinList;
        this.homologDict = homologDict;
        this.times = times;
        this.refTable = refTable;
        this.testTimes = testTimes;
        this.refDict = refDict;
        this.infoMapRef = informationMapRef;
        this.infoMapComp = informationMapComp;

    }

    public Aligner(double[] compTimes, Set<String> proteinList, Map<String,String> homologDict,
                   double[] times, Table<String, Double, Double> refTable, double[] testTimes){
        this.compTimes= compTimes;
        this.proteinList = proteinList;
        this.homologDict = homologDict;
        this.times = times;
        this.refTable = refTable;
        this.testTimes = testTimes;


    }


    public Aligner(double[] compTimes, Set<String> proteinList, Map<String,String> homologDict,
                   double[] times, Table<String, Double, Double> refTable, double[] testTimes,
                   Map<String,Double> informationMapRef,
                   Map<String,Double> informationMapComp ){
        this.compTimes= compTimes;
        this.proteinList = proteinList;
        this.homologDict = homologDict;
        this.times = times;
        this.refTable = refTable;
        this.testTimes = testTimes;
        this.infoMapRef = informationMapRef;
        this.infoMapComp = informationMapComp;

    }



    public static double largestSequence(double[] seq, int size){
        double max = StatUtils.sum(Arrays.copyOfRange(seq, 0, size));
        double current = StatUtils.sum(Arrays.copyOfRange(seq,0,size));
        for(int i=size;i<seq.length;i++){
            current = current+seq[i]-seq[i-size];
            if(current>max) max = current;
        }
        return max;
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
                    splineDict.getSpline(this.homologDict.get(protein)).evaluate(e)).toArray();
            total+= spearman_threshold(refVals,compVals,threshold);

        }
        return total;
    }

    public double align_polynomial_pearson(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        double[] trTimes;
        int[] indices;
        double total = 0.0;

        //create time surface
        double[] sTimes = new double[testTimes.length];
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return Double.NEGATIVE_INFINITY;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return Double.NEGATIVE_INFINITY;
        double[] mTimes;
        if(refDict==null){
            double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
            mTimes = Arrays.copyOfRange(this.times, indices[0], indices[1]);

        }
        else {
            indices = getBoundIndices(sTimes);
            trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
            //double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            //indices = getBoundIndices(transformedTimes);
            //trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
            mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);
        }
        double[] refVals;
        for(String protein:this.proteinList){

            if(refDict==null) {
                refVals = Doubles.toArray(this.refTable.row(protein).values());
                System.out.println("hello");
                refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            }
            else {

                MathFunction refSpline = this.refDict.getSpline(protein);
                refVals = Arrays.stream(mTimes).map(e ->
                        refSpline.evaluate(e)).toArray();
            }
            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            if(infoMapRef != null) {
                double info = 0.0;
                if(infoMapComp.get(homologDict.get(protein))>1.0 && infoMapRef.get(protein)>1.0) info = 1.0;

                total += pearson_threshold(refVals, compVals,threshold)*info;
                /*
                double corr = pearson.correlation(refVals, compVals);
                if(!Double.isNaN(corr)) {
                    total += pearson.correlation(refVals, compVals) * ((infoMapRef.get(protein) + infoMapComp.get(homologDict.get(protein))) / 2.0);
                }*/

            }else{

                total += pearson_threshold(refVals, compVals, threshold);
            }
        }
        return total;
    }


    public double align_polynomial_sp(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        double[] trTimes;
        int[] indices;
        double total = 0.0;

        //create time surface
        double[] sTimes = new double[testTimes.length];
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return Double.NEGATIVE_INFINITY;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return Double.NEGATIVE_INFINITY;
        double[] mTimes;

        indices = getBoundIndices(sTimes);
        trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
        //double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //indices = getBoundIndices(transformedTimes);
        //trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
        mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);

        double[] refVals;
        for(String protein:this.proteinList){


            MathFunction refSpline = splineDict.getSpline(protein);
            refVals = Arrays.stream(mTimes).map(e ->
                    refSpline.evaluate(e)).toArray();

            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            if(infoMapRef != null) {
                double info = 0.0;
                if(infoMapComp.get(homologDict.get(protein))>1.0 && infoMapRef.get(protein)>1.0) info = 1.0;

                total += pearson_threshold(refVals, compVals,threshold)*info;
                /*
                double corr = pearson.correlation(refVals, compVals);
                if(!Double.isNaN(corr)) {
                    total += pearson.correlation(refVals, compVals) * ((infoMapRef.get(protein) + infoMapComp.get(homologDict.get(protein))) / 2.0);
                }*/

            }else{

                total += pearson_threshold(refVals, compVals, threshold);
            }
        }
        return total;
    }



    public List<String> align_polynomial_sp_list(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        double[] trTimes;
        int[] indices;
        double total = 0.0;

        //create time surface
        double[] sTimes = new double[testTimes.length];
        LinkedList<String> passed = new LinkedList<>();
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return null;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return null;
        double[] mTimes;

        indices = getBoundIndices(sTimes);
        trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
        //double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //indices = getBoundIndices(transformedTimes);
        //trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
        mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);

        double[] refVals;
        for(String protein:this.proteinList){


            MathFunction refSpline = splineDict.getSpline(protein);
            refVals = Arrays.stream(mTimes).map(e ->
                    refSpline.evaluate(e)).toArray();

            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            if(infoMapRef != null) {
                double info = 0.0;
                if(infoMapComp.get(homologDict.get(protein))>1.0 && infoMapRef.get(protein)>1.0) info = 1.0;

                if(pearson_threshold(refVals, compVals,threshold)> 0.5 && info>0.5) passed.add(protein);
                /*
                double corr = pearson.correlation(refVals, compVals);
                if(!Double.isNaN(corr)) {
                    total += pearson.correlation(refVals, compVals) * ((infoMapRef.get(protein) + infoMapComp.get(homologDict.get(protein))) / 2.0);
                }*/

            }else{

                total += pearson_threshold(refVals, compVals, threshold);
            }
        }
        return passed;
    }

    public double align_polynomial_info(SplineDictionary splineDict, UnivariateFunction poly){
        double[] trTimes;
        int[] indices;
        double total = 0.0;

        //create time surface
        double[] sTimes = new double[testTimes.length];
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return Double.NEGATIVE_INFINITY;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return Double.NEGATIVE_INFINITY;
        double[] mTimes;
        if(refDict==null){
            double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
            mTimes = Arrays.copyOfRange(this.times, indices[0], indices[1]);

        }
        else {
            indices = getBoundIndices(sTimes);
            trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
            //double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            //indices = getBoundIndices(transformedTimes);
            //trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
            mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);
        }
        double[] refVals;
        double[] proteinInformation = new double[this.proteinList.size()];
        double[] proteinCorrelation = new double[this.proteinList.size()];
        int protIndex = 0;
        for(String protein:this.proteinList){

            if(refDict==null) {
                refVals = Doubles.toArray(this.refTable.row(protein).values());
                System.out.println("hello");
                refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            }
            else {

                MathFunction refSpline = this.refDict.getSpline(protein);
                refVals = Arrays.stream(mTimes).map(e ->
                        refSpline.evaluate(e)).toArray();
            }
            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            proteinCorrelation[protIndex] = pearson.correlation(refVals,compVals);
            if(Double.isNaN(proteinCorrelation[protIndex])){
                proteinCorrelation[protIndex] =0;
            }
            proteinInformation[protIndex] = (infoMapRef.get(protein)+
                    infoMapComp.get(homologDict.get(protein)))/2.0;
            protIndex++;
        }
        return spearman.correlation(proteinCorrelation,proteinInformation);
    }

    public double align_polynomial_difference(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        double[] trTimes;
        int[] indices;
        double total = 0.0;

        //create time surface
        double[] sTimes = new double[testTimes.length];
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return Double.NEGATIVE_INFINITY;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return Double.NEGATIVE_INFINITY;

        if(refDict==null){
            double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);

        }
        else {
            //indices = getBoundIndices(sTimes);
            //trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
            double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
        }
        double[] refVals;
        double[] mTimes = Arrays.copyOfRange(this.times, indices[0], indices[1]);
        for(String protein:this.proteinList){

            if(refDict==null) {
                refVals = Doubles.toArray(this.refTable.row(protein).values());
                System.out.println("hello");
                refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            }
            else {

                MathFunction refSpline = this.refDict.getSpline(protein);
                refVals = Arrays.stream(mTimes).map(e ->
                        refSpline.evaluate(e)).toArray();
            }
            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            if(infoMapRef != null) {
                total += difference(refVals, compVals)*((infoMapRef.get(protein)+
                        infoMapComp.get(homologDict.get(protein)))/2.0);
            }else{
                total += difference(refVals, compVals);
            }
        }
        return (-1.0*total);
    }


    public double align_polynomial_pearson_all(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){


        double total = 0.0;
        double[] transformedTimes;
        // = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //create time surface
        double[] sTimes = new double[testTimes.length];
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
        int[] indices ;
        double[] refVals;
        double[] compVals;
        double[] trTimes;
        double[] mTimes;
        if(refDict==null){
            transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
            mTimes = Arrays.copyOfRange(this.times, indices[0], indices[1]);

        }
        else {
            indices = getBoundIndices(sTimes);
            trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
            //double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
            //indices = getBoundIndices(transformedTimes);
            //trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
            mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);
        }

        //double[] trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);
        //double[] refTimes = Arrays.copyOfRange(this.times, indices[0], indices[1]);
        //double[] tTimes = Arrays.copyOfRange(transformedTimes,indices[0], indices[1]);
        double[] allCorr = new double[mTimes.length];
        List<String> proteinSortedList = new LinkedList<>();
        proteinSortedList.addAll(this.proteinList);
        int corrIndex = 0;


        for(double time:mTimes) {
            int i = 0;

            //refVals = new double[this.proteinList.size()];
            compVals = new double[this.proteinList.size()];

            refVals = new double[this.proteinList.size()];

            for (String protein : proteinSortedList) {
                if(refDict==null) {
                    refVals[i] = refTable.get(protein, time);
                }
                else{
                    refVals[i] = refDict.getSpline(protein).evaluate(time);
                }

                compVals[i] = splineDict.getSpline(homologDict.get(protein)).evaluate(poly.value(trTimes[corrIndex]));
                i++;
            }
            allCorr[corrIndex] = pearson.correlation(refVals,compVals);
            corrIndex++;
        }
        return StatUtils.sum(allCorr);

    }






    public Set<String> align_pearson_set(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        double[] trTimes;
        int[] indices;
        Set<String> chosen_proteins = new HashSet<>();

        double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //create time surface
        double[] sTimes = new double[testTimes.length];
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return null;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return null;

        if(refDict==null){
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);

        }
        else {
            indices = getBoundIndices(sTimes);
            trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
        }
        double[] refVals;
        for(String protein:this.proteinList){
            double[] mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);
            if(refDict==null) {
                refVals = Doubles.toArray(this.refTable.row(protein).values());
                refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            }
            else {

                MathFunction refSpline = this.refDict.getSpline(protein);
                refVals = Arrays.stream(mTimes).map(e ->
                        refSpline.evaluate(e)).toArray();
            }
            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            if(pearson_threshold(refVals, compVals, threshold)>0.5){
                chosen_proteins.add(this.homologDict.get(protein));
            }
        }
        return chosen_proteins;
    }

    public Set<String> align_pearson_refset(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        double[] trTimes;
        int[] indices;
        Set<String> chosen_proteins = new HashSet<>();

        double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //create time surface
        double[] sTimes = new double[testTimes.length];
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return null;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return null;

        if(refDict==null){
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);

        }
        else {
            indices = getBoundIndices(sTimes);
            trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
        }
        double[] refVals;
        for(String protein:this.proteinList){
            double[] mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);
            if(refDict==null) {
                refVals = Doubles.toArray(this.refTable.row(protein).values());
                refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            }
            else {

                MathFunction refSpline = this.refDict.getSpline(protein);
                refVals = Arrays.stream(mTimes).map(e ->
                        refSpline.evaluate(e)).toArray();
            }
            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            if(pearson_threshold(refVals, compVals, threshold)>0.5){
                chosen_proteins.add(protein);
            }
        }
        return chosen_proteins;
    }


    public Map<Double,String> align_polynomial_refDifMap(SplineDictionary splineDict, UnivariateFunction poly
            , double threshold){
        Map<Double,String> refMap =  new HashMap();
        double[] trTimes;
        int[] indices;
        double total = 0.0;
        double[] transformedTimes = Arrays.stream(this.times).map(e->poly.value(e)).toArray();
        //create time surface
        double[] sTimes = new double[testTimes.length];
        int index =0;
        for(double sTime: this.testTimes){
            sTimes[index] = poly.value(sTime);
            if (!(index==0) && sTimes[index]<sTimes[index-1]){
                return null;
            }
            index++;
        }
        if(!checkBounds(sTimes)) return null;

        if(refDict==null){
            indices = getBoundIndices(transformedTimes);
            trTimes = Arrays.copyOfRange(transformedTimes, indices[0], indices[1]);

        }
        else {
            indices = getBoundIndices(sTimes);
            trTimes = Arrays.copyOfRange(sTimes, indices[0], indices[1]);
        }
        double[] refVals;
        for(String protein:this.proteinList){
            double[] mTimes = Arrays.copyOfRange(this.testTimes, indices[0], indices[1]);
            if(refDict==null) {
                refVals = Doubles.toArray(this.refTable.row(protein).values());
                refVals = Arrays.copyOfRange(refVals,indices[0],indices[1]);
            }
            else {

                MathFunction refSpline = this.refDict.getSpline(protein);
                refVals = Arrays.stream(mTimes).map(e ->
                        refSpline.evaluate(e)).toArray();
            }
            MathFunction spline = splineDict.getSpline(this.homologDict.get(protein));
            double[] compVals = Arrays.stream(trTimes).map(e ->
                    spline.evaluate(e)).toArray();

            refMap.put(difference(refVals, compVals), protein);
        }
        return refMap;
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
        if ((Math.min(refTimes[refTimes.length-1], compTimes[compTimes.length-1]
                -Math.max(compTimes[0], refTimes[0])) < ((compTimes[compTimes.length-1]-compTimes[0])*0.75))){
                return false;
        }
        int tot = 0;
        for(double time: refTimes){
            if(time>= compTimes[0] &&  time<= compTimes[compTimes.length-1]) tot++;
        }
        //System.out.println(tot);
        //System.out.println(((double)refTimes.length)*0.75);

        return (tot >= ((((double)refTimes.length)*0.75)));

    }

    public double difference(double[] a, double[] b){
        double total = 0;
        for(int i=0;i<a.length;i++){
            total+=Math.abs(a[i]-b[i]);
        }
        return total;
    }
}
