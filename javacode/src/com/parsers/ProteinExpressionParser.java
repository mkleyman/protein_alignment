package com.parsers;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by mkleyman on 7/15/2016.
 */
public class ProteinExpressionParser {
    /**
     *
     * @param filename
     * @return
     * @throws FileNotFoundException
     */
    public static TreeBasedTable<String,Double,Double> parseFile(String filename) throws FileNotFoundException{
        File file = new File(filename);
        Scanner scan = new Scanner(file);
        HashBasedTable<String,Double,LinkedList<Double>> proteinExpression = HashBasedTable.create();
        TreeBasedTable<String,Double,Double> proteinExpressionAvg = TreeBasedTable.create();
        List<String> proteins = new LinkedList<>();
        String[] header = scan.nextLine().split(",");
        double[] times = new double[header.length-1];
        //read in times for column names
        for(int i=1;i<header.length;i++){
            times[i-1] = Double.parseDouble(header[i]);
        }
        //read in protein expression
        while(scan.hasNextLine()){
            String[] row = scan.nextLine().split(",");
            proteins.add(row[0]);
            for(int i=1;i<row.length;i++){
                if(!proteinExpression.contains(row[0],times[i-1])){
                    proteinExpression.put(row[0],times[i-1],new LinkedList());
                }
                proteinExpression.get(row[0],times[i-1]).add(Double.parseDouble(row[i]));
            }
        }
        scan.close();
        double[] uniqueTimes = Arrays.stream(times).distinct().toArray();
        //take median of table values
        Median med = new Median();
        for(String protein: proteins){
            for(double time: uniqueTimes){
                proteinExpressionAvg.put(protein, time,
                        med.evaluate(Doubles.toArray(proteinExpression.get(protein, time))));
            }
        }

        return proteinExpressionAvg;

    }


    public static TreeBasedTable<String,Double,Double> parseFileClean(String filename) throws FileNotFoundException{
        File file = new File(filename);
        Scanner scan = new Scanner(file);
        HashBasedTable<String,Double,LinkedList<Double>> proteinExpression = HashBasedTable.create();
        TreeBasedTable<String,Double,Double> proteinExpressionAvg = TreeBasedTable.create();
        List<String> proteins = new LinkedList<>();
        String[] header = scan.nextLine().split(",");
        double[] times = new double[header.length-1];
        //read in times for column names
        for(int i=1;i<header.length;i++){
            times[i-1] = Double.parseDouble(header[i]);
        }
        //read in protein expression
        while(scan.hasNextLine()){
            String[] row = scan.nextLine().split(",");
            proteins.add(row[0]);
            for(int i=1;i<row.length;i++){
                if(!proteinExpression.contains(row[0],times[i-1])){
                    proteinExpression.put(row[0],times[i-1],new LinkedList());
                }
                proteinExpression.get(row[0],times[i-1]).add(Double.parseDouble(row[i]));
            }
        }
        scan.close();
        double[] uniqueTimes = Arrays.stream(times).distinct().toArray();
        //take median of table values
        Median med = new Median();
        for(String protein: proteins){
            if(cleanCheck( proteinExpression.rowMap().get(protein))) {
                for (double time : uniqueTimes) {
                    proteinExpressionAvg.put(protein, time,
                            med.evaluate(Doubles.toArray(proteinExpression.get(protein, time))));
                }
            }
        }

        return proteinExpressionAvg;

    }

    public static boolean cleanCheck(Map<Double,LinkedList<Double>> proteinValMap){
        Median med = new Median();
        int size = proteinValMap.values().size();
        double[] vars = new double[size];
        double[] meds = new double[size];

        int i = 0;
        for(LinkedList<Double> singleTime: proteinValMap.values()){
            meds[i] = med.evaluate(Doubles.toArray(singleTime));
            vars[i] = StatUtils.variance(Doubles.toArray(singleTime));
            i++;
        }
        return (StatUtils.variance(meds)>StatUtils.mean(vars));
    }

    public static TreeBasedTable<String,Double,Double> normalize(TreeBasedTable<String,Double,Double> protTable){
        for(Double time: protTable.columnKeySet()){
            double avg = StatUtils.mean(Doubles.toArray(protTable.columnMap().get(time).values()));
            for(String prot: protTable.rowKeySet()){
                double val = protTable.get(prot,time)/avg;
                protTable.put(prot,time,val);
            }
        }
        return protTable;
    }
}


