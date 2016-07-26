package com.parsers;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
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
}


