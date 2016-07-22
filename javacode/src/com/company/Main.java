package com.company;

import com.Opitmizers.*;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import com.parsers.HomologParser;
import com.parsers.ProteinExpressionParser;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

public class Main {

    public static void main(String[] args){

        try{
            String human = "../processed/human_proteins.csv";
            String mouse = "../processed/mouse_proteins.csv";
            //String result_file = "../processed/java_results.csv";
            String homologs = "../raw_data/HOM_MouseHumanSequence.rpt";
            TreeBasedTable human_table = ProteinExpressionParser.parseFile(human);
            System.out.println(human_table.rowKeySet().size());
            System.out.println(human_table.columnKeySet().toString());
            TreeBasedTable mouse_table = ProteinExpressionParser.parseFile(mouse);
            System.out.println(mouse_table.rowKeySet().size());
            System.out.println(mouse_table.columnKeySet().toString());
            Map<String,String> homologMap = HomologParser.parse(homologs,"mouse, laboratory");
            homologMap = HomologParser.trim(homologMap, mouse_table.rowKeySet(),human_table.rowKeySet());
            Set<String> human_proteins = new HashSet<>();
            human_proteins.addAll(homologMap.values());
            System.out.println(human_proteins.size());
            SplineDictionary spDict = new SplineDictionary(human_table,new AkimaSplineInterpolator(),human_proteins);
            double[] compTimes= Doubles.toArray(human_table.columnKeySet());
            double[] refTimes = Doubles.toArray(mouse_table.columnKeySet());
            double[] checkTimes = new double[50];
            int index = 0;
            for(double num = refTimes[0]; num<=refTimes[refTimes.length-1];
                num+=(refTimes[refTimes.length-1]-refTimes[0])/50.0){
                checkTimes[index]= num;
                index++;

            }
            //double[] compBounds =  {compTimes[0], compTimes[compTimes.length-1]};
            Aligner aligner = new Aligner(compTimes,homologMap.keySet(), homologMap, refTimes,mouse_table,
                    checkTimes);
           // Aligner.setCompTimes(compBounds);
            //LinearOptimizer linear = new LinearOptimizer();
            LinkedHashMap<File, Optimizer> optMap = new LinkedHashMap<File, Optimizer>();
            File sqrtFile = new File("../processed/sqrt.csv");
            optMap.put(sqrtFile, new SqrtOptimizer());
            File quadFile = new File("../processed/quadratic.csv");
            optMap.put(quadFile,new QuadraticOptimizer());
            File logitFile = new File("../processed/logit.csv");
            optMap.put(logitFile,new LogitOptimizer());
            File logFile = new File("../processed/logarithm.csv");
            optMap.put(logFile, new LogarithmOptimizer());
            File linear = new File("../processed/linear.csv");
            optMap.put(linear, new LinearOptimizer());
            File expFile = new File("../processed/exp.csv");
            optMap.put(expFile, new ExponentialOptimizer());
            File cubicFile = new File("../processed/cubic.csv");
            optMap.put(cubicFile,new CubicOptmizer());
            FDRCalculator.calculateAllFDR(optMap,new double[]{0.6,0.7,0.8},aligner,human_table,
                    new AkimaSplineInterpolator(), human_proteins, 500);

           // linear.optimizePearson(aligner, spDict, 0.8);




        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
