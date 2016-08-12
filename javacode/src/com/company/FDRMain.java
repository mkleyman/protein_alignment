package com.company;

import com.Opitmizers.*;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import com.parsers.HomologParser;
import com.parsers.ProteinExpressionParser;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

public class FDRMain {

    public static void run(String mouse, String human, String out, String homologs, char mode){

        try{

            TreeBasedTable human_table = ProteinExpressionParser.parseFile(human);
            System.out.println(human_table.rowKeySet().size());
            System.out.println(human_table.columnKeySet().toString());
            TreeBasedTable mouse_table = ProteinExpressionParser.parseFile(mouse);
            System.out.println(mouse_table.rowKeySet().size());
            System.out.println(mouse_table.columnKeySet().toString());
            Map<String,Double> informationMapMouse = ProteinExpressionParser.informationMap(mouse);
            Map<String,Double> informationMapHuman = ProteinExpressionParser.informationMap(human);
            //mouse_table = ProteinExpressionParser.normalize(mouse_table);
            //human_table = ProteinExpressionParser.normalize(human_table);
            Map<String,String> homologMap = HomologParser.parse(homologs,"mouse, laboratory");
            homologMap = HomologParser.trim(homologMap, mouse_table.rowKeySet(),human_table.rowKeySet());
            Set<String> human_proteins = new HashSet<>();
            human_proteins.addAll(homologMap.values());
            System.out.println(human_proteins.size());
            SplineDictionary refDict = new SplineDictionary(mouse_table, new AkimaSplineInterpolator(),
                    homologMap.keySet());
            double[] compTimes= Doubles.toArray(human_table.columnKeySet());
            double[] refTimes = Doubles.toArray(mouse_table.columnKeySet());
            double[] checkTimes = new double[50];
            int index = 0;
            for(double num = refTimes[0]; num<=refTimes[refTimes.length-1];
                num+=(refTimes[refTimes.length-1]-refTimes[0])/50.0){
                checkTimes[index]= num;
                index++;

            }
            Set<String> proteinSet = new HashSet<>();
            for(String prot: homologMap.keySet()){
                if(informationMapMouse.get(prot)>=1.0 && informationMapHuman.get(homologMap.get(prot))>= 1.0){
                    proteinSet.add(prot);

                }
            }
            //double[] compBounds =  {compTimes[0], compTimes[compTimes.length-1]};
            Aligner aligner = new Aligner(compTimes,proteinSet, homologMap, refTimes,mouse_table,
                    checkTimes,refDict);

            LinkedHashMap<File, Optimizer> optMap = new LinkedHashMap<File, Optimizer>();

            //File quadFile = new File(out+"/quadratic.csv");
            //quadFile.createNewFile();
            //optMap.put(quadFile,new QuadraticOptimizer());

            File sqrtFile = new File(out+"/sqrt.csv");
            sqrtFile.createNewFile();
            optMap.put(sqrtFile, new SqrtOptimizer(mode));

            File expFile = new File(out+"/exp.csv");
            expFile.createNewFile();
            optMap.put(expFile, new ExponentialOptimizer(mode));

            File logFile = new File(out+"/logarithm.csv");
            logFile.createNewFile();
            optMap.put(logFile, new LogarithmOptimizer(mode));

            File linear = new File(out+"/linear.csv");
            linear.createNewFile();
            optMap.put(linear, new LinearOptimizer(mode));



            FDRCalculator.calculateAllFDR(optMap,new double[]{0.7},aligner,human_table,
                    new AkimaSplineInterpolator(), human_proteins, 500, mouse_table,homologMap);





        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
