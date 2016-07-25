package com.company;

import com.Opitmizers.*;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import com.parsers.HomologParser;
import com.parsers.ProteinExpressionParser;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

public class ParamMain {

    public static void run(String mouse, String human, String result_file, String homologs){

        try{

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
            LinkedList<Optimizer> optList = new LinkedList<>();
            optList.add(new LinearOptimizer());
            optList.add(new LogarithmOptimizer());
            optList.add(new SqrtOptimizer());
            //optList.add(new LogitOptimizer());
            optList.add(new QuadraticOptimizer());
            //optList.add(new CubicOptmizer());
            optList.add(new ExponentialOptimizer());
            ParameterSearch parSearch = new ParameterSearch(aligner, spDict, new double[]{0.6,0.7,0.8},optList);
            parSearch.recordResults(result_file);
           // Aligner.setCompTimes(compBounds);
            //LinearOptimizer linear = new LinearOptimizer();

           // linear.optimizePearson(aligner, spDict, 0.8);




        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
