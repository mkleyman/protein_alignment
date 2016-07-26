package com.company;

import com.Opitmizers.Optimizer;
import com.Opitmizers.QuadraticOptimizer;
import com.TransformationFunctions.ExpTransform;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import com.parsers.GOParser;
import com.parsers.HomologParser;
import com.parsers.ProteinExpressionParser;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class GoMain {

    public static void main(String[] args){

        try{
            String human = "../processed/human_proteins.csv";
            String mouse = "../processed/mouse_proteins.csv";
            String goTable = "../processed/go_table_exp.csv";
            String homologs = "../raw_data/HOM_MouseHumanSequence.rpt";
            String goFile = "../raw_data/HUMAN_9606_idmapping_selected.tab";
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

            Multimap<String,String> goMap = GOParser.parseGO(goFile, human_proteins);
            System.out.println("go parsed in");

            //double[] compBounds =  {compTimes[0], compTimes[compTimes.length-1]};
            Aligner aligner = new Aligner(compTimes,homologMap.keySet(), homologMap, refTimes,mouse_table,
                    checkTimes);
            UnivariateFunction funct = new ExpTransform(-240.0,0.17240000000000208);

            GOTools.createGOTable(aligner,spDict,funct,.7,goMap,goTable);
            System.out.println("done");



           // Aligner.setCompTimes(compBounds);
            //LinearOptimizer linear = new LinearOptimizer();

           // linear.optimizePearson(aligner, spDict, 0.8);




        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
