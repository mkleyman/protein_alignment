package com.company;

import com.Opitmizers.*;
import com.TransformationFunctions.ExpTransform;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import com.parsers.HomologParser;
import com.parsers.ProteinExpressionParser;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Exp;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import java.util.*;

public class TestMain {

    public static void main(String[] args){

        try{
            String human = "../processed/human_proteins.csv";
            String mouse = "../processed/mouse_proteins.csv";
            String result_file = "../processed/java_results.csv";
            String homologs = "../raw_data/HOM_MouseHumanSequence.rpt";
            TreeBasedTable<String,Double,Double> human_table = ProteinExpressionParser.parseFile(human);

            System.out.println(human_table.rowKeySet().size());
            System.out.println(human_table.columnKeySet().toString());
            TreeBasedTable<String,Double,Double> mouse_table = ProteinExpressionParser.parseFile(mouse);
            //mouse_table = ProteinExpressionParser.normalize(mouse_table);
            //human_table = ProteinExpressionParser.normalize(human_table);
            Map<String,Double> informationMapMouse = ProteinExpressionParser.informationMap(mouse);
            Map<String,Double> informationMapHuman = ProteinExpressionParser.informationMap(human);
            System.out.println(mouse_table.rowKeySet().size());
            System.out.println(mouse_table.columnKeySet().toString());
            Map<String,String> homologMap = HomologParser.parse(homologs,"mouse, laboratory");
            homologMap = HomologParser.trim(homologMap, mouse_table.rowKeySet(),human_table.rowKeySet());
            Set<String> human_proteins = new HashSet<>();
            human_proteins.addAll(homologMap.values());
            System.out.println(human_proteins.size());

            SplineDictionary spDict = new SplineDictionary(human_table,new AkimaSplineInterpolator(),
                    human_proteins);

            spDict.addTable(mouse_table,homologMap.keySet());

            SplineDictionary spDictRand = new SplineDictionary(human_table,new AkimaSplineInterpolator(),
                    human_proteins,117);
            spDictRand.addTable(mouse_table,homologMap.keySet(),3435 );

            SplineDictionary refDict = new SplineDictionary(mouse_table, new AkimaSplineInterpolator(),
                    homologMap.keySet());
            double[] compTimes= Doubles.toArray(human_table.columnKeySet());
            double[] refTimes = Doubles.toArray(mouse_table.columnKeySet());
            double[] checkTimes = new double[50];
            int index = 0;
            for(double num = refTimes[0]; num<=refTimes[refTimes.length-1];
                num+=(refTimes[refTimes.length-1]-refTimes[0])/50){
                checkTimes[index]= num;
                index++;

            }
            //double[] compBounds =  {compTimes[0], compTimes[compTimes.length-1]};
            Aligner aligner = new Aligner(compTimes,homologMap.keySet(), homologMap, refTimes,mouse_table,
                    checkTimes,refDict,informationMapMouse,informationMapHuman);


            String line;
            //ParamMain.run(mouse,human,result_file,homologs,'s');
            /*
            for(String protein: homologMap.keySet()){
                if(informationMapMouse.get(protein)>1.0 &&
                        informationMapHuman.get(homologMap.get(protein))>1.0)
                System.out.println(protein);
            }*/
            /*
            for(String protein : homologMap.keySet()){
               if(homologMap.values().contains(protein)) System.out.println(protein);
            }*/


            Optimizer opt = new LinearOptimizer('s');
            double[] randresult= opt.optimizePearson(aligner, spDictRand, 0.7);
            System.out.println(Doubles.join(",", randresult));

            double[] result= opt.optimizePearson(aligner, spDict, 0.7);
            System.out.println(Doubles.join(",", result));


            /*
            UnivariateFunction lin = new PolynomialFunction(new double[]{-1080.0,93.12560000000015});
            Set<String> protSet = aligner.align_pearson_set(spDict,lin,0.6);*/
            /*
            for(String pro:protSet){
                System.out.println(pro);
            }*/




            /*

            List<String> shuffled = new LinkedList();
            shuffled.addAll(homologMap.keySet());
            Collections.shuffle(shuffled);
            int i = 0;

            for(String chosen: shuffled){
                if(i==5) break;
                i++;
                for(Double mTime: mouse_table.columnKeySet()){
                    line = "real,"+chosen+","+mTime+","+mouse_table.get(chosen,mTime);
                    System.out.println(line);
                }
                for(Double sTime:checkTimes){
                    line = "spline,"+chosen+","+sTime+","+refDict.getSpline(chosen).evaluate(sTime);
                    System.out.println(line);
                }

            }*/


            UnivariateFunction lin = new PolynomialFunction(new double[]{-856.0,95.62360000000015});
            List<String> chosen_proteins= aligner.align_polynomial_sp_list(spDict,lin,0.7);
            //List<String> chosen =
            //};
            //Set<String> chosen_proteins =aligner.align_pearson_refset(spDict,lin,0.7);
            int num_chosen = 0;


            for(String chosen:chosen_proteins){
                if(num_chosen==5) break;
                for(Double mTime: checkTimes){
                    //line = chosen +","+homologMap.get(chosen)+",mouse,"+mTime+","+mouse_table.get(chosen,mTime);
                    //System.out.println(line);
                    line = chosen +","+homologMap.get(chosen)+",mouse_spline,"+mTime+","+spDict.getSpline(chosen).evaluate(mTime);

                    System.out.println(line);
                    //if(lin.value(mTime)>=1.0 && lin.value(mTime)<=3285) {
                        line = chosen +","+homologMap.get(chosen)+ ",mouse_transformed," + lin.value(mTime) + "," + spDict.getSpline(chosen).evaluate(mTime);
                        System.out.println(line);
                        line = chosen +","+homologMap.get(chosen)+ ",human_spline," + lin.value(mTime) + "," +
                                spDict.getSpline(homologMap.get(chosen)).evaluate(lin.value(mTime));
                        System.out.println(line);
                   // }
                }
                for(Double hTime: human_table.columnKeySet()){
                    line = chosen +","+homologMap.get(chosen)+",human,"+hTime+","+human_table.get(homologMap.get(chosen),hTime);
                    System.out.println(line);
                }
                num_chosen++;

            }
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
