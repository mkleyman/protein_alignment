package com.parsers;

import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.collect.TreeBasedTable;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Map;

/**
 * Created by mkleyman on 8/8/2016.
 */
public class SplineWriter {
    public static void writeFile(String filename,  TreeBasedTable<String,Double,Double> human_table,
                          TreeBasedTable<String,Double,Double> mouse_table, SplineDictionary splineDictionary,
                          double[] sTimes, Map<String,String> homologDict,
                          UnivariateFunction funct, Aligner aligner) throws IOException{
        PrintWriter pw = new PrintWriter(new File(filename));
        String humanProtein;
        String line;
        double time_transformed;
        line = "mouse,human,type,time,expression\n";
        pw.print(line);
        //double[] transformedTimes =  Arrays.stream(mouse_table.columnKeySet().toArray()).map(e->funct.value((double)e)).toArray();
        for(String mouseProtein:homologDict.keySet()){
            humanProtein = homologDict.get(mouseProtein);
            for(Double time: mouse_table.columnKeySet()){
                time_transformed = funct.value(time);
                line = mouseProtein+","+humanProtein+",mouse,"+time+","+mouse_table.get(mouseProtein,time)+"\n";
                pw.print(line);
                line = mouseProtein+","+humanProtein+",mouse_transformed,"+time_transformed+","+mouse_table.get(mouseProtein,time)+"\n";
                pw.print(line);
            }
            for(Double time: human_table.columnKeySet()){
                line = mouseProtein+","+humanProtein+",human,"+time+","+human_table.get(humanProtein,time)+"\n";
                pw.print(line);
            }

            for(double time: sTimes){
                line = mouseProtein+","+humanProtein+",mouse_spline,"+time+","+
                        splineDictionary.getSpline(mouseProtein).evaluate(time)+"\n";
                pw.print(line);
                time_transformed = funct.value(time);
                try {
                    line = mouseProtein + "," + humanProtein + ",mouse_spline_transformed," + time_transformed + "," +
                            splineDictionary.getSpline(mouseProtein).evaluate(time) + "\n";
                    pw.print(line);

                    line = mouseProtein + "," + humanProtein + ",human_spline," + time_transformed + "," +
                            splineDictionary.getSpline(humanProtein).evaluate(time_transformed) + "\n";
                    pw.print(line);
                }catch(Exception e){}

            }

        }
        pw.close();
    }
}
