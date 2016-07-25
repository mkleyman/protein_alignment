package com.company;

import com.Opitmizers.Optimizer;
import com.google.common.primitives.Doubles;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

/**
 * Created by mkleyman on 7/20/2016.
 */
public class ParameterSearch {
    private Aligner aligner;
    private SplineDictionary spDict;
    private double[] thresholds;
    private List<Optimizer> optimizerList;

    public ParameterSearch(Aligner aligner, SplineDictionary splines, double[] threshs,
                           List<Optimizer> optimizers){
        this.aligner = aligner;
        this.spDict = splines;
        this.thresholds = threshs;
        this.optimizerList = optimizers;
    }

    public void recordResults(String fileName) throws IOException {
        File outFile = new File(fileName);
        FileWriter writer = new FileWriter(outFile);
        for(Double threshold:thresholds){
            for(Optimizer opt: optimizerList){
                System.out.println(opt.getName()+","+Double.toString(threshold));
                double[] result = opt.optimizePearson(this.aligner,this.spDict,threshold);
                String line = opt.getName()+","+ Doubles.join(",",result)+"\n";
                System.out.println(line);
                writer.write(line);
            }
        }
        writer.close();

    }
}
