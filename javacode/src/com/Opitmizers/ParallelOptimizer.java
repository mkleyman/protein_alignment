package com.Opitmizers;

import com.company.Aligner;
import com.company.SplineDictionary;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import java.io.File;
import java.nio.charset.Charset;

/**
 * Created by mkleyman on 7/19/2016.
 */
public class ParallelOptimizer implements Runnable {
    private Optimizer opt;
    private Aligner aligner;
    private SplineDictionary splineDict;
    private double threshold;
    private File outFile;
    public ParallelOptimizer(Optimizer opt, Aligner aligner, SplineDictionary splineDict,
                             double threshold, File outFile){
        this.opt = opt;
        this.aligner = aligner;
        this.splineDict = splineDict;
        this.threshold = threshold;
        this.outFile = outFile;

    }

    public void run(){
        double[] result = opt.optimizePearson(this.aligner, this.splineDict, this.threshold);
        synchronized (this.outFile){
            try {
                Files.append(Doubles.join(",", result)+"\n", this.outFile, Charset.defaultCharset());
            }catch(Exception e){
                e.printStackTrace();
            }
        }
    }

}
