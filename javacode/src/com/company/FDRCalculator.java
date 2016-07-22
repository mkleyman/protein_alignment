package com.company;

import com.Opitmizers.Optimizer;
import com.Opitmizers.ParallelOptimizer;
import com.google.common.collect.TreeBasedTable;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;

import java.io.File;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by mkleyman on 7/19/2016.
 */
public class FDRCalculator {
    public static void calculateAllFDR(Map<File,Optimizer> optimizerMap, double[] thresholds, Aligner aligner,
                                       TreeBasedTable<String,Double,Double> expressionTable, AkimaSplineInterpolator splineMaker,
                                       Set<String> proteinSet, int numTrials){
        ExecutorService executor = Executors.newScheduledThreadPool(Runtime.getRuntime().availableProcessors());
        Random rand = new Random(numTrials);
        SplineDictionary[] randomSplineDicts = new SplineDictionary[numTrials];
        for(int i=0;i<numTrials;i++){
            randomSplineDicts[i] = new SplineDictionary(expressionTable,splineMaker,proteinSet, rand.nextLong());
        }
        for(double threshold: thresholds){
            for(Map.Entry<File,Optimizer> pair: optimizerMap.entrySet()){
                for(SplineDictionary splineDict: randomSplineDicts) {
                    ParallelOptimizer pOpt = new ParallelOptimizer(pair.getValue(), aligner, splineDict,threshold,
                            pair.getKey());
                    executor.execute(pOpt);
                }
            }
        }
        try {
            executor.awaitTermination(7, TimeUnit.DAYS);
        }catch (Exception e){
            e.printStackTrace();
        }


    }
}
