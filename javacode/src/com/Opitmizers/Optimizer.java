package com.Opitmizers;

import com.company.Aligner;
import com.company.SplineDictionary;

/**
 * Created by mkleyman on 7/19/2016.
 */
public interface Optimizer {

    double[] optimizePearson(Aligner aligner, SplineDictionary splineDict,
                                           double threshold);
    String getName();
}
