package com.Opitmizers;

import com.company.Aligner;
import com.company.SplineDictionary;
import org.apache.commons.math3.analysis.UnivariateFunction;

/**
 * Created by mkleyman on 7/19/2016.
 */
public abstract class Optimizer {


    abstract public void setMode(char mode);

    abstract public double[] optimizePearson(Aligner aligner, SplineDictionary splineDict,
                                           double threshold);
    abstract public String getName();

    double optVal(Aligner aligner,SplineDictionary splineDict, UnivariateFunction poly,
                  double threshold, char mode){
        if(mode=='a' ) {
            return aligner.align_polynomial_pearson_all(splineDict,poly, threshold);
        }
        else if(mode=='d'){
            return aligner.align_polynomial_difference(splineDict,poly,threshold);
        }
        else if(mode=='i'){
            return aligner.align_polynomial_info(splineDict,poly);
        }
        else if(mode=='s'){
            return aligner.align_polynomial_sp(splineDict,poly,threshold);
        }

        else{
            return aligner.align_polynomial_pearson(splineDict, poly, threshold);
        }
    }
}
