package com.company;

import com.google.common.collect.Multimap;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.Set;

/**
 * Created by mkleyman on 7/25/2016.
 */
public class GOTools {

    public static void createGOTable(Aligner aligner, SplineDictionary spDict,UnivariateFunction funct,
                                     double threshold, Multimap<String,String> GOmap){
        Set<String> proteins = aligner.align_pearson_set(spDict,funct,threshold);
        for(String GoCat: GOmap.keys()){
            int inter;
        }
    }
}
