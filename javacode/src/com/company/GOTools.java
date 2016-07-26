package com.company;

import com.google.common.collect.Multimap;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by mkleyman on 7/25/2016.
 */
public class GOTools {

    public static void createGOTable(Aligner aligner, SplineDictionary spDict,UnivariateFunction funct,
                                     double threshold, Multimap<String,String> GOmap, String filename)
            throws IOException {
        FileWriter fwriter = new FileWriter(new File(filename));
        Set<String> proteins = aligner.align_pearson_set(spDict,funct,threshold);
        for(String GoCat: GOmap.keys()){
            HashSet<String> goProteins = new HashSet<>();
            goProteins.addAll(GOmap.get(GoCat));
            goProteins.retainAll(proteins);
            int inter = goProteins.size();
            int total = GOmap.get(GoCat).size();
            String line = GoCat+","+inter+","+total+"\n";
            fwriter.write(line);
        }
        fwriter.close();
    }
}
