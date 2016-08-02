package com.company;

import java.util.*;

/**
 * Created by mkleyman on 7/26/2016.
 */
public class playground {
    public static void main(String[] args){
        String human = "../processed/human_proteins.csv";
        String mouse = "../processed/mouse_proteins.csv";
        String result_file = "../processed";
        String homologs = "../raw_data/HOM_MouseHumanSequence.rpt";
        Random rand = new Random(500);
        List<Integer> seqList;
        Integer[] sequence = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        for(int i= 0;i<21;i++){
            long lrand = rand.nextLong();
            seqList = Arrays.asList(sequence.clone());
            Collections.shuffle(seqList, new Random(lrand));
            //System.out.println(seqList);
            seqList = Arrays.asList(sequence.clone());
            Collections.shuffle(seqList, new Random(lrand));
            //System.out.println(Arrays.toString(sequence));
            //System.out.println(seqList);
        }
        FDRMain.run(mouse,human,result_file,homologs,'p');


    }
}
