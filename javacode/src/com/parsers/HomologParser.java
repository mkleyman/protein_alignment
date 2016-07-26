package com.parsers;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by mkleyman on 7/15/2016.
 */
public class HomologParser {

    public static Map<String,String> parse(String filename, String reference) throws FileNotFoundException{
        File file = new File(filename);
        Map<String, Integer> refMap = new HashMap<>();
        Map<Integer,String> compMap = new HashMap<>();
        Map<String,String> homologMap = new HashMap<>();
        Scanner scan = new Scanner(file);
        scan.nextLine();
        while(scan.hasNextLine()){

            String[] row = scan.nextLine().split("\t");
            if(row.length==13){
                if(row[1].equals(reference)){
                    refMap.put(row[12],Integer.parseInt(row[0]));
                }
                else{
                    compMap.put(Integer.parseInt(row[0]),row[12]);
                }
            }
        }
        scan.close();
        for(String homolog: refMap.keySet()){
            if(compMap.containsKey(refMap.get(homolog))){
                homologMap.put(homolog, compMap.get(refMap.get(homolog)));
            }
        }
        return homologMap;
    }

    public static Map<String,String> trim(Map<String,String> homologMap, Set ref, Set comp){
        Iterator<Map.Entry<String,String>> it = homologMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry<String,String> pair = it.next();
            if (!(ref.contains(pair.getKey()) && comp.contains(pair.getValue()))){
                it.remove();
            }
        }
        return homologMap;
    }



}
