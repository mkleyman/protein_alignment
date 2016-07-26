package com.parsers;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Set;

/**
 * Created by mkleyman on 7/25/2016.
 */
public class GOParser {
    public static Multimap<String,String> parseGO(String filename, Set<String> proteins) throws FileNotFoundException {
        HashMultimap<String,String> goMap = HashMultimap.create();
        Scanner scan = new Scanner(new File(filename));
        int proteinCol = 0;
        int goCol = 6;
        while(scan.hasNextLine()){
            String[] row = scan.nextLine().split("\t");
            if(proteins.contains(row[proteinCol])) {
                String[] goCats = row[goCol].split(";");
                for (String go : goCats) {
                    goMap.put(go, row[proteinCol]);
                }
            }
        }
        scan.close();
        return goMap;
    }
}
