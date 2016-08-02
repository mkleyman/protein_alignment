package com.company;

import com.parsers.HomologParser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by mkleyman on 7/26/2016.
 */
public class HomologSignificance {
    public static void main(String args[]){
        String result_file = "../processed/significance_results.csv";
        String homologs = "../raw_data/HOM_MouseHumanSequence.rpt";
        String human = "../raw_data/human_all.txt";
        String mouse = "../raw_data/mouse_all.txt";;

        try {
            FileWriter fwriter = new FileWriter(new File(result_file));
            Map<String,String> homologMap = HomologParser.parse(homologs, "mouse, laboratory") ;
            List<String> humanProteins = readSwissProt(human);
            System.out.println("done reading human");
            List<String> mouseProteins = readSwissProt(mouse);
            System.out.println("done reading mouse");
            for(int i=0;i<10000;i++){
                int total = 0;
                Collections.shuffle(humanProteins);
                Collections.shuffle(mouseProteins);
                List<String> mouseSublist= mouseProteins.subList(0,1020);
                Set<String> humanSubset = new HashSet(humanProteins.subList(0,848));
                for(String mouseProtein: mouseSublist) {
                    if (homologMap.containsKey(mouseProtein) &&
                            humanSubset.contains(homologMap.get(mouseProtein))) {

                        total++;
                    }
                }
                String line = "random,"+total+"\n";
                fwriter.write(line);

            }
            fwriter.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static List<String> readSwissProt(String filename) throws FileNotFoundException {
        LinkedList<String> proteins = new LinkedList<>();
        Scanner scan = new Scanner(new File(filename));
        scan.nextLine();
        while(scan.hasNextLine()){
            proteins.add(scan.nextLine());
        }
        scan.close();
        return proteins;
    }


}
