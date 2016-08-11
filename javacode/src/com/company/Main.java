package com.company;

import com.parsers.ProteinExpressionParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;


import java.io.File;
import java.util.*;

public class Main {

    public static void main(String[] args){
        System.out.println(Arrays.toString(args));
        try{
            Options options = new Options();
            options.addOption("f",true, "function");
            options.addOption("r",true, "reference file");
            options.addOption("c",true, "comparison file");
            options.addOption("o",true, "output file or folder");
            options.addOption("h",true, "homolog file");
            //options.addOption("m",true, "mode");
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options,args);
            String ref = cmd.getOptionValue("r");
            System.out.println(ref);
            String comp = cmd.getOptionValue("c");
            System.out.println(comp);
            String outfile = cmd.getOptionValue("o");
            System.out.println(outfile);
            String homolog = cmd.getOptionValue("h");
            System.out.println(homolog);
            char mode = 'p';
            System.out.println(mode);
            if(cmd.getOptionValue("f").equalsIgnoreCase("fdr")) {
                FDRMain.run(ref, comp, outfile, homolog, mode);
            }
            else if(cmd.getOptionValue("f").equalsIgnoreCase("fdrq")){
                    FDRQuadMain.run(ref,comp,outfile,homolog,mode);
                }
            else{
                ParamMain.run(ref, comp, outfile, homolog, mode);
            }





        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
