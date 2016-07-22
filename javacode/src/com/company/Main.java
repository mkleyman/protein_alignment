package com.company;

import com.Opitmizers.*;
import com.google.common.collect.TreeBasedTable;
import com.google.common.primitives.Doubles;
import com.parsers.HomologParser;
import com.parsers.ProteinExpressionParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

public class Main {

    public static void main(String[] args){

        try{
            Options options = new Options();
            options.addOption("f", "Function");
            options.addOption("r", "reference file");
            options.addOption("c", "comparison file");
            options.addOption("o", "output file or folder");
            options.addOption("h", "homolog file");
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options,args);
            String ref = cmd.getOptionValue("r");
            String comp = cmd.getOptionValue("comp");
            String outfile = cmd.getOptionValue("o");
            String homolog = cmd.getOptionValue("h");
            if(cmd.getOptionValue("f").equalsIgnoreCase("fdr")){
                FDRMain.run(ref,comp,outfile,homolog);
            }
            else{
                ParamMain.run(ref,comp,outfile,homolog);
            }





        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
