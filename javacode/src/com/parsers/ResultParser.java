package com.parsers;

import com.TransformationFunctions.ExpTransform;
import com.TransformationFunctions.LogarithmTransform;
import com.TransformationFunctions.LogitTransform;
import com.TransformationFunctions.SqrtTransform;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;

/**
 * Created by mkleyman on 7/25/2016.
 */
public class ResultParser {

    public static void createTable(String resultFile, String tableFile,
                                   double start, double stop, double increment) throws IOException {
        Scanner scan = new Scanner(new File(resultFile));
        File table =  new File(tableFile);
        table.createNewFile();
        FileWriter fwriter = new FileWriter(table);
        while(scan.hasNext()){
            String[] row = scan.nextLine().split(",");
            UnivariateFunction funct = rowToFunction(row);
            for(double i=start;i<stop;i+=increment){
                System.out.println(row[0]);
                System.out.println(i);
                if(funct instanceof PolynomialFunction){
                    PolynomialFunction x = (PolynomialFunction)funct;
                    System.out.println(Arrays.toString(x.getCoefficients()));
                }
                System.out.println(funct.value(i)/365.0);
                String[] tableRow = {row[0],row[1],row[2],Double.toString(i), Double.toString(funct.value(i))};
                fwriter.write(String.join(",", tableRow)+"\n");
            }
        }
        scan.close();
        fwriter.close();

    }

    public static UnivariateFunction rowToFunction(String[] row){
        if(row[0].equals("Linear")){
            double[] poly_args = {Double.parseDouble(row[3]),Double.parseDouble(row[4])};
            return  new PolynomialFunction(poly_args);
        }
        else if(row[0].equals("Quadratic")){
            double[] poly_args = {Double.parseDouble(row[3]),Double.parseDouble(row[4]),
                    Double.parseDouble(row[5])};
            return  new PolynomialFunction(poly_args);
        }
        else if(row[0].equals("Cubic")){
            double[] poly_args = {Double.parseDouble(row[3]),Double.parseDouble(row[4]),
                    Double.parseDouble(row[5]),Double.parseDouble(row[6])};
            return  new PolynomialFunction(poly_args);
        }
        else if(row[0].equals("Exp")){
            return  new ExpTransform(Double.parseDouble(row[3]), Double.parseDouble(row[4]));
        }
        else if(row[0].equals("Log")){
            return  new LogarithmTransform(Double.parseDouble(row[3]), Double.parseDouble(row[4]));
        }
        else if(row[0].equals("Logit")){
            return  new LogitTransform(Double.parseDouble(row[3]), Double.parseDouble(row[4]),
                    Double.parseDouble(row[5]) );
        }
        else if(row[0].equals("Sqrt")){
            return  new SqrtTransform(Double.parseDouble(row[3]), Double.parseDouble(row[4]));
        }

        return null;
    }

    public static void main(String[] args){
        String resultFile = "../processed/search_results.txt";
        String tableFile = "../processed/search_table.csv";
        try {
            createTable(resultFile,tableFile,16,48,0.1);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }



}
