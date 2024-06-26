package org.cnr.datanalysis.ecomod.test;

import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;

import org.cnr.datanalysis.ecomod.habitat.HabitatRepresentativenessScore;

public class Main {
	
	
public static void main(String[] args) throws Exception{
		
		File t1 = new File("sampletable1.csv");
		File t2 = new File("sampletable2.csv");
		int numberOfBins = 10;
		
		System.out.println("Running HRS");
		
		if (args!=null && args.length>1) {
			
			t1 = new File(args[0]);
			t2 = new File(args[1]);
			numberOfBins = Integer.parseInt(args[2]);
		}
		
		
		double[][] ref = HabitatRepresentativenessScore.loadFeatures(t2);
		double[][] input = HabitatRepresentativenessScore.loadFeatures(t1);
		HabitatRepresentativenessScore hrs = new HabitatRepresentativenessScore(numberOfBins);
		
		//NOTE: the two dataset should be normalized
		hrs.calcHRS(input, ref);
		
		System.out.println("Habitat vector score: "+Arrays.toString(hrs.HRS_VECTOR));
		System.out.println("HRS NFEAT-SUM (HRSVEC): "+hrs.HRS_NFEAT_SUM);
		System.out.println("HRS: "+hrs.HRS_PURE);
		System.out.println("HRS SIMILARITY PERCENTAGE (HRS_PERC): "+hrs.HRS_PERC);
		
		File output = new File("hrs_summary.txt");
		FileWriter fw = new FileWriter(output);
		fw.write("HRS_VECTOR="+Arrays.toString(hrs.HRS_VECTOR)+"\n");
		fw.write("HRS="+hrs.HRS_PURE+"\n");
		fw.write("HRS_NFEATURES_MINUS_HRS="+hrs.HRS_NFEAT_SUM+"\n");
		fw.write("HRS_SIMILARITY_PERC="+hrs.HRS_PERC+"%\n");
		
		fw.close();
		
	}


}
