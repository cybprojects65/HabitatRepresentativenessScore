package org.cnr.datanalysis.ecomod.test;

import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;

import org.cnr.datanalysis.ecomod.habitat.HabitatRepresentativenessScore;

public class Main {
	
	
public static void main(String[] args) throws Exception{
		
		File t1 = new File("sampletable1.csv");
		File t2 = new File("sampletable2.csv");

		if (args!=null && args.length>1) {
			
			t1 = new File(args[0]);
			t2 = new File(args[1]);
			
		}
		
		
		double[][] ref = HabitatRepresentativenessScore.loadFeatures(t1);
		double[][] input = HabitatRepresentativenessScore.loadFeatures(t2);
		HabitatRepresentativenessScore hrs = new HabitatRepresentativenessScore();
		
		//NOTE: the two dataset should be normalized
		hrs.calcHRS(input, ref);
		
		System.out.println("Habitat vector score "+Arrays.toString(hrs.HRS_VECTOR));
		System.out.println("HRS "+hrs.HRS);
		
		File output = new File("hrs_summary.txt");
		FileWriter fw = new FileWriter(output);
		fw.write("HRS_VECTOR="+Arrays.toString(hrs.HRS_VECTOR)+"\n");
		fw.write("HRS="+hrs.HRS+"\n");
		
		fw.close();
		
	}


}
