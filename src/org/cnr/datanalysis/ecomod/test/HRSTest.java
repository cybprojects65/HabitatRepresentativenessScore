package org.cnr.datanalysis.ecomod.test;

import java.io.File;
import java.util.Arrays;

import org.cnr.datanalysis.ecomod.habitat.HabitatRepresentativenessScore;

public class HRSTest {

	
	public static void main(String[] args) throws Exception{
		
		File t1 = new File("sampletable1.csv");
		//File t2 = new File("sampletable2.csv");
		File t2 = new File("sampletable3.csv");
		//File t2 = new File("sampletable1.csv");
		
		double[][] ref = HabitatRepresentativenessScore.loadFeatures(t1);
		double[][] input = HabitatRepresentativenessScore.loadFeatures(t2);
		HabitatRepresentativenessScore hrs = new HabitatRepresentativenessScore();
		
		//NOTE: the two dataset should be normalized
		hrs.calcHRS(input, ref);
		
		System.out.println("Habitat vector score "+Arrays.toString(hrs.HRS_VECTOR));
		System.out.println("HRS "+hrs.HRS);
	}
	
	
	
}
