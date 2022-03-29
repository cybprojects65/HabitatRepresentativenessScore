package org.cnr.datanalysis.ecomod.test;

import java.io.File;
import java.util.Arrays;

import org.cnr.datanalysis.ecomod.habitat.HabitatComparator;
import org.cnr.datanalysis.ecomod.habitat.HabitatRepresentativenessScore;
import org.cnr.datanalysis.ecomod.utils.Operations;

public class HabitatComparisonTwoEuropeanSeas {

public static void main (String args[]) throws Exception{
		
		String habitat1Path = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS input data\\Western_Mediterranean_Sea\\2020";
		String habitat2Path = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS input data\\Levantine_Sea\\2020";
		
		String nonNullexceptions [] = {"Sea_Ice_Concentration"};
		
		File habitat1 = new File(habitat1Path);
		File habitat2 = new File(habitat2Path);
		HabitatComparator ch = new HabitatComparator();
		double [][] fmatrix1 = ch.extractAlignFeatures(habitat1, nonNullexceptions);
		double [][] fmatrix2 = ch.extractAlignFeatures(habitat2, nonNullexceptions);
		
		fmatrix1 = Operations.normalizeMatrix(fmatrix1);
		fmatrix2 = Operations.normalizeMatrix(fmatrix2);
		
		HabitatRepresentativenessScore hrs = new HabitatRepresentativenessScore();
		
		//NOTE: the two dataset should be normalized
		hrs.calcHRS(fmatrix1, fmatrix2);
		String fnames = Arrays.toString(cleanupStrings(ch.featureNames.get(0)));
		System.out.println("Habitat vector score\n"+fnames+"\n"+Arrays.toString(hrs.HRS_VECTOR));
		System.out.println("HRS "+hrs.HRS);
		
		double meanHRS = 0;
		int nComparisons = 0;
		for (int p1=0;p1<fmatrix1[0].length;p1++) {
		
			for (int p2=p1+1;p2<fmatrix1[0].length;p2++) {
				double f1[][] = Operations.permuteColumns(fmatrix1, p1, p2);
				double f2[][] = Operations.permuteColumns(fmatrix2, p1, p2);
				hrs = new HabitatRepresentativenessScore();
				hrs.calcHRS(f1, f2);
				meanHRS = meanHRS+hrs.HRS;
				nComparisons++;
			}
			
		}
		meanHRS = meanHRS/(double)nComparisons;
		
		System.out.println("HRS after permutation "+meanHRS);
		
	}

	public static String [] cleanupStrings(String[] ss) {
		int i = 0;
		for (String s:ss) {
			
			s = s.substring(0,s.indexOf("_res"));
			s = s.replace("_", " ");
			ss[i] = s;
			i++;
		}
		
		return ss;
	}

}
