package org.cnr.datanalysis.ecomod.test;

import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;

import org.cnr.datanalysis.ecomod.habitat.HabitatComparator;
import org.cnr.datanalysis.ecomod.habitat.HabitatRepresentativenessScore;
import org.cnr.datanalysis.ecomod.utils.Operations;

public class HabitatComparisonAllEuropeanSeas {

	public static void main(String args[]) throws Exception {

		String basePath = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS input data";
		//NOTE: the ASC files should be aligned
		String habitats[] = { "Adriatic_Sea", "Aegean_Sea", "Baltic_Sea", "Bay_of_Biscay", "Black_Sea", "Levantine_Sea",
				"North_Sea", "Western_Mediterranean_Sea" };
		int year = 2020;
		int bins = 200;
		
		boolean doPermutations = false;
		
		String nonNullexceptions[] = { "Sea_Ice_Concentration" };
		double HRSMatrix[][] = new double[habitats.length][habitats.length];
		int hi = 0;
		for (String habitat1 : habitats) {
			String habitat1Path = basePath + "\\" + habitat1 + "\\" + year;
			int hj = 0;
			for (String habitat2 : habitats) {
				/*
				if (habitat1.equals(habitat2)) {
					HRSMatrix[hi][hj] = 1;
					hj++;
					continue;
				}
				*/
				String habitat2Path = basePath + "\\" + habitat2 + "\\" + year;

				File habitat1f = new File(habitat1Path);
				File habitat2f = new File(habitat2Path);

				HabitatComparator ch = new HabitatComparator();
				double[][] fmatrix1 = ch.extractAlignFeatures(habitat1f, nonNullexceptions);
				double[][] fmatrix2 = ch.extractAlignFeatures(habitat2f, nonNullexceptions);

				fmatrix1 = Operations.normalizeMatrix(fmatrix1);
				fmatrix2 = Operations.normalizeMatrix(fmatrix2);

				HabitatRepresentativenessScore hrs = new HabitatRepresentativenessScore(bins);

				// NOTE: the two dataset should be normalized
				hrs.calcHRS(fmatrix1, fmatrix2);
				String fnames = Arrays.toString(cleanupStrings(ch.featureNames.get(0)));
				System.out.println("Habitat vector score\n" + fnames + "\n" + Arrays.toString(hrs.HRS_VECTOR));
				System.out.println("HRS of "+habitat1+" vs "+habitat2 +": "+ hrs.HRS);
				double meanHRS = 0;
				meanHRS = hrs.HRS;
				if (!hrs.isGoodfit()) {
					System.out.println("#########HRS does not distinguish dimensions using "+bins+" bins - Unfair comparison");
					System.exit(0);
				}
				if (doPermutations) {
					
					int nComparisons = 0;
					for (int p1 = 0; p1 < fmatrix1[0].length; p1++) {

						for (int p2 = p1 + 1; p2 < fmatrix1[0].length; p2++) {
							double f1[][] = Operations.permuteColumns(fmatrix1, p1, p2);
							double f2[][] = Operations.permuteColumns(fmatrix2, p1, p2);
							hrs = new HabitatRepresentativenessScore(bins);
							hrs.calcHRS(f1, f2);
							meanHRS = meanHRS + hrs.HRS;
							nComparisons++;
						}

					}
					meanHRS = meanHRS / (double) nComparisons;
					System.out.println("HRS after permutation " + meanHRS);
				}
				
				HRSMatrix[hi][hj] = meanHRS;
				hj++;
			}
			hi++;
		}
		
		
		System.out.println("############################DONE############################");
		
		HRSMatrix = new Operations().standardize(HRSMatrix);
		//HRSMatrix = Operations.symmetrizeByMean(HRSMatrix);
		HRSMatrix = Operations.symmetrizeByMax(HRSMatrix);
		
		StringBuffer sb = new StringBuffer();
		String header = ","+Arrays.toString(habitats).replace("[", "").replace("]", "");
		sb.append(header+"\n");
		for (int i=0;i<HRSMatrix.length;i++) {
			sb.append(habitats[i]+",");
			for (int j=0;j<HRSMatrix[0].length;j++) {
				sb.append(HRSMatrix[i][j]);
				if (j<HRSMatrix[0].length-1)
					sb.append(",");
				else
					sb.append("\n");
			}
		}
		
		
		
		String outFile = basePath+"\\"+"HRSComparisonMatrix"+"_Perm_"+doPermutations+".csv";
		File outFileF = new File(outFile);
		FileWriter fw = new FileWriter(outFileF);
		fw.write(sb.toString());
		fw.close();
	}

	public static String[] cleanupStrings(String[] ss) {
		int i = 0;
		for (String s : ss) {

			s = s.substring(0, s.indexOf("_res"));
			s = s.replace("_", " ");
			ss[i] = s;
			i++;
		}

		return ss;
	}

}
