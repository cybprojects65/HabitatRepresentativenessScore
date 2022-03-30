package org.cnr.datanalysis.ecomod.test;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.cnr.datanalysis.ecomod.habitat.HabitatComparator;
import org.cnr.datanalysis.ecomod.habitat.HabitatRepresentativenessScore;
import org.cnr.datanalysis.ecomod.utils.Operations;

public class HabitatComparisonTwoEuropeanSeas {

	public static void main(String args[]) throws Exception {

		String basePath = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS input data";
		//NOTE: the ASC files should be aligned
		String habitats[] = { "Adriatic_Sea","Aegean_Sea"};
		int year = 2018;
		int bins = 200;
		
		boolean doPermutations = false;
		
		String nonNullexceptions[] = { "Sea_Ice_Concentration" };
		double HRSMatrix[][] = new double[habitats.length][habitats.length];
		int hi = 0;
		List<double []>  hrsvectors = new ArrayList<>();
		
		for (String habitat1 : habitats) {
			String habitat1Path = basePath + "\\" + habitat1 + "\\" + year;
			int hj = 0;
			for (String habitat2 : habitats) {
				
				if (habitat1.equals(habitat2)) {
					HRSMatrix[hi][hj] = 1;
					hrsvectors.add(new double[0]);
					hj++;
					continue;
				}
				
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
				
				hrsvectors.add(hrs.HRS_VECTOR);
				double[] firsteigen = Operations.normalize(hrs.pca.getEigenvector(0));
				
				System.out.println("E0:"+Arrays.toString(firsteigen));
				//System.exit(0);
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
				
				
				double covariance [][] = new double[fmatrix1[0].length][fmatrix1[0].length];
				double[][] fmatrix1T = Operations.traspose(fmatrix1);
				double[][] fmatrix2T = Operations.traspose(fmatrix2);
				
				for (int f1=0;f1<fmatrix1[0].length;f1++) {
					double [] v1 = fmatrix1T[f1];
					
					for (int f2=0;f2<fmatrix1[0].length;f2++) {
						
						double [] v2= fmatrix2T[f2];
						if (f1==f2) {
							double c = Operations.covariance(v1, v2);
							covariance[f1][f2] = c;
						}
					}
					
				}
				
				System.out.println("COVARIANCE\n");
				for (int ff=0;ff<covariance.length;ff++) {
					System.out.println(Arrays.toString(covariance[ff]));
				}
				HRSMatrix[hi][hj] = meanHRS;
				hj++;
			}
			hi++;
		}
		
		
		System.out.println("############################DONE############################");
		
		//HRSMatrix = new Operations().standardize(HRSMatrix);
		//HRSMatrix = Operations.symmetrizeByMean(HRSMatrix);
		//HRSMatrix = Operations.symmetrizeByMax(HRSMatrix);
		
		StringBuffer sb = new StringBuffer();
		StringBuffer sb2 = new StringBuffer();
		
		String header = ","+Arrays.toString(habitats).replace("[", "").replace("]", "");
		sb.append(header+"\n");
		sb2.append(header+"\n");
		int arrayIdx = 0;
		for (int i=0;i<HRSMatrix.length;i++) {
			sb.append(habitats[i]+",");
			sb2.append(habitats[i]+",");
			for (int j=0;j<HRSMatrix[0].length;j++) {
				sb.append(HRSMatrix[i][j]);
				sb2.append("\""+Arrays.toString(hrsvectors.get(arrayIdx)).replace("[", "").replace("]", "")+"\"");
				
				if (j<HRSMatrix[0].length-1) {
					sb.append(",");
					sb2.append(",");
				}
				else {
					sb.append("\n");
					sb2.append("\n");
				}
				arrayIdx++;
			}
		}
		
		
		
		String outFile = basePath+"\\"+"HRSComparisonMatrix_"+habitats[0]+"_vs_"+habitats[1]+"_Perm_"+doPermutations+".csv";
		File outFileF = new File(outFile);
		FileWriter fw = new FileWriter(outFileF);
		fw.write(sb.toString());
		fw.close();
		
		outFile = basePath+"\\"+"HRSComparisonMatrix_"+habitats[0]+"_vs_"+habitats[1]+"_Perm_"+doPermutations+"_VECTORS.csv";
		outFileF = new File(outFile);
		fw = new FileWriter(outFileF);
		fw.write(sb2.toString());
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
