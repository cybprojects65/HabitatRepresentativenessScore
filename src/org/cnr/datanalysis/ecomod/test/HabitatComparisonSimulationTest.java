package org.cnr.datanalysis.ecomod.test;

import java.io.File;
import java.util.Arrays;
import java.util.Random;

import org.cnr.datanalysis.ecomod.habitat.HabitatComparator;
import org.cnr.datanalysis.ecomod.habitat.HabitatRepresentativenessScore;
import org.cnr.datanalysis.ecomod.utils.Operations;

public class HabitatComparisonSimulationTest {

	public static void main(String args[]) throws Exception {

		String basePath = "C:\\Users\\Utente\\Ricerca\\Experiments\\Lago Massaciuccoli\\Step by step workflow and results\\Step 05 bis - Reprojections - final environmental data\\forHRS\\";
		// NOTE: the ASC files should be aligned in space
		// String habitats[] = { "2100", "2016"};
		String habitats[] = { "2016", "2016" }; // RIGHT CONFIGURATION : 2016 AS THE REFERENCE
		int bins = 300;

		boolean doStandardize = true;

		String nonNullexceptions[] = {};
		double HRSMatrix[][] = new double[habitats.length][habitats.length];
		int hi = 0;
		for (String habitat1 : habitats) {
			String habitat1Path = basePath + "\\" + habitat1 + "\\";
			int hj = 0;
			for (String habitat2 : habitats) {

				
				String habitat2Path = basePath + "\\" + habitat2 + "\\";

				File habitat1f = new File(habitat1Path);
				File habitat2f = new File(habitat2Path);

				HabitatComparator ch = new HabitatComparator();
				double[][] fmatrix1 = ch.extractAlignFeatures(habitat1f, nonNullexceptions);
				
				for (int i=0;i<fmatrix1.length;i++) {
					//fmatrix1 [i][0] = fmatrix1 [i][0] + 10*new Random().nextDouble();
					//fmatrix1 [i][1] = fmatrix1 [i][1] + 10*new Random().nextDouble();
					//fmatrix1 [i][2] = fmatrix1 [i][2] + 10*new Random().nextDouble();
					fmatrix1 [i][3] = fmatrix1 [i][3] + 0.01*new Random().nextDouble();
					//fmatrix1 [i][4] = fmatrix1 [i][4] + 10*new Random().nextDouble();
				}
				
				double[][] fmatrix2 = ch.extractAlignFeatures(habitat2f, nonNullexceptions);

				HabitatRepresentativenessScore hrs = new HabitatRepresentativenessScore(bins);

				// THE SECOND HABITAT IS THE REFERENCE, THE FIRST IS THE TEST!
				boolean standardise = true;
				hrs.calcHRS(fmatrix1, fmatrix2, standardise);
				String fnames = Arrays.toString(ch.featureNames.get(0));
				System.out.println("Habitat vector score\n" + fnames + "\n" + Arrays.toString(hrs.HRS_VECTOR));
				System.out.println("HRS of " + habitat1 + " vs " + habitat2 + ": " + hrs.HRS_NFEAT_SUM);
				System.out.println("HRS NFEAT-SUM (HRSVEC): " + hrs.HRS_NFEAT_SUM);
				System.out.println("HRS: " + hrs.HRS_PURE);
				System.out.println("HRS SIMILARITY PERCENTAGE (HRS_PERC): " + hrs.HRS_PERC+"%");

				double hrsNoStand = hrs.HRS_PERC;
				double meanHRS = hrsNoStand;

				if (!hrs.isGoodfit()) {
					System.out.println(
							"#########HRS does not distinguish dimensions using " + bins + " bins - Unfair comparison");
					System.exit(0);
				}

				HRSMatrix[hi][hj] = meanHRS;
				hj++;
			}
			hi++;
			break;
		}

		if (doStandardize) {
			HRSMatrix = new Operations().standardize(HRSMatrix);
			HRSMatrix = Operations.symmetrizeByMean(HRSMatrix);
		}

		StringBuffer sb = new StringBuffer();
		String header = "," + Arrays.toString(habitats).replace("[", "").replace("]", "");
		sb.append(header + "\n");
		for (int i = 0; i < HRSMatrix.length; i++) {
			sb.append(habitats[i] + ",");
			for (int j = 0; j < HRSMatrix[0].length; j++) {
				sb.append(HRSMatrix[i][j]);
				if (j < HRSMatrix[0].length - 1)
					sb.append(",");
				else
					sb.append("\n");
			}
		}

		System.out.println("Final matrix:\n" + sb.toString());

	}

}
