package org.cnr.datanalysis.ecomod.habitat;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.cnr.datanalysis.ecomod.utils.Operations;
import org.cnr.datanalysis.ecomod.utils.PrincipalComponentAnalysis;

public class HabitatRepresentativenessScore {


	public double[] HRS_VECTOR = null;
	public double HRS_NFEAT_SUM = 0d;
	public double HRS_PURE = 0d;
	public double HRS_PERC = 0d;
	public int maxNofBins = 100;
	public double[][] pcaComponents;
	public PrincipalComponentAnalysis pca;
	
	public boolean goodfit = false;
	
	public boolean isGoodfit() {
		return goodfit;
	}

	public HabitatRepresentativenessScore(int numberOfComparisonBins) {
		this.maxNofBins = numberOfComparisonBins;
	}

	public HabitatRepresentativenessScore() {
		this.maxNofBins = 10;
	}
	
	List<double[]> frequencyDistrib;
	List<double[]> intervals; // uniform subdivision of the numeric ranges
	
	//tables:
	//one feature array for each row
	//one row for each sampling point
	
	public static double[][] loadFeatures(File csvfeatures) throws Exception{
		
		List<String> allines = Files.readAllLines(csvfeatures.toPath());
		double [][] feat = null;
		int rowCounter = 0;
		for (String line:allines) {
		
			String[] elements = line.split(",");
			if (rowCounter==0) {
				feat = new double[allines.size()][elements.length];
			}
			int colCounter = 0;
			for (String e:elements) {
				
				feat[rowCounter][colCounter] = Double.parseDouble(e);
				
				colCounter++;
			}
			rowCounter++;
		}
		
		feat = Operations.normalizeMatrix(feat);
		
		return feat;
		
	}
	
	public void calcHRS(double[][] firstHabitat,double[][] secondHabitat) throws Exception {
		calcHRS(firstHabitat,secondHabitat,true);
	}
	
	//NOTE: standardisation = false should be used only for very small areas to discover micro-variations
	public void calcHRS(double[][] firstHabitat,double[][] secondHabitat, boolean standardise) throws Exception{
		
		Operations operations = new Operations();
		System.out.println("STANDARDIZING REFERENCE HABITAT");
		// standardize the matrix
		double[][] rightAreaPoints = secondHabitat;
		
		if (standardise)
			rightAreaPoints = operations.standardize(secondHabitat);
			
		System.out.println("CALCULATING PCA");
		//calculate PCA
		pca = new PrincipalComponentAnalysis();
		pca.calcPCA(rightAreaPoints);
		// get the pca components for all the vector
		pcaComponents = pca.getProjectionsMatrix(rightAreaPoints);
		
		
		System.out.println("EXTRACTED "+pcaComponents[0].length+" COMPONENTS");
		System.out.println("CALCULATING HISTOGRAMS");
		// calculate the frequency distributions for all the pca: each row will be a frequency distribution for a pca component associated to uniform divisions of the range
		calcFrequenciesDistributionsForComponents(pcaComponents);
		
		System.out.println("STANDARDIZING INPUT HABITAT");
		//double[][] leftAreaPoints = operations.standardize(firstHabitat);
		double[][] leftAreaPoints = firstHabitat;
		if (standardise) {
			System.out.println(operations.means);
			leftAreaPoints = operations.standardize(firstHabitat,operations.means,operations.variances);
			//leftAreaPoints = operations.standardize(firstHabitat);
		}
		System.out.println("PROJECTING INPUT HABITAT INTO THE PCA SPACE");
		//project
		double[][] leftAreaComponents = pca.getProjectionsMatrix(leftAreaPoints);
		
		//calculate frequencies distributions for each component, respect to previous intervals
		int components = leftAreaComponents[0].length;
		
		System.out.println("EVALUATING HISTOGRAM DISCREPANCIES");
		//calculate absolute differences and sum -> obtain a hrs for each PCA component = for each feature
		HRS_VECTOR = new double[components];
		System.out.println("\tNUMBER OF COMPONENTS TO TEST "+components);
		double[][] habitatPcaPointsMatrix = Operations.traspose(leftAreaComponents);
		int ncomparisons = 0;
		for (int i = 0; i < components; i++) {
			double[] habitatPcaPoints = habitatPcaPointsMatrix[i];
			// calculate frequency distributions respect to previous intervals
			double intervalsI [] = intervals.get(i);
			double frequenciesI [] = frequencyDistrib.get(i);
			double[] habitatPcafrequencies = Operations.calcFrequencies(intervalsI , habitatPcaPoints);
			habitatPcafrequencies = Operations.normalizeFrequencies(habitatPcafrequencies, habitatPcaPoints.length);
			System.out.println("\tCOMPONENT "+(i+1)+" - INTERVALS "+intervalsI.length);
			System.out.println("\tHISTO TRAIN "+Arrays.toString(frequenciesI));
			System.out.println("\tHISTO TEST "+Arrays.toString(habitatPcafrequencies));
			double[] absdifference = Operations.vectorialAbsoluteDifference(habitatPcafrequencies, frequenciesI);
			ncomparisons = ncomparisons + absdifference.length; 
			System.out.println("\tDISCREPANCY "+Arrays.toString(absdifference));
			HRS_VECTOR[i] = Operations.sumNormalVector(absdifference,2d);
			System.out.println("\tSUMMED DISCREPANCY "+HRS_VECTOR[i]);
		}
		
		System.out.println("");
		//obtain hrsScore by weighted sum of hrs respect to inverse eigenvalues - too variable, substituted with the sum of the scores
		//HRS = 1-(Operations.sumVector(HRS_VECTOR)/(double)ncomparisons);
		HRS_PURE = Operations.sumVector(HRS_VECTOR);
		HRS_NFEAT_SUM = (double) (HRS_VECTOR.length)-(double) Operations.sumVector(HRS_VECTOR);
		HRS_PERC = 100d-Operations.roundDecimal(HRS_PURE*100d / (double) HRS_VECTOR.length,2);
		
		//alternative HRSs
			//HRS= Operations.scalarProduct(HRS_VECTOR, pca.getInverseNormalizedEigenvalues());
			//HRS = 1-(Operations.sumVector(HRS_VECTOR)/(double)HRS_VECTOR.length);
		
		
		System.out.println("TOTAL HRS SCORE: "+HRS_PURE);
		System.out.println("TOTAL HRS PURE SCORE DIFFERENCE: "+HRS_NFEAT_SUM);
		
	}
	
	
	// calculate a frequency distribution for each component
	public void calcFrequenciesDistributionsForComponents(double[][] pcaComponents) {
		frequencyDistrib = null;
		if (pcaComponents.length > 0) {
			int sizeDistrib = pcaComponents[0].length;
			frequencyDistrib = new ArrayList<double[]>();
			intervals = new ArrayList<double[]>();
			double[][] pcaColumns = Operations.traspose(pcaComponents);
			goodfit = true;
			for (int i = 0; i < sizeDistrib; i++) {
				double[] pcaPoints = pcaColumns[i];
				double[] interval = Operations.uniformDivide(Operations.getMax(pcaPoints), Operations.getMin(pcaPoints), pcaPoints,maxNofBins);
				if (interval.length==maxNofBins) {
					goodfit = false;
					System.out.println("The maximum number of bins is also the optimal - there is margin for improvement");
				}
				double[] frequencies = Operations.calcFrequencies(interval, pcaPoints);
				frequencies = Operations.normalizeFrequencies(frequencies, pcaPoints.length);
				intervals.add(interval);
				frequencyDistrib.add(frequencies);
			}
			
		}
	}
}
