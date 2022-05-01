package org.cnr.datanalysis.ecomod.featureselection;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import org.cnr.datanalysis.ecomod.utils.Operations;
import org.cnr.datanalysis.ecomod.utils.PrincipalComponentAnalysis;

import it.cnr.raster.asc.AscRaster;
import it.cnr.raster.asc.AscRasterReader;

public class PCAInspectorMultiVariate {

	public List<String[]> featureNames = new ArrayList<>();
	
	public double[][] extractAlignFeatures(File [] ascfeaturefiles, String featuresWithNullEqualToZero []) throws Exception{
		
		featureNames = new ArrayList<>();
		
		double resolution = -1;
		int nrows = 0;
		int ncols = 0;
		double nodata = -9999;
		
		//1 column per feature = 1 feature array per row
		LinkedHashMap<Integer, Double[]> featurearray = new LinkedHashMap<>();
		int featureCounter = 0;
		int nFeatures = ascfeaturefiles.length;
		String [] fnames = new String[nFeatures];
		
		for (File f:ascfeaturefiles) {
			//for each input file = spatial feature
			System.out.println("Reading file "+f.getName());
			fnames[featureCounter] = f.getName();
			//read the ASC file
			AscRaster asc =  new AscRasterReader().readRaster(f.getAbsolutePath());
			nodata = Double.parseDouble(asc.getNDATA());
			//if it is the first file, extract resolution information->note that all files should be at the same spatial resolution
			if (resolution == -1) {
				resolution = asc.getCellsize();
				System.out.println("Spatial data resolution of the comparison "+resolution);
				 nrows = asc.getRows();
				 ncols = asc.getCols();
			}
			
			//get the ASC data
			double [][] data = asc.getData();
			
			int featurerowcounter = 0;
			System.out.println("Populating feature array with "+(nrows*ncols)+" data and "+nFeatures+" features");
			//for the initial grid of elements populate the feature matrix by adding elements to the rightmost column 
			for (int i=0;i<nrows;i++) {
				for (int j=0;j<ncols;j++) {
					//get the previously populated row
					Double [] dd = featurearray.get(featurerowcounter);
					//build up the row if there was any
					if (dd==null) {
						dd = new Double[nFeatures];
					}
					//take the current ASC element
					double dat = data [i][j];
					//if the ASC is among those with null=0 add 0
					for (String nonNull:featuresWithNullEqualToZero) {
						if (f.getName().startsWith(nonNull) && (dat == nodata))
							dat = 0;
					}
					//add the element at the end of the array
					dd[featureCounter] = dat;
					//update the row
					featurearray.put(featurerowcounter,dd);
					//pass to the next row
					featurerowcounter++;
				}
			}
			//pass to the next feature
			featureCounter ++;
		}
		
		System.out.println("Removing NODATA");
		
		//remove lines with nodata
		int idxmax = featurearray.size();
		int removedLines = 0;
		for (int idx=0;idx<idxmax;idx++) {
			Double [] features = featurearray.get(idx);
			for (Double f:features) {
				 if (f == nodata) {
					featurearray.remove(idx);
					removedLines++;
					break;
				}
			}
			
		}
		
		//build a matrix with the all-populated rows
		System.out.println("Removed "+removedLines+" over "+idxmax+"; Remaining "+(idxmax-removedLines)+" ("+Math.round((idxmax-removedLines)*100/idxmax)+"%)");
		ArrayList<Double[]> aslist = new ArrayList<Double[]>(featurearray.values());
		
		double [][] outputMatrix = new double [aslist.size()][nFeatures];
		for (int i=0;i<outputMatrix.length;i++) {
			for (int j=0;j<nFeatures;j++) {
				outputMatrix [i][j] = aslist.get(i)[j];
			}
			
		}
		
		//update the feature list names
		featureNames.add(fnames);
		return outputMatrix;
	}
	
	
	double sigma1;
	double sigma2;
	double MSE1;
	double MSE2;
	double percweight1;
	double percweight2;
	double weight1;
	double weight2;
	public List<Integer> valuableFeatures = new ArrayList<>();
	public List<Double> valuablePercentages = new ArrayList<>();
	public List<String> valuableNames = new ArrayList<>();
	
	//feature matrix->matrix with one feature for each column
	public void PCACompare (double [][] featureMatrix, double cumulativePercentContributionEigenvalueThreshold, double individualPercentContributionVariableThreshold) throws Exception{
		
		//run PCA
		Operations operations = new Operations();
		System.out.println("Standardizing Matrix");
		// standardize the matrix
		double[][] stdmatrix = operations.standardize(featureMatrix);
		System.out.println("CALCULATING PCA");
		//calculate PCA
		PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
		pca.calcPCA(stdmatrix);
		//estimate the load from first eigen
		double e [][] =  new double[featureMatrix[0].length][featureMatrix[0].length];
		double load [][] =  new double[featureMatrix[0].length][featureMatrix[0].length];
		double lambda [] =  new double[featureMatrix[0].length];
		double weights [] =  new double[featureMatrix[0].length];
		double percweights [] =  new double[featureMatrix[0].length];
		
		System.out.println("Eigenval X Eigenvector");
		double totalweighteigen = 0;
		for (int i=0;i<e.length;i++) {
			e[i] = pca.getEigenvector(i);
			lambda [i] = pca.getEigenvalue(i);
			totalweighteigen+=lambda [i];
			System.out.println(Operations.roundDecimal(lambda [i], 2)+" X "+Arrays.toString(e[i]));
		}
		int eigenMaxIndex = e.length-1;
		double cumulativeLambda = 0;
		for (int i=0;i<e.length;i++) {
			double lambdaperc = lambda [i]*100/totalweighteigen;
			cumulativeLambda += lambdaperc;
			
			System.out.print(Operations.roundDecimal(lambda [i], 2)+" = "+Operations.roundDecimal(lambdaperc, 2)+"% Cumulative "+Operations.roundDecimal(cumulativeLambda,2)+"%");
			//if (lambdaperc<percentContributionEigenvalueThreshold) {
			if (cumulativeLambda>cumulativePercentContributionEigenvalueThreshold) {
				eigenMaxIndex = i-1;
				System.out.println(" -> Excluded");
				break;
			}else
				System.out.println("");
		}
		
		System.out.println("Loadings");
		for (int i=0;i<e.length;i++) {
			for (int j=0;j<e[0].length;j++) {
				load [i][j]=e[i][j]*Math.sqrt(lambda[i]);
			}
			System.out.println(Arrays.toString(load[i]));
		}
		
		System.out.println("Average loading per variable -> cutting eigen threshold ");
		double totalweight = 0;
		for (int j=0;j<e[0].length;j++) {
			double weigthfeature = 0;
			for (int i=0;i<=eigenMaxIndex;i++) {
				weigthfeature = weigthfeature+ load [i][j];
			}
			weigthfeature = weigthfeature /(double) (eigenMaxIndex+1);
			weights[j] = weigthfeature;
			System.out.println("F"+(j+1)+" avg weight = "+Operations.roundDecimal(weigthfeature,2));
			if (weigthfeature>0)
				totalweight = totalweight+weigthfeature;
		}
		
		System.out.println("Keeping features with contribution over "+individualPercentContributionVariableThreshold);
				
		for (int i=0;i<weights.length;i++) {
			if (weights[i]>0) {
				percweights [i] = weights[i]*100d/totalweight;
			}else
				percweights [i] = 0;
			if (percweights [i]>individualPercentContributionVariableThreshold) {
				int position = 0;
				
				for (Double p:valuablePercentages) {
					if (p<percweights [i]) {
						break;
					}
					position ++;
				}
				valuableFeatures.add(position,i);
				valuablePercentages.add(position,percweights [i]);
				valuableNames.add(position,featureNames.get(0)[i]);
				
				System.out.println("w"+(i)+": "+Operations.roundDecimal(weights[i],2) + " ("+Operations.roundDecimal(percweights [i],2)+"%)"+" - "+featureNames.get(0)[i]+"->VALUABLE");
			}else
				System.out.println("w"+(i)+": "+Operations.roundDecimal(weights[i],2) + " ("+Operations.roundDecimal(percweights [i],2)+"%)"+" - "+featureNames.get(0)[i]);
		}
		
		System.out.println("Ranking:");
		System.out.println("Feature Names "+valuableNames.toString());
		System.out.println("Feature Indexes "+valuableFeatures.toString());
		System.out.println("Feature Percentages "+valuablePercentages.toString());
	}
	
	public static void main(String[] args) throws Exception{
		String nonNullexceptions[] = { "Sea_Ice_Concentration" };
		String year = "2017";
		//String area = "Adriatic_Sea";
		//String area = "Aegean_Sea";
		//String area = "Baltic_Sea"; //90%
		//String area = "Bay_of_Biscay";
		//String area = "Black_Sea";
		//String area = "Levantine_Sea";
		//String area = "North_Sea";
		String area = "Western_Mediterranean_Sea";
		double cumulativePercentContributionEigenvalueThreshold = 95;
		double individualPercentContributionVariableThreshold = 1;
		
		//String basepath = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS\\Adriatic_Sea\\"+year+"\\";
		String basepath = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS\\"+area+"\\"+year+"\\";
		int casef = 10;
		File [] featureFiles = new File[2]; 
		//take two features
		if (casef ==1) {
			featureFiles[0] = new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = new File(basepath+"Sea-surface_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		
		}
		if (casef ==2) {
			featureFiles[1] = new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[0] = new File(basepath+"Sea-surface_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}
		if (casef ==3) {
			
			featureFiles[0] = 	new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}

		if (casef == 4) {
			featureFiles[0] = 	new File(basepath+"Net_Primary_Production_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = 	new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}
		if (casef == 5) {
			featureFiles[0] = 	new File(basepath+"Sea_Ice_Concentration_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = 	new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}

		if (casef == 6) {
			featureFiles[0] = 	new File(basepath+"Sea-bottom_dissolved_oxygen_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = 	new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}
		
		if (casef == 7) {
			featureFiles[0] = 	new File(basepath+"Sea-surface_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = 	new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}
		
		if (casef == 8) {
			featureFiles[0] = 	new File(basepath+"Sea-bottom_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = 	new File(basepath+"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}
		
		if (casef == 9) {
			featureFiles[0] = 	new File(basepath+"Sea-bottom_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
			featureFiles[1] = 	new File(basepath+"Net_Primary_Production_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc");
		}
		if (casef == 10) {
			featureFiles = new File[7];
			featureFiles[0] = new File(basepath+"Sea-bottom_salinity_res_01_annual_years_"+year+"_Clim_scen_today_regional_"+area+".asc");
			featureFiles[1] = new File(basepath+"Sea-surface_temperature_res_01_annual_years_"+year+"_Clim_scen_today_regional_"+area+".asc");
			featureFiles[2] = 	new File(basepath+"Net_Primary_Production_res_01_annual_years_"+year+"_Clim_scen_today_regional_"+area+".asc");
			featureFiles[3] = 	new File(basepath+"Sea_Ice_Concentration_res_01_annual_years_"+year+"_Clim_scen_today_regional_"+area+".asc");
			featureFiles[4] = 	new File(basepath+"Sea-bottom_dissolved_oxygen_res_01_annual_years_"+year+"_Clim_scen_today_regional_"+area+".asc");
			featureFiles[5] = 	new File(basepath+"Sea-surface_salinity_res_01_annual_years_"+year+"_Clim_scen_today_regional_"+area+".asc");
			featureFiles[6] = 	new File(basepath+"Sea-bottom_temperature_res_01_annual_years_"+year+"_Clim_scen_today_regional_"+area+".asc");
		}
		
		//TODO: fix values to the same spatial grid
		//analyse the variability of the other seas
		//analyse the variability per year
		PCAInspectorMultiVariate inspector = new PCAInspectorMultiVariate();
		double [][] featureMatrix = inspector.extractAlignFeatures(featureFiles, nonNullexceptions);
		inspector.PCACompare(featureMatrix,cumulativePercentContributionEigenvalueThreshold,individualPercentContributionVariableThreshold);
		
		int counter = 0;
		
		for (String variable:inspector.valuableNames) {
			
			String v = variable.substring(0,variable.indexOf("_res"));
			v = v.replace("_", " ");
			double val = inspector.valuablePercentages.get(counter);
			
			System.out.println(v+" ("+Operations.roundDecimal(val, 1)+"%)");
			counter++;
		}
	}	
	
}
