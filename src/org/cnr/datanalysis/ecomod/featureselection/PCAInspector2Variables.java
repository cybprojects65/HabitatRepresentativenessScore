package org.cnr.datanalysis.ecomod.featureselection;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import org.cnr.datanalysis.ecomod.utils.Operations;
import org.cnr.datanalysis.ecomod.utils.PrincipalComponentAnalysis;

import it.cnr.raster.asc.filemanagement.AscRaster;
import it.cnr.raster.asc.filemanagement.AscRasterReader;

public class PCAInspector2Variables {

	public List<String[]> featureNames = new ArrayList<>();
	
	public double[][] extractAlignFeatures(File [] ascfeaturefiles, String nonNullexceptions []) throws Exception{
		
		double resolution = -1;
		double xll = 0;
		double yll = 0;
		int nrows = 0;
		int ncols = 0;
		double nodata = -9999;
		
		//1 column per feature = 1 feature array per row
		LinkedHashMap<Integer, Double[]> featurearray = new LinkedHashMap<>();
		int featureCounter = 0;
		int nFeatures = ascfeaturefiles.length;
		String [] fnames = new String[nFeatures];
		
		for (File f:ascfeaturefiles) {
			System.out.println("Reading file "+f.getName());
			
			fnames[featureCounter] = f.getName();
			
			AscRaster asc =  new AscRasterReader().readRaster(f.getAbsolutePath());
			
			if (resolution == -1) {
				resolution = asc.getCellsize();
				System.out.println("Spatial data resolution "+resolution);
				 xll = asc.getXll();
				 yll = asc.getYll();
				 nrows = asc.getRows();
				 ncols = asc.getCols();
				 nodata = Double.parseDouble(asc.getNDATA());
			}
			double [][] data = asc.getData();
			int featurerowcounter = 0;
			
			System.out.println("Populating feature array with "+(nrows*ncols)+" data and "+nFeatures+" features");
			
			for (int i=0;i<nrows;i++) {
				
				for (int j=0;j<ncols;j++) {
					Double [] dd = featurearray.get(featurerowcounter);
					if (dd==null) {
						dd = new Double[nFeatures];
					}
					double dat = data [i][j];
					for (String nonNull:nonNullexceptions) {
						if (f.getName().startsWith(nonNull) && (dat == nodata))
							dat = 0;
					}
						dd[featureCounter] = dat;
					featurearray.put(featurerowcounter,dd);
					
					featurerowcounter++;
				}
			}
			
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
		
		System.out.println("Removed "+removedLines+" over "+idxmax+"; Remaining "+(idxmax-removedLines)+" ("+Math.round((idxmax-removedLines)*100/idxmax)+"%)");
		ArrayList<Double[]> aslist = new ArrayList<Double[]>(featurearray.values());
		
		double [][] outputMatrix = new double [aslist.size()][nFeatures];
		for (int i=0;i<outputMatrix.length;i++) {
			
			for (int j=0;j<nFeatures;j++) {
				outputMatrix [i][j] = aslist.get(i)[j];
			}
			
		}
		
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
	
	public void PCACompare (File [] featureFiles) throws Exception{
		
		String nonNullexceptions[] = { "Sea_Ice_Concentration" };
		//one feature array for each row
		double [][] featureMatrix = extractAlignFeatures(featureFiles, nonNullexceptions);
		
		//run PCA
		Operations operations = new Operations();
		System.out.println("Standardizing Matrix");
		// standardize the matrix
		double[][] stdmatrix = operations.standardize(featureMatrix);
		System.out.println("CALCULATING PCA");
		//calculate PCA
		PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
		pca.calcPCA(stdmatrix);
		// get the pca components for all the vector
		double [][] featuresProjected = pca.getProjectionsMatrix(stdmatrix);
				
		//estimate the load from first eigen
		
		double e0 [] = pca.getEigenvector(0);
		double lambda0 = pca.getEigenvalue(0);
		
		double e1 [] = pca.getEigenvector(1);
		double lambda1 = pca.getEigenvalue(1);
		
		System.out.println("Eigenval X Eigenvector");
		System.out.println(Operations.roundDecimal(lambda0, 2)+" X "+Arrays.toString(e0));
		System.out.println(Operations.roundDecimal(lambda1, 2)+" X "+Arrays.toString(e1));
		
		double load0 [] = new double[e0.length];
		double load1 [] = new double[e1.length];
		for (int i=0;i<load0.length;i++) {
			
			load0 [i]=e0[i]*Math.sqrt(lambda0);
			load1 [i]=e1[i]*Math.sqrt(lambda1);
		}
		
		System.out.println("Loadings");
		System.out.println(Arrays.toString(load0));
		System.out.println(Arrays.toString(load1));
		
		 weight1 = (load0[0]+load1[0])/2;
		 weight2 = (load0[1]+load1[1])/2;
		double totalweight = Math.abs(weight1)+Math.abs(weight2);
		 percweight1 = weight1*100d/totalweight;
		 percweight2 = weight2*100d/totalweight;
		
		System.out.println("PCA loadings");
		System.out.println("w1: "+Operations.roundDecimal(weight1,2) + " ("+Operations.roundDecimal(percweight1,2)+"%)"+" - "+featureNames.get(0)[0]);
		System.out.println("w2: "+Operations.roundDecimal(weight2,2) + " ("+Operations.roundDecimal(percweight2,2)+"%)"+" - "+featureNames.get(0)[1]);

		//cut off one eigenvalue
		double [][] featuresProjectedReduced = new double[featuresProjected.length][1];
		for (int i=0;i<featuresProjected .length;i++) {
			featuresProjectedReduced [i][0] = featuresProjected[i][0];
		}

		//centroid in the projected space
		double centroid [][] = new double [1][1];
		double mean = 0;
		for (int i=0;i<featuresProjectedReduced .length;i++) {
			mean = mean + featuresProjectedReduced[i][0];
		}
		
		centroid [0] [0]= mean/(double)featuresProjectedReduced .length;
		//reproject back the centroid
		double [][] eigenMatrix = new double [e0.length][1];
		
		for (int i=0;i<eigenMatrix.length;i++) {
			eigenMatrix [i] [0]= e0[i];	
		}
		 
		double[][] centroidReprojection = Operations.multiplyMatrices(eigenMatrix,centroid);
		centroidReprojection = Operations.traspose(centroidReprojection);
		double[] centroidReprojected =  centroidReprojection[0];
		
		System.out.println("Reprojected centroid: ");
		System.out.println("c: "+Arrays.toString(centroidReprojected));
		
		double diffs1 [] = new double [featureMatrix.length];
		double diffs2 [] = new double [featureMatrix.length];
		
		for (int i=0;i<featureMatrix.length;i++) {
			diffs1 [i] = (featureMatrix[i][0] - centroidReprojected[0])*(featureMatrix[i][0] - centroidReprojected[0]);
			diffs2 [i] = (featureMatrix[i][1] - centroidReprojected[1])*(featureMatrix[i][1] - centroidReprojected[1]);
		}
		
		MSE1 = Operations.mean(diffs1)/(double)featureMatrix.length;
		MSE2 = Operations.mean(diffs2)/(double)featureMatrix.length;
		System.out.println("MSE1: "+Operations.roundDecimal(MSE1, 2)+" - "+featureNames.get(0)[0]);		
		System.out.println("MSE2: "+Operations.roundDecimal(MSE2, 2)+" - "+featureNames.get(0)[1]);
		
		double f1 [] = new double [featureMatrix.length];
		double f2 [] = new double [featureMatrix.length];
		
		for (int i=0;i<featureMatrix.length;i++) {
			f1 [i] = featureMatrix[i][0];
			f2 [i] = featureMatrix[i][1];
		}
		
		 sigma1 = Operations.variance(f1);
		 sigma2 = Operations.variance(f2);
		
		//this check demonstrates that an important aspect is also the interdependency between the variables - covariance
		System.out.println("Sigma1: "+Operations.roundDecimal(sigma1, 2)+" - "+featureNames.get(0)[0]);		
		System.out.println("Sigma2: "+Operations.roundDecimal(sigma2, 2)+" - "+featureNames.get(0)[1]);

		
	}
	public static void main(String[] args) throws Exception{
		String basepath = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS\\Adriatic_Sea\\2020\\";
		int casef = 9;
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
		
		PCAInspector2Variables inspector = new PCAInspector2Variables();
		inspector.PCACompare(featureFiles);
	}	
	
}
