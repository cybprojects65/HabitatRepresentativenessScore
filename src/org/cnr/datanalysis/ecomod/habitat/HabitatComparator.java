package org.cnr.datanalysis.ecomod.habitat;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import org.cnr.datanalysis.ecomod.utils.Operations;

import it.cnr.raster.asc.filemanagement.AscRaster;
import it.cnr.raster.asc.filemanagement.AscRasterReader;

public class HabitatComparator {

	public List<String[]> featureNames = new ArrayList<>();
	
	public double[][] extractAlignFeatures(File habitat, String nonNullexceptions []) throws Exception{
		
		File [] listFilesh1=habitat.listFiles();
		double resolution = -1;
		double xll = 0;
		double yll = 0;
		int nrows = 0;
		int ncols = 0;
		double nodata = -9999;
		//1 column per feature = 1 feature array per row
		LinkedHashMap<Integer, Double[]> featurearray = new LinkedHashMap<>();
		int featureCounter = 0;
		int nFeatures = listFilesh1.length;
		String [] fnames = new String[nFeatures];
		
		
		System.out.println("Analysing "+habitat.getAbsolutePath());
		for (File f:listFilesh1) {
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
		String fnames = Arrays.toString(ch.featureNames.get(0));
		System.out.println("Habitat vector score\n"+fnames+"\n"+Arrays.toString(hrs.HRS_VECTOR));
		System.out.println("HRS "+hrs.HRS);
	}
	
	
}
