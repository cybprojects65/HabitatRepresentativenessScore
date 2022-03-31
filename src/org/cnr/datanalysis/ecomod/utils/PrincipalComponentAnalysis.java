package org.cnr.datanalysis.ecomod.utils;

import com.rapidminer.example.ExampleSet;
import com.rapidminer.operator.IOContainer;
import com.rapidminer.operator.IOObject;
import com.rapidminer.operator.features.transformation.PCA;
import com.rapidminer.operator.features.transformation.PCAModel;
import com.rapidminer.tools.OperatorService;

public class PrincipalComponentAnalysis {

	static boolean RapidMinerInitialised = false; 
	public PrincipalComponentAnalysis(){
		if (!RapidMinerInitialised) {
			new RapidMinerConfiguration().initRapidMiner();
			RapidMinerInitialised = true;
		}
	}
	
	
	
	PCAModel innermodel;
	int numberOfComponents;
	
	public PCAModel getModel(){
		
		return innermodel;
	}
	
	public double[] getEigenvector (int index){
		
		return innermodel.getEigenvector(index);
	}
	
	public double getEigenvalue (int index){
		
		return innermodel.getEigenvalue(index);
	}
	
	
	public double [] getEigenvalues (){
		double [] values = new double[numberOfComponents];
		for (int i=0;i<numberOfComponents;i++){
			values[i] = getEigenvalue(i);
		}

		return values;
	}

	public double [] getNormalizedEigenvalues (){
		double [] values = new double[numberOfComponents];
		
		for (int i=0;i<numberOfComponents;i++){
			values[i] = getEigenvalue(i);
		}

		double sumEigen = Operations.sumVector(values);
		
		for (int i=0;i<numberOfComponents;i++){
			values[i] = values[i]/sumEigen;
		}
		
		return values;
	}

	
	public double [] getInverseEigenvalues (){
		double [] values = new double[numberOfComponents];
		for (int i=0;i<numberOfComponents;i++){
			values[i] = 1d/getEigenvalue(i);
		}
		return values;
	}

	public double [] getInverseNormalizedEigenvalues (){
		double [] values = new double[numberOfComponents];
		double[] weightedEigens = getNormalizedEigenvalues();
		for (int i=0;i<numberOfComponents;i++){
			values[i] = 1d/weightedEigens[i];
		}
		return values;
	}
	
	public double[][] getProjectionsMatrix(double[][] vectors) throws Exception{
		
		int nsamples=vectors.length;
		double[][] projected = new double[nsamples][];
		
		for (int i=0;i<nsamples;i++){
			projected[i] = getProjection(vectors[i]);
		}
		
		return projected;
	}

	
	public double[] getProjection(double[] vector) throws Exception{
		
		double [] projected = new double[numberOfComponents];
		for (int i=0;i<numberOfComponents;i++){
			projected[i] = Operations.scalarProduct(vector, getEigenvector(i));
		}
		return projected;
	}
	
	protected double[][] getPCA(double[][] sampleVectors) throws Exception{
		
		ExampleSet set = Transformations.matrix2ExampleSet(sampleVectors);
		ExampleSet outset = innermodel.apply(set);
		return Transformations.exampleSet2Matrix(outset);
		
	}
	
	public void calcPCA(double [][] sampleVectors) throws Exception{
		
		System.out.println("STARTING PCA COMPUTATION");
		
		PCA pca = (PCA) OperatorService.createOperator("PCA");
		pca.setParameter("variance_threshold", "0.95");
		//pca.setParameter("dimensionality_reduction", "keep variance");
		pca.setParameter("dimensionality_reduction", "none");
		pca.setParameter("number_of_components", "-1");
		
		ExampleSet set = Transformations.matrix2ExampleSet(sampleVectors);
		
		IOContainer innerInput = new IOContainer(set);
		IOContainer output = pca.apply(innerInput);
		IOObject[] outputvector = output.getIOObjects();
		innermodel = (PCAModel) outputvector[1];
		numberOfComponents = innermodel.getMaximumNumberOfComponents();
		System.out.println("MODEL APPLIED");
	}
	
	
}
