package org.cnr.datanalysis.ecomod.utils;

import com.rapidminer.example.ExampleSet;
import com.rapidminer.example.table.DataRow;
import com.rapidminer.example.table.ExampleTable;

public class Transformations {

	public static ExampleSet matrix2ExampleSet(double[][] sampleVectors) {

		int m = sampleVectors.length;

		BigSamplesTable samples = new BigSamplesTable();

		for (int k = 0; k < m; k++)
			samples.addSampleRow("sample", sampleVectors[k]);

		return samples.generateExampleSet();

	}

	public static double[][] exampleSet2Matrix(ExampleSet set) {

		int m = set.size();
		ExampleTable table = set.getExampleTable();
		int n = table.getAttributeCount();
		double[][] matrix = new double[m][n - 1];
		for (int i = 0; i < m; i++) {
			DataRow row = table.getDataRow(i);
			for (int j = 0; j < n - 1; j++) {
				if (!table.getAttribute(j).isNominal()) {
					double d = row.get(table.getAttribute(j));
					matrix[i][j] = d;
				}
			}
		}

		return matrix;

	}


}
