package org.cnr.datanalysis.ecomod.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.rapidminer.example.ExampleSet;
import com.rapidminer.example.table.DataRow;
import com.rapidminer.example.table.ExampleTable;

public class Transformations {

	public static List<String> parseCSVString(String row, String delimiter) throws Exception {

		List<String> elements = new ArrayList<String>();
		String phrase = row;
		int idxdelim = -1;
		boolean quot = false;
		phrase = phrase.trim();
		while ((idxdelim = phrase.indexOf(delimiter)) >= 0) {
			quot = phrase.startsWith("\"");
			if (quot) {
				phrase = phrase.substring(1);
				String quoted = "";
				if (phrase.startsWith("\""))
					phrase = phrase.substring(1);
				else{
					Pattern pattern = Pattern.compile("[^\\\\]\"", Pattern.CASE_INSENSITIVE);
					Matcher regexp = pattern.matcher(phrase);
			    
					boolean matching = regexp.find();

					if (matching) {
						int i0 = regexp.start(0);
						quoted = phrase.substring(0, i0 + 1).trim();
						phrase = phrase.substring(i0 + 2).trim();
					}
				}

				if (phrase.startsWith(delimiter))
					phrase = phrase.substring(1);

				elements.add(quoted);

			} else {
				elements.add(phrase.substring(0, idxdelim));
				phrase = phrase.substring(idxdelim + 1).trim();
			}
		}
		if (phrase.startsWith("\""))
			phrase = phrase.substring(1);

		if (phrase.endsWith("\""))
			phrase = phrase.substring(0, phrase.length() - 1);

		elements.add(phrase);

		return elements;
	}
	
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
