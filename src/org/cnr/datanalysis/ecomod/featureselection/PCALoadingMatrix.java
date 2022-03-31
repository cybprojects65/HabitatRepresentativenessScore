package org.cnr.datanalysis.ecomod.featureselection;

import java.io.File;
import java.io.FileWriter;

import org.cnr.datanalysis.ecomod.test.HabitatComparisonAllEuropeanSeas;
import org.cnr.datanalysis.ecomod.utils.Operations;

public class PCALoadingMatrix extends PCAInspector{

	
	public static void main(String[] args) throws Exception{
		
		String basepath = "C:\\Users\\Utente\\Ricerca\\Experiments\\Q-Quatics Climatic and AquaMaps data\\HRS input data\\Aegean_Sea\\2020\\";
		
		int year = 2020;
		/*
		String sea = "Adriatic_Sea";
		
		//+1 is assigned is weight(row) < weight(col) -1 otherwise; 0 if equal
		String featureList []= {
				"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc",
				"Sea-surface_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc",
				"Net_Primary_Production_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc",
				"Sea_Ice_Concentration_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc",
				"Sea-bottom_dissolved_oxygen_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc",
				"Sea-surface_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc",
				"Sea-bottom_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Adriatic_Sea.asc"
		};
		*/
		String sea = "Aegean_Sea";
		
		//+1 is assigned is weight(row) < weight(col) -1 otherwise; 0 if equal
		String featureList []= {
				"Sea-bottom_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Aegean_Sea.asc",
				"Sea-surface_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Aegean_Sea.asc",
				"Net_Primary_Production_res_01_annual_years_2020_Clim_scen_today_regional_Aegean_Sea.asc",
				"Sea_Ice_Concentration_res_01_annual_years_2020_Clim_scen_today_regional_Aegean_Sea.asc",
				"Sea-bottom_dissolved_oxygen_res_01_annual_years_2020_Clim_scen_today_regional_Aegean_Sea.asc",
				"Sea-surface_salinity_res_01_annual_years_2020_Clim_scen_today_regional_Aegean_Sea.asc",
				"Sea-bottom_temperature_res_01_annual_years_2020_Clim_scen_today_regional_Aegean_Sea.asc"
		};
		
		int [][] bonus = new int[featureList.length][featureList.length];
		
		for (int i=0;i<featureList.length;i++) {
			
			for (int j=0;j<featureList.length;j++) {
				System.out.println("Comparing "+featureList[i]+" vs "+featureList[j]);
				File [] featureFiles = {
						new File(basepath+featureList[i]),
						new File(basepath+featureList[j]),
				};
				PCAInspector inspector = new PCAInspector();
				inspector.PCACompare(featureFiles);
				int bon = 0;
				if (Operations.roundDecimal(inspector.weight1,2)>Operations.roundDecimal(inspector.weight2,2))
					bon = -1;
				else if (Operations.roundDecimal(inspector.weight1,2)<Operations.roundDecimal(inspector.weight2,2))
					bon = 1;
				
				bonus [i][j] = bon;
				System.out.println("BONUS "+bon+" w1 "+inspector.weight1+" vs w2 "+inspector.weight2);
			}
		}
		
		
		
		StringBuffer sb = new StringBuffer();
		sb.append(",");
		String fnames [] = HabitatComparisonAllEuropeanSeas.cleanupStrings(featureList);
		for (int i=0;i<bonus.length;i++) {
			sb.append(fnames[i]);
			if (i<bonus.length-1)
				sb.append(",");
			else
				sb.append("\n");
		}
		
		for (int i=0;i<bonus.length;i++) {
			sb.append(fnames[i]+",");
			
			for (int j=0;j<bonus.length;j++) {
				sb.append(bonus[i][j]);
				if (j<bonus.length-1)
					sb.append(",");
				else
					sb.append("\n");
			}	
			
		}
		
		String outFile = new File(basepath).getParentFile().getAbsolutePath()+"\\"+"PCALoadings_y"+year+"_"+sea+".csv";
		File outFileF = new File(outFile);
		FileWriter fw = new FileWriter(outFileF);
		fw.write(sb.toString());
		fw.close();
		
	}	
	
}
