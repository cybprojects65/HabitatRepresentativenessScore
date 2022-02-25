package org.cnr.datanalysis.ecomod.utils;

import com.rapidminer.RapidMiner;

public class RapidMinerConfiguration{

	
	public void initRapidMiner(){
		System.setProperty("rapidminer.init.operators", "./cfg/operators.xml");
		RapidMiner.init();
		System.out.println("Rapid Miner initialized");
	}
	

}
