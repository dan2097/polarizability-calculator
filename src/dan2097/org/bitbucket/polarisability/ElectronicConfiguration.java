package dan2097.org.bitbucket.polarisability;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class ElectronicConfiguration {
	private ElectronicConfiguration(int atomicNumber, List<Map<SlaterShell,Double>> electronPopulations) {
		this.atomicNumber =atomicNumber;
		this.electronPopulations = electronPopulations;
	}

	/**
	 * A list for the principal quantum number
	 * Then a map with keys SlaterShell.sp, SlaterShell.d or SlaterShell.f and double indicating the number of electrons in that shell (can be fractional after application of PEOE)
	 */
	private final List<Map<SlaterShell,Double>> electronPopulations;
	
	private final int atomicNumber;
	
	int getAtomicNumber() {
		return atomicNumber;
	}

	List<Map<SlaterShell, Double>> getElectronPopulations() {
		return electronPopulations;
	}

	static Map<Integer, ElectronicConfiguration> mapToAtomElectronConfig = new HashMap<Integer, ElectronicConfiguration>();
	
	static{
		List<Map<SlaterShell,Double>> electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		Map<SlaterShell, Double> map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 1d);
		electronPopulations.add(map);
		ElectronicConfiguration config = new ElectronicConfiguration(1 , electronPopulations);
		mapToAtomElectronConfig.put(1, config);//hydrogen

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 3d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(5 , electronPopulations);
		mapToAtomElectronConfig.put(5, config);//boron
		
		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 4d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(6 , electronPopulations);
		mapToAtomElectronConfig.put(6, config);//carbon

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 5d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(7 , electronPopulations);
		mapToAtomElectronConfig.put(7, config);//nitrogen

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 6d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(8 , electronPopulations);
		mapToAtomElectronConfig.put(8, config);//oxygen
		

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 7d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(9 , electronPopulations);
		mapToAtomElectronConfig.put(9, config);//fluorine

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 4d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(14 , electronPopulations);
		mapToAtomElectronConfig.put(14, config);//silicon
		
		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 5d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(15 , electronPopulations);
		mapToAtomElectronConfig.put(15, config);//phosphorus

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 6d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(16 , electronPopulations);
		mapToAtomElectronConfig.put(16, config);//sulfur

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 7d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(17 , electronPopulations);
		mapToAtomElectronConfig.put(17, config);//chlorine

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		map.put(SlaterShell.d, 10d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 7d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(35 , electronPopulations);
		mapToAtomElectronConfig.put(35, config);//bromine

		electronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 2d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		map.put(SlaterShell.d, 10d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 8d);
		map.put(SlaterShell.d, 10d);
		electronPopulations.add(map);
		map  = new HashMap<SlaterShell, Double>();
		map.put(SlaterShell.sp, 7d);
		electronPopulations.add(map);
		config = new ElectronicConfiguration(53 , electronPopulations);
		mapToAtomElectronConfig.put(53, config);//iodine
	}

	/**
	 * Given an atomicNumber generates its electronicConfig.
	 * Null is returned if the atomicNumber given does not correspond to an atom with encoded electronic structure
	 * @param atomicNumber
	 * @param partialCharge 
	 * @return
	 */
	public static ElectronicConfiguration getElectronicConfig(int atomicNumber, Double partialCharge) {
		ElectronicConfiguration config = mapToAtomElectronConfig.get(atomicNumber);
		if (config == null){
			throw new IllegalArgumentException("The electronic configuration of an atom with atomic number: " + atomicNumber + " is unknown");
		}
		List<Map<SlaterShell,Double>> electronPopulations = config.getElectronPopulations();
		List<Map<SlaterShell,Double>> clonedElectronPopulations = new ArrayList<Map<SlaterShell, Double>>();
		for (Map<SlaterShell, Double> shellToElectronCountMap : electronPopulations) {
			Map<SlaterShell, Double> shellToElectronCountMapClone = new HashMap<SlaterShell, Double>();
			for (Entry<SlaterShell, Double> entry : shellToElectronCountMap.entrySet()) {
				shellToElectronCountMapClone.put(entry.getKey(), entry.getValue());
			}
			clonedElectronPopulations.add(shellToElectronCountMapClone);
		}
		adjustValenceShellElectronOccupancy(clonedElectronPopulations, partialCharge);
		return new ElectronicConfiguration(config.getAtomicNumber(), clonedElectronPopulations);
	}

	private static void adjustValenceShellElectronOccupancy(List<Map<SlaterShell, Double>> clonedElectronPopulations, Double partialCharge) {
		Map<SlaterShell, Double> shellToElectronCountMap = clonedElectronPopulations.get(clonedElectronPopulations.size()-1);
		if (shellToElectronCountMap.get(SlaterShell.f)!=null){
			shellToElectronCountMap.put(SlaterShell.f, shellToElectronCountMap.get(SlaterShell.f) - partialCharge);
		}
		else if (shellToElectronCountMap.get(SlaterShell.d)!=null){
			shellToElectronCountMap.put(SlaterShell.d, shellToElectronCountMap.get(SlaterShell.d) - partialCharge);
		}
		else {
			shellToElectronCountMap.put(SlaterShell.sp, shellToElectronCountMap.get(SlaterShell.sp) - partialCharge);
		}
	}
}
