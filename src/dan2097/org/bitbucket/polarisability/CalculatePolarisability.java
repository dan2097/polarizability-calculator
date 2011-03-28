package dan2097.org.bitbucket.polarisability;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.ggasoftware.indigo.Indigo;
import com.ggasoftware.indigo.IndigoObject;

public class CalculatePolarisability {

	private static final double a0Cubed = 1.48184534*Math.pow(10,-31);
	private static Indigo indigo = new Indigo();
	
	public static void main(String[] args) throws IOException {
		System.err.println("Input is line seperated SMILES. Output format is tab delimited: Smiles, polarisability (with PEOE), polarisability (without PEOE)");
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		String smiles = stdinReader.readLine().trim();
		while(smiles !=null) {
			if (smiles.equals("") || smiles.startsWith("//")){
				System.out.println("");
			}
			else {
				Double polarisabilityWithoutPEOE = null;
				Double polarisabilityWithPEOE = null;
				try{
					polarisabilityWithoutPEOE = CalculatePolarisability.smilesToPolarisability(smiles, false);
					try{
						polarisabilityWithPEOE = CalculatePolarisability.smilesToPolarisability(smiles, true);
					}
					catch (Exception e) {
						System.err.println("Failed to calculate polarisability for: " + smiles +" (with PEOE)");
					}
				}
				catch (Exception e) {
					System.err.println("Failed to calculate polarisability for: " + smiles +" (without PEOE)");
				}
				System.out.println(smiles +"\t" + (polarisabilityWithPEOE != null ? polarisabilityWithPEOE : "")+"\t" + (polarisabilityWithoutPEOE != null ? polarisabilityWithoutPEOE : ""));
			}
			System.out.flush();
			smiles = stdinReader.readLine();
		}
	}
	
	/**
	 * Uses the method described in "A fast empirical method for the calculation of molecular polarizability"
	 * to rapidly calculate molecular polarizability. DOI: 10.1007/BF00125380
	 * 
	 * Results are ouput in units of 10^-25 cm3
	 * 
	 * PEOE (partial equalization of orbital electronegativity) should be enabled for the most accurate result
	 * @param smiles
	 * @param usePEOE
	 * @return
	 */
	public static double smilesToPolarisability(String smiles, boolean usePEOE){
		IndigoObject molecule = indigo.loadMolecule(smiles);
		return calculatePolarisability(molecule, usePEOE);
	}

	/**
	 * Calculates molecular polarisability in units of 10^-25 cm3
	 * @param molecule
	 */
	private static double calculatePolarisability(IndigoObject molecule, boolean peoe) {
		molecule.unfoldHydrogens();
		
		/*Maps between atom ids and partial charges*/
		Map<Integer, Double> partialCharges = new HashMap<Integer, Double>();
		for (Iterator<IndigoObject> iterator = molecule.iterateAtoms(); iterator.hasNext();) {
			IndigoObject atom = iterator.next();
			partialCharges.put(atom.index(), Double.valueOf(atom.charge()));//use the formal charge as a first approximation
		}
		if (peoe){
			determinePartialChargesUsingPEOE(molecule, partialCharges);
		}

		double moleculePolarisability = 0;
		for (Iterator<IndigoObject> iterator = molecule.iterateAtoms(); iterator.hasNext();) {
			IndigoObject atom = iterator.next();
			ElectronicConfiguration config = ElectronicConfiguration.getElectronicConfig(atom.atomicNumber(), partialCharges.get(atom.index()));
			Double polarisability = 4d/9d * determineSumOfMeanSquareRadiusSquared(config) * a0Cubed;
			moleculePolarisability += polarisability;
		}
		moleculePolarisability*=1000000;//convert to cm3
		return moleculePolarisability;
	}

	private static void determinePartialChargesUsingPEOE(IndigoObject molecule, Map<Integer, Double> partialCharges) {
		for (int i = 1; i <= 5; i++) {//perform 5 iterations
			for (Iterator<IndigoObject> iterator = molecule.iterateBonds(); iterator.hasNext();) {
				IndigoObject bond = iterator.next();
				IndigoObject sourceAtom = bond.source();
				IndigoObject destinationAtom = bond.destination();
				double parentAtomElectronegativity = determineElectronegativity(sourceAtom, partialCharges.get(sourceAtom.index()));
				double neighbourAtomElectronegativity = determineElectronegativity(destinationAtom, partialCharges.get(destinationAtom.index()));
				IndigoObject moreElectronegativeAtom = parentAtomElectronegativity >=neighbourAtomElectronegativity ? sourceAtom : destinationAtom;
				double moreElectroNegValue = parentAtomElectronegativity >=neighbourAtomElectronegativity ? parentAtomElectronegativity : neighbourAtomElectronegativity;
				IndigoObject lessElectronegativeAtom = parentAtomElectronegativity >= neighbourAtomElectronegativity ? destinationAtom : sourceAtom;
				double lessElectroNegValue = parentAtomElectronegativity >=neighbourAtomElectronegativity ? neighbourAtomElectronegativity : parentAtomElectronegativity;
				double lessElectroNegPositiveIonValue = determineElectronegativity(lessElectronegativeAtom, partialCharges.get(lessElectronegativeAtom.index())+1);
				double chargeTransferred = (( moreElectroNegValue - lessElectroNegValue)/lessElectroNegPositiveIonValue)*Math.pow(0.5, i);
				partialCharges.put(moreElectronegativeAtom.index(), partialCharges.get(moreElectronegativeAtom.index()) -chargeTransferred);//make the more electronegative atom more negative
				partialCharges.put(lessElectronegativeAtom.index(), partialCharges.get(lessElectronegativeAtom.index()) +chargeTransferred);
			}
		}
	}

	private static double determineElectronegativity(IndigoObject atom, Double partialCharge) {
		int neighbourCount =0;
		for (Iterator<IndigoObject> neighbourIterator = atom.iterateNeighbors(); neighbourIterator.hasNext();) {
			neighbourIterator.next();
			neighbourCount++;
		}
		int atomicNumber =atom.atomicNumber();
		if (atomicNumber==1){
			return 7.17 +6.24 * partialCharge + -0.56*Math.pow(partialCharge, 2);
		}
		if (atomicNumber==6){
			if (neighbourCount==4){
				return 7.98 + 9.18 * partialCharge + 1.88 * Math.pow(partialCharge, 2);
			}
			else if (neighbourCount==3){
				return 8.79 + 9.32 * partialCharge + 1.51 * Math.pow(partialCharge, 2);
			}
			else if (neighbourCount==2){
				return 10.39 + 9.45 * partialCharge + 0.73 * Math.pow(partialCharge, 2);
			}
			else{
				throw new IllegalArgumentException("Carbon in unusual valency");
			}
		}
		else if (atomicNumber==7){
			if (neighbourCount==3 && atom.valence()==3){
				return 11.54 + 10.82 * partialCharge + 1.36 * Math.pow(partialCharge, 2);
			}
			else if (neighbourCount==2 && atom.valence()==3){
				return 12.87 + 11.15 * partialCharge + 0.85 * Math.pow(partialCharge, 2);
			}
			else if (neighbourCount==1 && atom.valence()==3){
				return 15.68 + 11.7 * partialCharge + -0.27 * Math.pow(partialCharge, 2);
			}
			else{
				throw new IllegalArgumentException("Nitrogen in unusual valency");
			}
		}
		else if (atomicNumber==8){
			if (neighbourCount==2){
				return 14.18 + 12.92 * partialCharge + 1.39 * Math.pow(partialCharge, 2);
			}
			else if (neighbourCount==1){
				return 17.07 + 13.79 * partialCharge + 0.47 * Math.pow(partialCharge, 2);
			}
			else{
				throw new IllegalArgumentException("Oxygen in unusual valency");
			}
		}
		else if (atomicNumber==9){
			if (neighbourCount==1){
				return 14.66 + 13.85 * partialCharge + 2.31 * Math.pow(partialCharge, 2);
			}
			else{
				throw new IllegalArgumentException("Fluorine in unusual valency");
			}
		}
		else if (atomicNumber==16){
			return 10.14 + 9.13 * partialCharge + 1.38 * Math.pow(partialCharge, 2);
		}
		else if (atomicNumber==17){
			if (neighbourCount==1){
				return 11.00 + 9.69 * partialCharge +1.35 * Math.pow(partialCharge, 2);
			}
			else{
				throw new IllegalArgumentException("Chlorine in unusual valency");
			}
		}
		else if (atomicNumber==35){
			if (neighbourCount==1){
				return 10.08 + 8.47 * partialCharge + 1.16 * Math.pow(partialCharge, 2);
			}
			else{
				throw new IllegalArgumentException("Bromine in unusual valency");
			}
		}
		else if (atomicNumber==53){
			if (neighbourCount==1){
				return 9.90 + 7.96 * partialCharge +0.96 * Math.pow(partialCharge, 2);
			}
			else{
				throw new IllegalArgumentException("Iodine in unusual valency");
			}
		}
		else{
			throw new IllegalArgumentException("Non-organic atoms are not currently supported");
		}
	}

	private static Double determineSumOfMeanSquareRadiusSquared(ElectronicConfiguration config) {
		List<Map<SlaterShell, Double>> electronPopulations = config.getElectronPopulations();
		Double sumOfMeanSquareRadiusSquared = 0d;
		for (int i = electronPopulations.size() -1; i >=0; i--) {
			Map<SlaterShell, Double> shellToElectronCountMap = electronPopulations.get(i);
			for (SlaterShell subshell : shellToElectronCountMap.keySet()) {
				Double effectiveNuclearCharge = config.getAtomicNumber() - getScreeningConstantForSubshell(electronPopulations, i, subshell);
				Double effectivePrincipleQuantumNumber = getEffectivePrincipleQuantumNumber(i+1);
				Double prefactor = Math.pow(effectivePrincipleQuantumNumber/(2*effectiveNuclearCharge), 2);
				prefactor*=(2*effectivePrincipleQuantumNumber+1);
				prefactor*=(2*effectivePrincipleQuantumNumber+2);
				Double meanSquareRadius = prefactor;
				sumOfMeanSquareRadiusSquared += (shellToElectronCountMap.get(subshell) * Math.pow(meanSquareRadius, 2));
			}
		}
		return sumOfMeanSquareRadiusSquared;
	}

	/**
	 * Given the electronPopulations map, the principle equantum number (-1) and the subshell determines the screening constant
	 * @param electronPopulations
	 * @param principleQuantumNumberMinusOne
	 * @param subshell
	 * @return
	 */
	private static double getScreeningConstantForSubshell(List<Map<SlaterShell, Double>> electronPopulations, int principleQuantumNumberMinusOne, SlaterShell valentSubShell) {
		double screeningConstant =0;
		for (int i = principleQuantumNumberMinusOne; i >=0; i--) {
			Map<SlaterShell, Double> shellToElectronCountMap = electronPopulations.get(i);
			for (SlaterShell subshell : shellToElectronCountMap.keySet()) {
				double electronCount = shellToElectronCountMap.get(subshell);
				if (i==principleQuantumNumberMinusOne){//same principle quantum number
					if (subshell.equals(valentSubShell)){
						electronCount--;
						if (i==0){//1s is a special case
							screeningConstant += electronCount *0.3;
						}
						else{
							screeningConstant += electronCount *0.35;
						}
					}
					else if ((valentSubShell.equals(SlaterShell.d) && subshell.equals(SlaterShell.sp)) || 
								valentSubShell.equals(SlaterShell.f)){
						screeningConstant += electronCount;
					}
				}
				else if (i == principleQuantumNumberMinusOne -1 && valentSubShell.equals(SlaterShell.sp)){//one shell down 
					screeningConstant += electronCount *0.85;
				}
				else{
					screeningConstant += electronCount;//all other cases we assume the inner electrons shield exactly one charge each
				}
			}
		}
		return screeningConstant;
	}

	/**
	 * Given a principal quantum number returns n*, a value that empirically gave the best values
	 * @param i
	 * @return 
	 */
	private static double getEffectivePrincipleQuantumNumber(int i) {
		if(i==1){
			return 0.94;
		}
		else if(i==2){
			return 1.88;
		}
		else if(i==3){
			return 2.65;
		}
		else if(i==4){
			return 3.15;
		}
		else if(i==5){
			return 3.35;
		}
		else{
			throw new IllegalArgumentException("Only principal quantum numbers 1 through 5 are currently supported");
		}
	}
}
