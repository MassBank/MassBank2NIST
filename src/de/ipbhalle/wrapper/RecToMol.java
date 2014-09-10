package de.ipbhalle.wrapper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import com.sun.xml.internal.messaging.saaj.util.Base64;

public class RecToMol {

	private String id;
	private String name;
	private String file;
	
	public  enum InstType {IT, TQ, Q, TOF, ICR, FTMS, ESI};	// ESI == ESI-TOF
	public enum IoniMethod {EI, CI, APCI, ESI, nano, TS, MALDI, CAESIUM, APMALDI, APPI};	// nano == nano-ESI
	public enum IoniPolarity {neg, pos, undefined, both};		// neg == 0, pos == 1, undefined == -1, both == 2
	
	private static final String outDir = "MBinternal";
	static final String HEXES = "0123456789ABCDEF";

	public RecToMol(String id, String name, String file) {
		this.id = id;
		this.name = name;
		this.file = file;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File rectomol = new File("/home/mgerlich/mbinternal.txt");
		FileReader fr = new FileReader(rectomol);
		BufferedReader br = new BufferedReader(fr);
		String line = "";
		
		String prevId = "";
		String prevFile = "";
		String prevName = "";
		List<String> molfiles = new ArrayList<String>();
		
		String sep = System.getProperty("file.separator");
		String outputDir = System.getProperty("user.dir");
		outputDir += sep + outDir;
		File dir = new File(outputDir);
		dir.mkdir();
		
		while((line = br.readLine()) != null) {
			String[] split = line.split("\t");	// id	name	file
			if(prevId.equals(split[0]) && !prevFile.equals(split[2])) {
				molfiles.add("/vol/www/msbi/DocumentRoot/MassBankInternal/DB/molfile/MassBankInternal/" + split[2] + ".mol");
			}
			else if(prevId.isEmpty() && prevFile.isEmpty()) {	// first compound
				prevId = split[0];
				prevName = split[1];
				prevName = prevName.replaceAll("\\s", "_");
				prevFile = split[2];
				molfiles.add("/vol/www/msbi/DocumentRoot/MassBankInternal/DB/molfile/MassBankInternal/" + split[2] + ".mol");
			}
			else {
				Container compare = compareMols(molfiles, prevId);
				if(compare.isIsomorph()) {	// all molfiles are equal - write out corresponding .spectrum file with Base64 encoded molfile
					System.out.println("\tcreating spectrum file for [" + prevId + "]");
					//File spectrum = new File(outputDir, prevId + "_" + prevName + ".spectrum");
					//if(!spectrum.exists()) {
						boolean success = writeMassBankToSpectrum("/vol/www/msbi/DocumentRoot/MassBankInternal/DB/annotation/MassBankInternal/" + prevId + ".txt",
								outputDir, compare);
						if(success)
							System.out.println("\tconversion successful for [" + prevId + "]");
						else System.err.println("\terror while converting [" + prevId + "]");
					//}
					//else System.out.println("\tspectrum file already present -> " + spectrum.getAbsolutePath());	
					
					prevId = split[0];
					prevName = split[1];
					prevFile = split[2];
					molfiles = new ArrayList<String>();
				}
				else {			// molfiles were different
					for (String s : molfiles) {
						System.out.print(s + "  ");
					}
					System.err.println("\nmolfiles were different for id " + prevId);
					prevId = split[0];
					prevName = split[1];
					prevFile = split[2];
					molfiles = new ArrayList<String>();
				}
			}
		}
		br.close();
	}

	public static boolean writeMassBankToSpectrum(String recordFile, String dir, Container c) throws IOException {
		File record = new File(recordFile);
		if(!record.exists()) {
			System.err.println("file " + record.getAbsolutePath() + " does not exist!");
			return false;
		}
		
		BufferedReader br = new BufferedReader(new FileReader(record));
		String accession = "", title = "", date = "", authors = "", copyright = "", formula = "", mass = "", nominal_mass = "", 
			smiles = "", iupac = "", instrument = "", inst_type = "";
		int numPeaks = 0, peaksPerRow = 10;
		List<String> names = new ArrayList<String>();
		List<String> links = new ArrayList<String>();
		List<String> acs = new ArrayList<String>();
		List<String> peaks = new ArrayList<String>();
		//String contributor = "";
		String colEnergy = "";
		String MSMS = "";
		String polarity = "";
		String IoniMethod = "";		// ionization method
		//String structure = convertStringToHex(c.getMoldata());	// convert mol data into hex string
		String structure = getHex(c.getMoldata().getBytes());
		
		String kegg = "", cid = "", cas = "";
		
		String line = "";
		while((line = br.readLine()) != null) {
			if(line.startsWith("ACCESSION:"))
				accession = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("RECORD_TITLE:")) {
				String[] split = line.split(";");
				title = split[0].substring(split[0].indexOf(":") + 1).trim();
				
				if(split[1].trim().equals("MS/MS"))
					MSMS = "2";
				else MSMS = "";
				
				if(split[2].trim().contains("ESI-QTOF")) {
					inst_type = InstType.ESI.toString();
				}
				else if(split[2].trim().contains("ESI")) {
					inst_type = InstType.ESI.toString();
				}
				else if(split[2].trim().contains("IT")) {
					inst_type = InstType.IT.toString();
				}
				else if(split[2].trim().contains("FT")) {
					inst_type = InstType.FTMS.toString();
				}
				else if(split[2].trim().equals("QqTOF")) {
					inst_type = InstType.TOF.toString();
				}
				else if(split[2].trim().equals("QqQ")) {
					inst_type = InstType.TQ.toString();
				}
				else inst_type = "";
				
				if(split[3].contains("CE")) {
					split[3] = split[3].trim();
					colEnergy = split[3].substring(split[3].indexOf(":") + 1, split[3].lastIndexOf(" "));
				}
				
				if(split[4].trim().contains("+"))
					polarity = IoniPolarity.pos.toString();
				else if(split[4].trim().contains("-"))
					polarity = IoniPolarity.neg.toString();
				else polarity = "-1";
			}
			else if(line.startsWith("DATE:")) {
				date = line.substring(line.indexOf(":") + 1).trim();
				date = date.replaceAll("\\.", "-");
				date = date + " 14:48:36+02:00";
			}
			else if(line.startsWith("AUTHORS:"))
				authors = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("COPYRIGHT:"))
				copyright = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$NAME:"))
				names.add(line.substring(line.indexOf(":") + 1).trim());
			else if(line.startsWith("CH$FORMULA:"))
				formula = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$EXACT_MASS:")) {
				mass = line.substring(line.indexOf(":") + 1).trim();
				nominal_mass = mass.substring(0, mass.indexOf("."));
			}
			else if(line.startsWith("CH$SMILES:"))
				smiles = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$IUPAC:"))
				iupac = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$LINK:")) {
				links.add(line.substring(line.indexOf(":") + 1).trim());
				if(line.contains("KEGG"))
					kegg = line.substring(line.indexOf("KEGG") + 4).trim();
				else if(line.contains("CID"))
					cid = line.substring(line.indexOf("CID:") + 4).trim();
				else if(line.contains("CAS"))
					cas = line.substring(line.indexOf("CAS") + 3).trim();
			}
			else if(line.startsWith("AC$INSTRUMENT:"))
				instrument = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("AC$ANALYTICAL_CONDITION:")) {
				acs.add(line.substring(line.indexOf(":") + 1).trim());
				if(line.contains("IONIZATION"))
					IoniMethod = line.substring(line.indexOf("IONIZATION") + 10).trim();
				if(line.contains("COLLISION_ENERGY"))
					colEnergy = line.substring(line.indexOf("COLLISION_ENERGY") + 17, line.indexOf(" ", line.indexOf("COLLISION_ENERGY") + 17)).trim();
				if(line.contains("MODE")) {
					String temp = line.substring(line.indexOf("MODE") + 4).trim();
					if(temp.equals("POSITIVE"))
						polarity = "pos";
					else polarity = "neg";
				}
					
			}
			else if(line.startsWith("PK$NUM_PEAK:"))
				numPeaks = Integer.parseInt(line.substring(line.indexOf(":") + 1).trim());
			else if(line.startsWith("PK$PEAK:")) {
				for (int i = 0; i < numPeaks; i++) {
					line = br.readLine();
					String[] split = line.trim().split(" ");
					peaks.add(split[0] + " " + split[2]);
				}
			}
				
		}
		br.close();
		
		boolean simpleName = false;
		for (int i = 0; i < names.size(); i++) {
			if(title.contains(",") && !names.get(i).contains(",")) {
				title = names.get(i);
				simpleName = true;
				break;
			}
		}
		if(!simpleName) {	// no simple/trivial name found - convert special characters to _
			title = title.replaceAll(",", "_");
			title = title.replaceAll("\\s", "_");
		}
		
		/**
		 * for MB internal only
		 */
		if(IoniMethod.isEmpty())
			IoniMethod = "ESI";
		
		/**
		 * write spectrum file
		 */
		File spectrum = new File(dir, accession + "_" + title + ".spectrum");
		if(spectrum.exists()) {
			System.out.println("\tspectrum file already present -> " + spectrum.getAbsolutePath());	
			return true;
		}
		
		FileWriter output = new FileWriter(spectrum);
		output.write("Name: " + title + "\n");
		if(!cas.isEmpty())
			output.write("CAS: " + cas + "\n");
		output.write("MW: " + nominal_mass + "\n");
		if(formula != null && !formula.isEmpty())
			output.write("Formula: " + formula + "\n");
		if(!kegg.isEmpty())
			output.write("Synonym: " + kegg + "\n");
		if(!cid.isEmpty())
			output.write("Synonym: CID: " + cid + "\n");
		if(names != null && names.size() > 0) {
			output.write("Synonym: ");
			StringBuffer sb = new StringBuffer();
			for (String s : names) {
				sb.append(s + ", ");
			}
			String temp = sb.toString();
			if(temp.endsWith(", "))
				temp = temp.substring(0, temp.length() - 2);
			
			output.write(temp + "\n");
		}
		output.write("Comment: " + accession + "\n");
		if(structure != null && !structure.isEmpty())
			output.write("Structure: " + structure + "\n");
		output.write("Contributor: " + authors + "\n");
		if(inst_type != null && !inst_type.isEmpty())
			output.write("InstType: " + inst_type + "\n");
		if(instrument != null && !instrument.isEmpty())
			output.write("InstName: " + instrument + "\n");
		if(IoniMethod != null && !IoniMethod.isEmpty())
			output.write("IoniMethod: " + IoniMethod + "\n");
		if(polarity != null && !polarity.isEmpty())
			output.write("IonPolarity: " + polarity + "\n");
		if(MSMS != null && !MSMS.isEmpty())
			output.write("MSMS: " + MSMS + "\n");
		if(colEnergy != null && !colEnergy.isEmpty())
			output.write("ColEnergy: " + colEnergy + "\n");
		if(date != null && !date.isEmpty())
			output.write("Date: " + date + "\n");
		output.write("Num: " + numPeaks + "\n");
		for (int i = 0; i < peaks.size(); i++) {
			if(i != 0 && i % peaksPerRow == 0) {
				output.write("\n" + peaks.get(i) + " ");
			}
			else {
				output.write(peaks.get(i) + " ");
			}
		}
		output.write("\n");
		
		output.flush();
		output.close();
		
		
		// write/append library file
		String sep = System.getProperty("file.separator");
		String outputDir = System.getProperty("user.dir");
		outputDir += sep + outDir;
		FileWriter lib = new FileWriter(new File(outputDir, "MassBankInternal.library"), true);
		lib.write("Name: " + title + "\n");
		if(!cas.isEmpty())
			lib.write("CAS: " + cas + "\n");
		lib.write("MW: " + nominal_mass + "\n");
		lib.write("Formula: " + formula + "\n");
		if(!kegg.isEmpty())
			lib.write("Synonym: " + kegg + "\n");
		if(!cid.isEmpty())
			lib.write("Synonym: CID: " + cid + "\n");
		if(names != null && names.size() > 0) {
			lib.write("Synonym: ");
			StringBuffer sb = new StringBuffer();
			for (String s : names) {
				sb.append(s + ", ");
			}
			String temp = sb.toString();
			if(temp.endsWith(","))
				temp = temp.substring(0, temp.length() - 1);
			
			lib.write(temp + "\n");
		}
		lib.write("Comment: " + accession + "\n");
		if(structure != null && !structure.isEmpty())
			lib.write("Structure: " + structure + "\n");
		lib.write("Contributor: " + authors + "\n");
		if(inst_type != null && !inst_type.isEmpty())
			lib.write("InstType: " + inst_type + "\n");
		if(instrument != null && !instrument.isEmpty())
			lib.write("InstName: " + instrument + "\n");
		if(IoniMethod != null && !IoniMethod.isEmpty())
			lib.write("IoniMethod: " + IoniMethod + "\n");
		if(polarity != null && !polarity.isEmpty())
			lib.write("IonPolarity: " + polarity + "\n");
		if(MSMS != null && !MSMS.isEmpty())
			lib.write("MSMS: " + MSMS + "\n");
		if(colEnergy != null && !colEnergy.isEmpty())
			lib.write("ColEnergy: " + colEnergy + "\n");
		if(date != null && !date.isEmpty())
			lib.write("Date: " + date + "\n");
		lib.write("Num: " + numPeaks + "\n");
		for (int i = 0; i < peaks.size(); i++) {
			if(i != 0 && i % peaksPerRow == 0) {
				lib.write("\n" + peaks.get(i) + " ");
			}
			else {
				lib.write(peaks.get(i) + " ");
			}
		}
		lib.write("\n\n");	// end line of peaks and add newline to mark end of record
		
		lib.flush();
		lib.close();
	
		return true;
	}
	
	public static String convertStringToHex(String str) {

		char[] chars = str.toCharArray();

		StringBuffer hex = new StringBuffer();
		for (int i = 0; i < chars.length; i++) {
			hex.append(Integer.toHexString((int) chars[i]));
		}

		return hex.toString();
	}

	public static String convertHexToString(String hex) {

		StringBuilder sb = new StringBuilder();
		StringBuilder temp = new StringBuilder();

		// 49204c6f7665204a617661 split into two characters 49, 20, 4c...
		for (int i = 0; i < hex.length() - 1; i += 2) {

			// grab the hex in pairs
			String output = hex.substring(i, (i + 2));
			// convert hex to decimal
			int decimal = Integer.parseInt(output, 16);
			// convert the decimal to character
			sb.append((char) decimal);

			temp.append(decimal);
		}
		//System.out.println("Decimal : " + temp.toString());

		return sb.toString();
	}
	
	public static String getHex( byte [] raw ) {
	    if ( raw == null ) {
	      return null;
	    }
	    final StringBuilder hex = new StringBuilder( 2 * raw.length );
	    for ( final byte b : raw ) {
	      hex.append(HEXES.charAt((b & 0xF0) >> 4))
	         .append(HEXES.charAt((b & 0x0F)));
	    }
	    return hex.toString();
	  }
	
	public static String encodeBase64(String filename) throws IOException {
		File f = new File(filename);
		FileReader fr = new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		String line = "";
		StringBuffer sb = new StringBuffer();
		while((line = br.readLine()) != null) {
			sb.append(line);
		}
		br.close();
		fr.close();
		
		String toEncode = sb.toString();
		byte[] encoded = Base64.encode(toEncode.getBytes());
		String encodedString = new String(encoded);
		System.out.println(encodedString);
		
		return encodedString;
	}
	
	public static double computeNominalMass(double exactMass) {
		String temp = String.valueOf(exactMass);
		temp = temp.substring(0, temp.indexOf("."));
		return Double.parseDouble(temp);
	}
	
	public static double computeNaturalMass(String formula) {
		IMolecularFormula mf = MolecularFormulaManipulator.getMolecularFormula(formula, NoNotificationChemObjectBuilder.getInstance());
		return MolecularFormulaManipulator.getNaturalExactMass(mf);
	}
	
	public static double computeExactMass(String formula) {
		IMolecularFormula mf = MolecularFormulaManipulator.getMolecularFormula(formula, NoNotificationChemObjectBuilder.getInstance());
		return MolecularFormulaManipulator.getTotalExactMass(mf);
	}
	
	public static Container compareMols(List<String> molfiles, String id) {
		List<IAtomContainer> containers = new ArrayList<IAtomContainer>();
		
		//String linesep = System.getProperty("line.separator");
		String moldata = "";
		boolean wasread = false;
		for (String string : molfiles) {
			File f = new File(string);
			if(!f.exists())
				continue;
			
			if(!wasread) {
				try {
					BufferedReader br = new BufferedReader(new FileReader(f));
					StringBuffer sb = new StringBuffer();
					String line = "";
					while((line = br.readLine()) != null) {
						sb.append(line).append("\r\n");		// append windows like line separator
					}
					br.close();
					moldata = sb.toString();
					wasread = true;
				} catch (FileNotFoundException e) {
					wasread = false;
				} catch (IOException e) {
					wasread = false;
				}
			}
			
			FileInputStream fis = null;
			try {
				fis = new FileInputStream(f);
			} catch (FileNotFoundException e1) {
				System.err.println("File not found - " + f);
			}
			InputStream is = fis;
			MDLReader reader = new MDLReader(is);
			IChemFile chemFile = new ChemFile();
			IAtomContainer container = null;
			try {			
				chemFile = (IChemFile) reader.read(chemFile);
				container = ChemFileManipulator.getAllAtomContainers(chemFile).get(0);
				containers.add(container);	// add container to list
				
				// remove hydrogens
				//container = AtomContainerManipulator.removeHydrogens(container);
			} catch (java.lang.NumberFormatException e) {
				System.err.println("NumberFormatException occured while parsing mol file - " + f);
			} catch (CDKException e) {
				System.err.println("CDKException occured for mol file - " + f);
			}
		}
		
		if(containers.size() != molfiles.size())
			System.err.println("sizes of files [" + molfiles.size() + "] and containers [" + containers.size() + "] varies!");
		
		boolean check = true;
		if(containers.size() > 1) {		// only perform check if more than one molfile available
			for (int i = 1; i < containers.size(); i++) {
				try {
					check = UniversalIsomorphismTester.isIsomorph(containers.get(0), containers.get(i));
				} catch (CDKException e) {
					System.err.println("error for " + id);
				} 
			}
		}
		System.out.println("are isomorph -> " + check);
		
		Container c =  new Container(id, moldata, containers, check);
		return c;		// check
	}
	
	public void setId(String id) {
		this.id = id;
	}


	public String getId() {
		return id;
	}


	public void setName(String name) {
		this.name = name;
	}


	public String getName() {
		return name;
	}


	public void setFile(String file) {
		this.file = file;
	}


	public String getFile() {
		return file;
	}

}
