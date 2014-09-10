package de.ipbhalle.converter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.codec.DecoderException;
import org.apache.commons.codec.binary.Hex;
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

import de.ipbhalle.wrapper.Container;

public class MassBankToLibrary {

	private enum InstType {IT, TQ, Q, TOF, ICR, FTMS, ESI};	// ESI == ESI-TOF
	private enum IoniMethod {EI, CI, APCI, ESI, nano, TS, MALDI, CAESIUM, APMALDI, APPI};	// nano == nano-ESI
	private enum IoniPolarity {neg, pos, undefined, both};		// neg == 0, pos == 1, undefined == -1, both == 2
	
	static final String HEXES = "0123456789ABCDEF";
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) {
		String sep = System.getProperty("file.separator");
		
		if(args == null || args.length == 0 || args.length != 3) {
			System.err.println("Missing arguments!");
			System.err.println("Usage: /path/to/spectra /path/to/moldata /output/path");
			System.exit(-1);
		}
		
		String wd = args[0];	//"EI+GC";
		String molPath = args[1];
		String outpath = args[2];
		
		String prevId = "";
		File dir = new File(wd);
		
		/**
		 * implement own FilenameFilter
		 */
		FilenameFilter ff = new FilenameFilter() {
			
			@Override
			public boolean accept(File dir, String name) {
				if(!name.endsWith("txt"))
					return false;
				else return true;
			}
		};
		
		File[] records = dir.listFiles(ff);
		System.out.println("#records -> " + records.length);
		for (File f : records) {
			prevId = f.getName().substring(0, f.getName().indexOf("."));
			String molfile = prevId + ".mol";
			
			String moldata = "";
			File mol = new File(molPath, molfile);
			if(!mol.exists()) {
				System.out.println("molfile [" + mol.getName() + "] does not exist! -> continuing.");
				//continue;
			}
			
			try {
				BufferedReader br = new BufferedReader(new FileReader(mol));
				StringBuffer sb = new StringBuffer();
				String line = "";
				while((line = br.readLine()) != null) {
					sb.append(line).append("\n");		// append windows like line separator
				}
				br.close();
				moldata = sb.toString().trim();
			} catch (FileNotFoundException e) {
				System.err.println("FileNotFound exception for [" + prevId + "]");
				moldata = "";
			} catch (IOException e) {
				System.err.println("IOException for [" + prevId + "]");
				moldata = "";
			}
			Container c = new Container(prevId, moldata, null, true);
			
			System.out.println("\tcreating spectrum file for [" + prevId + "]");
			boolean success;
			try {
				success = writeMassBankToSpectrum(wd + sep + prevId + ".txt", outpath, c);
				if(success)
					System.out.println("\tconversion successful for [" + prevId + "]");
				else System.err.println("\terror while converting [" + prevId + "]");
			} catch (IOException e) {
				System.err.println("IOExeption for [" + prevId + "] -> line 92");
			}
		}
	}

	public static boolean writeMassBankToSpectrum(String recordFile, String dir, Container c) throws IOException {
		File outDir = new File(dir);
		if(!outDir.exists())
			outDir.mkdir();
		
		File record = new File(recordFile);
		if(!record.exists()) {
			System.err.println("file " + record.getAbsolutePath() + " does not exist!");
			return false;
		}
		
		BufferedReader br = new BufferedReader(new FileReader(record));
		String accession = "", title = "", date = "", authors = "", copyright = "", formula = "", mass = "", smiles = "", 
		iupac = "", instrument = "", inst_type = "";
		int numPeaks = 0, peaksPerRow = 10 , nominal_mass = 0;
		List<String> names = new ArrayList<String>();
		List<String> links = new ArrayList<String>();
		List<String> acs = new ArrayList<String>();
		List<String> peaks = new ArrayList<String>();
		//String contributor = "";
		String colEnergy = "";
		String MSMS = "";
		String polarity = "";
		String IoniMethod = "";		// ionization method
		
		String name = "", nist = "", un = "", mw = "";
		String year = "", link = "", comment = "";
		String inchi = "";
		String preion = "", prodion = "", trapdrive = "", skim = "", fragampl = "", isolwidth = "", targetgas = "", targetgaspres = "",
			reagention = "", reagentgaspres = "", peakwidth = "", reflector = "", psd = "", chargedecon = "", column = "", rettime = "",
			ssid = "", analid = "", analname = "", massrange = "";
		String contributor = "";
		String instrument_type = "ESI-TOF";
		String compound_class_unknown = "CH$COMPOUND_CLASS: unknown";
		String compound_class = "";
		String ev = "", ion = "", precursor = "";
		String ion_method = "ESI";
		String ion_mode = "pos";
		int ms = 1;
		
		// use apache codec function to encode/decode hex
//		String structure = String.valueOf(Hex.encodeHex(c.getMoldata().getBytes()));
//		try {
//			System.out.println("\n Hex decode -> \n" + Hex.decodeHex(structure.toCharArray()));
//		} catch (DecoderException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		String structure = convertStringToHex(c.getMoldata());	// convert mol data into hex string
		if(c.getMoldata().isEmpty())
			structure = "";
		if(c != null && !c.getMoldata().isEmpty())
			structure = getHex(c.getMoldata().getBytes());
		
		String kegg = "", cid = "", cas = "";
		StringBuilder comments = new StringBuilder();	// adding all AC$ and MS$
		
		String line = "";
		while((line = br.readLine()) != null) {
			if(line.startsWith("ACCESSION:")) {
				accession = line.substring(line.indexOf(":") + 1).trim();
				names.add("MassBank: " + accession);
			}
			else if(line.startsWith("RECORD_TITLE:")) {
				String[] split = line.split(";");
				title = split[0].substring(split[0].indexOf(":") + 1).trim();
				
				if(line.contains("EI-") || split.length == 3 || split.length == 4) {	// EI spectrum
					title = line.substring(line.indexOf(":") + 1, line.indexOf(";")).trim();
					MSMS = "";
					inst_type = de.ipbhalle.converter.MassBankToLibrary.IoniMethod.EI.toString();
					continue;
				}
				
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
				//date = date + " 14:48:36+02:00";
				/** require standard date for Bruker LibraryEditor */
				date = "2010-10-09 19:48:36+02:00";	
			}
			else if(line.startsWith("AUTHORS:"))
				authors = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("COPYRIGHT:"))
				copyright = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$NAME:") || line.contains("CH$NAME:"))
				names.add(line.substring(line.indexOf(":") + 1).trim());
			else if(line.startsWith("CH$COMPOUND_CLASS:"))
				comments.append("Comment: " + line.substring(line.indexOf("$") + 1).trim() + "\n");
			else if(line.startsWith("CH$FORMULA:"))
				formula = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$EXACT_MASS:")) {
				mass = line.substring(line.indexOf(":") + 1).trim();
				//nominal_mass = mass.substring(0, mass.indexOf("."));
				Double d = Double.parseDouble(mass);
				nominal_mass = d.intValue();
			}
			else if(line.startsWith("CH$SMILES:"))
				smiles = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$IUPAC:"))
				iupac = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("CH$LINK:")) {
				links.add(line.substring(line.indexOf(":") + 1).trim());
				if(line.contains("KEGG"))
					//kegg = line.substring(line.indexOf("KEGG") + 4, line.indexOf("</a>")).trim();
					kegg = line.substring(line.indexOf("KEGG") + 4).trim();
					//kegg = "KEGG: " + line.substring(line.indexOf("_blank\">") + 8, line.indexOf("</a>")).trim();
				else if(line.contains("CID"))
					cid = line.substring(line.indexOf("CID:") + 4).trim();
					//cid = line.substring(line.indexOf("_blank\">") + 8, line.indexOf("</a>")).trim();
				else if(line.contains("PUBCHEM"))
					cid = line.substring(line.indexOf(":") + 1).trim();
					//cid = line.substring(line.indexOf("_blank\">") + 8, line.indexOf("</a>")).trim();
				else if(line.contains("CAS"))
					cas = line.substring(line.indexOf("CAS") + 3).trim();
					//cas = line.substring(line.indexOf("_blank\">") + 8, line.indexOf("</a>")).trim();
			}
			else if(line.startsWith("AC$INSTRUMENT:") || line.contains("AC$INSTRUMENT:"))
				instrument = line.substring(line.indexOf(":") + 1).trim();
			else if(line.startsWith("AC$ANALYTICAL_CONDITION:")) {
				acs.add(line.substring(line.indexOf("$") + 1).trim());
				if(line.contains("IONIZATION"))
					IoniMethod = line.substring(line.indexOf("IONIZATION") + 10).trim();
				else if(line.contains("COLLISION_ENERGY"))
					colEnergy = line.substring(line.indexOf("COLLISION_ENERGY") + 17, line.indexOf(" ", line.indexOf("COLLISION_ENERGY") + 17)).trim();
				else if(line.contains("MODE")) {
					String temp = line.substring(line.indexOf("MODE") + 4).trim();
					if(temp.equals("POSITIVE"))
						polarity = "pos";
					else polarity = "neg";
				}
				else if(line.contains("MS_TYPE") || line.contains("MS/MS")) {
					String temp = line.substring(line.indexOf("MS_TYPE") + 8).trim();
					if(temp.equals("MS/MS"))
						ms = 2;
					else if(temp.equals("MS"))
						ms = 1;
				}
				else if(line.contains("IsolWidth:"))
					isolwidth = line.substring(line.indexOf("IsolWidth:") + 11).trim();
				else if(line.contains("TrapDrive:"))
					trapdrive = line.substring(line.indexOf("TrapDrive:") + 11).trim();
				else if(line.contains("Skim1:"))
					skim = line.substring(line.indexOf("Skim1:") + 7).trim();
				else if(line.contains("FragAmpl:"))
					fragampl = line.substring(line.indexOf("FragAmpl:") + 10).trim();
				else if(line.contains("IsolWidth:"))
					isolwidth = line.substring(line.indexOf("IsolWidth:") + 11).trim();
				else if(line.contains("TargetGas:"))
					targetgas = line.substring(line.indexOf("TargetGas:") + 11).trim();
				else if(line.contains("TargetGasPres"))
					targetgaspres = line.substring(line.indexOf("TargetGasPres") + 14).trim();
				else if(line.contains("ReagentIon:"))
					reagention = line.substring(line.indexOf("ReagentIon:") + 12).trim();
				else if(line.contains("ReagentGasPres"))
					reagentgaspres = line.substring(line.indexOf("ReagentGasPres") + 15).trim();
				else if(line.contains("PeakWidth:"))
					peakwidth = line.substring(line.indexOf("PeakWidth:") + 11).trim();
				else if(line.contains("Refl"))
					reflector = line.substring(line.indexOf("Refl") + 5).trim();
				else if(line.contains("PSD:"))
					psd = line.substring(line.indexOf("PSD:") + 5).trim();
				else if(line.contains("ChargeDecon"))
					chargedecon = line.substring(line.indexOf("ChargeDecon") + 12).trim();
				else if(line.contains("Column:"))
					column = line.substring(line.indexOf("Column:") + 8).trim();
				else if(line.contains("RetTime:"))
					rettime = line.substring(line.indexOf("RetTime:") + 9).trim();
				else if(line.contains("SSID:"))
					ssid = line.substring(line.indexOf("SSID:") + 6).trim();
				else if(line.contains("AnalID:"))
					analid = line.substring(line.indexOf("AnalID:") + 8).trim();
				else if(line.contains("AnalName:"))
					analname = line.substring(line.indexOf("AnalName:") + 10).trim();
				else if(line.contains("Mass") | line.contains("MassRange"))
					massrange = line.substring(line.indexOf("Mass") + 5).trim();
			}
			else if(line.startsWith("PK$NUM_PEAK:") || line.contains("PK$NUM_PEAK:"))
				numPeaks = Integer.parseInt(line.substring(line.indexOf(":") + 1).trim());
			else if(line.startsWith("PK$PEAK:")) {
				for (int i = 0; i < numPeaks; i++) {
					line = br.readLine();
					String[] split = line.trim().split(" ");
					peaks.add(split[0] + " " + split[2]);
				}
			}
			
			if(line.contains("AC$") || line.contains("MS$"))
				comments.append("Comment: " + line.substring(line.indexOf("$") + 1).trim() + "\n");
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
			//title = title.replaceAll(",", "_");
			//title = title.replaceAll("\\s", "_");
		}
		
		/**
		 * for MB internal only
		 */
//		if(IoniMethod.isEmpty())
//			IoniMethod = "ESI";
		if(IoniMethod.isEmpty())	// for all EI spetra
			IoniMethod = "EI";
		
		/**
		 * write spectrum file
		 */
		//File spectrum = new File(dir, accession + "_" + title + ".msp");
		File spectrum = new File(dir, accession + ".msp");
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
			//StringBuffer sb = new StringBuffer();
			for (String s : names) {
				output.write("Synonym: " + s + "\n");	//sb.append(s + ", ");
			}
			//String temp = sb.toString().trim();
			//if(temp.endsWith(","))
			//	temp = temp.substring(0, temp.lastIndexOf(","));
			
			//output.write(temp + "\n");
		}
//		output.write("Comment: " + accession + "\n");
//		output.write(comments.toString());	// add all entries from AC$ and MS$
		
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
		if(preion != null && !preion.isEmpty())
			output.write("PreIon: " + preion + "\n");
		if(prodion != null && !prodion.isEmpty())
			output.write("ProdIon: " + prodion + "\n");
		if(trapdrive != null && !trapdrive.isEmpty())
			output.write("TrapDrive: " + trapdrive + "\n");
		if(skim != null && !skim.isEmpty())
			output.write("Skim1: " + skim + "\n");
		if(fragampl != null && !fragampl.isEmpty())
			output.write("FragAmpl: " + fragampl + "\n");
		if(isolwidth != null && !isolwidth.isEmpty())
			output.write("IsolWidth: " + isolwidth + "\n");
		if(targetgas != null && !targetgas.isEmpty())
			output.write("TargetGas: " + targetgas + "\n");
		if(targetgaspres != null && !targetgaspres.isEmpty())
			output.write("TargetGasPres: " + targetgaspres + "\n");
		if(reagention != null && !reagention.isEmpty())
			output.write("ReagentIon: " + reagention + "\n");
		if(reagentgaspres != null && !reagentgaspres.isEmpty())
			output.write("ReagentGasPres: " + reagentgaspres + "\n");
		if(colEnergy != null && !colEnergy.isEmpty())
			output.write("ColEnergy: " + colEnergy + "\n");
		if(peakwidth != null && !peakwidth.isEmpty())
			output.write("PeakWidth: " + peakwidth + "\n");
		if(reflector != null && !reflector.isEmpty())
			output.write("Refl: " + reflector + "\n");
		if(psd != null && !psd.isEmpty())
			output.write("PSD: " + psd + "\n");
		if(chargedecon != null && !chargedecon.isEmpty())
			output.write("ChargeDecon: " + chargedecon + "\n");
		if(date != null && !date.isEmpty())
			output.write("Date: " + date + "\n");
		if(column != null && !column.isEmpty())
			output.write("Column: " + column + "\n");
		if(rettime != null && !rettime.isEmpty())
			output.write("RetTime: " + rettime + "\n");
		if(ssid != null && !ssid.isEmpty())
			output.write("SSID: " + ssid + "\n");
		if(analid != null && !analid.isEmpty())
			output.write("AnalID: " + analid + "\n");
		if(analname != null && !analname.isEmpty())
			output.write("AnalName: " + analname + "\n");
		if(massrange != null && !massrange.isEmpty())
			output.write("Mass: " + massrange + "\n");
		output.write("Num Peaks: " + numPeaks + "\n");
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
		//String sep = System.getProperty("file.separator");
		FileWriter lib = new FileWriter(new File(dir, "MassBank_converted.library"), true);
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
			for(String s : names) {
				lib.write("Synonym: " + s + "\n");
			}
//			lib.write("Synonym: ");
//			StringBuffer sb = new StringBuffer();
//			for (String s : names) {
//				sb.append(s + ", ");
//			}
//			String temp = sb.toString().trim();
//			if(temp.endsWith(","))
//				temp = temp.substring(0, temp.lastIndexOf(","));
//			
//			lib.write(temp + "\n");
		}
//		lib.write("Comment: " + accession + "\n");
//		lib.write(comments.toString());	// add all entries from AC$ and MS$
		
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
		if(preion != null && !preion.isEmpty())
			lib.write("PreIon: " + preion + "\n");
		if(prodion != null && !prodion.isEmpty())
			lib.write("ProdIon: " + prodion + "\n");
		if(trapdrive != null && !trapdrive.isEmpty())
			lib.write("TrapDrive: " + trapdrive + "\n");
		if(skim != null && !skim.isEmpty())
			lib.write("Skim1: " + skim + "\n");
		if(fragampl != null && !fragampl.isEmpty())
			lib.write("FragAmpl: " + fragampl + "\n");
		if(isolwidth != null && !isolwidth.isEmpty())
			lib.write("IsolWidth: " + isolwidth + "\n");
		if(targetgas != null && !targetgas.isEmpty())
			lib.write("TargetGas: " + targetgas + "\n");
		if(targetgaspres != null && !targetgaspres.isEmpty())
			lib.write("TargetGasPres: " + targetgaspres + "\n");
		if(reagention != null && !reagention.isEmpty())
			lib.write("ReagentIon: " + reagention + "\n");
		if(reagentgaspres != null && !reagentgaspres.isEmpty())
			lib.write("ReagentGasPres: " + reagentgaspres + "\n");
		if(colEnergy != null && !colEnergy.isEmpty())
			lib.write("ColEnergy: " + colEnergy + "\n");
		if(peakwidth != null && !peakwidth.isEmpty())
			lib.write("PeakWidth: " + peakwidth + "\n");
		if(reflector != null && !reflector.isEmpty())
			lib.write("Refl: " + reflector + "\n");
		if(psd != null && !psd.isEmpty())
			lib.write("PSD: " + psd + "\n");
		if(chargedecon != null && !chargedecon.isEmpty())
			lib.write("ChargeDecon: " + chargedecon + "\n");
		if(date != null && !date.isEmpty())
			lib.write("Date: " + date + "\n");
		if(column != null && !column.isEmpty())
			lib.write("Column: " + column + "\n");
		if(rettime != null && !rettime.isEmpty())
			lib.write("RetTime: " + rettime + "\n");
		if(ssid != null && !ssid.isEmpty())
			lib.write("SSID: " + ssid + "\n");
		if(analid != null && !analid.isEmpty())
			lib.write("AnalID: " + analid + "\n");
		if(analname != null && !analname.isEmpty())
			lib.write("AnalName: " + analname + "\n");
		if(massrange != null && !massrange.isEmpty())
			lib.write("Mass: " + massrange + "\n");
		lib.write("Num Peaks: " + numPeaks + "\n");
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
}
