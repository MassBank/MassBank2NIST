package de.ipbhalle.converter;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import net.sf.jniinchi.INCHI_RET;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smsd.algorithm.mcsplus.ExactMapping;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

// TODO: Auto-generated Javadoc
/**
 * The Class LibraryToMassBank.
 */
public class LibraryToMassBank {

	/** The prefix. */
	private String prefix = "XX";
	
	/** The output path. */
	private String outputPath = ""; //"/home/mgerlich/Desktop/bruker/";
	
	/** The mol path. */
	private String molPath = ""; 	//"/home/mgerlich/Desktop/bruker/moldata/";
	
	/** The rec path. */
	private String recPath = "";	//"/home/mgerlich/Desktop/bruker/recdata";
	
	/** The list mol. */
	File listMol;					// = new File("/home/mgerlich/Desktop/bruker/moldata/list.tsv");
	
	Properties prop;
	
	// These keys are present in the Bruker Library format
	private final String KEY_NAME = "Name:";
	private final String KEY_CAS = "CAS:";
	private final String KEY_PUBCHEM= "PubChem:";
	private final String KEY_NIST = "NIST:";
	private final String KEY_UN = "UN:";
	private final String KEY_MW = "MW:";
	private final String KEY_INCHI = "InChI:";
	private final String KEY_INCHIKEY = "InChIKey:";
	private final String KEY_FORMULA = "Form";	// allow both Form(ula)
	private final String KEY_SYNONYM = "Syn";	// allow both Syn(onym)
	private final String KEY_COMMENT = "Com";	// allow both Com(ment)
	private final String KEY_STRUCTURE = "Struc";	// allow both Struc(ture)
	private final String KEY_CONTRIBUTOR = "Cont";	// allow both Cont(ributor)
	private final String KEY_INSTTYPE = "InstType:";
	private final String KEY_INSTNAME = "InstName:";
	private final String KEY_IONIMETHOD = "IoniMethod:";
	private final String KEY_IONPOLARITY = "IonPol";	// allow both IonPol and IonPolarity
	private final String KEY_MSMS = "MSMS";		// allow both MSMS(Stage)
	private final String KEY_PREION = "PreIon:";
	private final String KEY_PRODION = "ProdIon:";
	private final String KEY_TRAPDRIVE = "TrapDrive:";
	private final String KEY_SKIM1 = "Skim1:";
	private final String KEY_FRAGAMPL = "FragAmpl:";
	private final String KEY_ISOLWIDTH = "IsolWidth:";
	private final String KEY_TARGETGAS = "TargetGas:";
	private final String KEY_TARGETGASPRES = "TargetGasPres";	// allow both TargetGasPres(sure)
	private final String KEY_REAGENTION = "ReagentIon:";
	private final String KEY_REAGENTGASPRES = "ReagentGasPres";	// allow both ReagentGasPres(sure)
	private final String KEY_COLENERGY = "ColEnergy:";
	private final String KEY_PEAKWIDTH = "PeakWidth:";
	private final String KEY_REFLECTOR = "Refl";	// allow both Refl(ector)
	private final String KEY_PSD = "PSD:";
	private final String KEY_CHARGEDECONVOLVED = "ChargeDecon";	// allow both ChargeDecon(voluted)
	private final String KEY_DATE = "Date:";
	private final String KEY_COLUMN = "Column:";
	private final String KEY_RETTIME = "RetTime:";
	private final String KEY_SSID = "SSID:";
	private final String KEY_ANALID = "AnalID:";
	private final String KEY_ANALNAME = "AnalName:";
	private final String KEY_MASSRANGE = "Mass";	// allow both Mass(Range)
	private final String KEY_NUMPEAKS = "Num";		// allo both Num(Peaks)

	
	public enum InstType {IT, TQ, Q, TOF, ICR, FTMS, ESI_TOF};
	public enum IoniMethod {EI, CI, APCI, ESI, nano_ESI, TS, MALDI, CAESIUM, APMALDI, APPI};
	
	
	/**
	 * 1.param library file
	 * 2.param output path
	 * 3.param prefix
	 *
	 * @param args the arguments
	 * @param prop 
	 */
	public static void main(String[] args, Properties prop) {
		// /home/mgerlich/workspace_new/Edda/test_lib/Metabo_Q_Lib_PubChemID.library /home/mgerlich/workspace_new/Edda/test_lib/ BD
		
		if(args == null || args.length == 0 || args.length != 3) {
			System.err.println("Missing arguments!");
			System.err.println("Usage: /path/to/library /output/path prefix");
			System.err.println("Note that the prefix are 2 or 3 upper case letters!");
			System.exit(-1);
		}
		
		String libfile = args[0];
		String outpath = args[1];
		String prefix = args[2];
		
		LibraryToMassBank lto = new LibraryToMassBank();
		lto.setOutputPath(outpath);
		lto.setPrefix(prefix);
		lto.prop = prop;
		
		File f = new File(libfile);
		try {
			lto.convertFile(f);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * Convert string to hex.
	 *
	 * @param str the str
	 * @return the string
	 */
	public String convertStringToHex(String str) {

		char[] chars = str.toCharArray();

		StringBuffer hex = new StringBuffer();
		for (int i = 0; i < chars.length; i++) {
			hex.append(Integer.toHexString((int) chars[i]));
		}

		return hex.toString();
	}

	/**
	 * Convert hex to string.
	 *
	 * @param hex the hex
	 * @return the string
	 */
	public String convertHexToString(String hex) {

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

	/**
	 * Parses the Bruker Library file.
	 *
	 * @param f the f
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public void parseBrukerLibraryFile(File f) throws IOException {
		if(!f.exists()) {
			System.err.println("File [" + f.getAbsolutePath() + "] not found!");
			return;
		}
		
		// generate moldata directory
		File dir = new File(outputPath, "moldata");
		dir.mkdir();
		this.molPath = dir.getAbsolutePath();
		
		// generate recdata directory
		dir = new File(outputPath, "recdata");
		dir.mkdir();
		this.recPath = dir.getAbsolutePath();
		
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = "";
		int counter = 1;
		String id = this.getPrefix() + formatCounter(counter);
		double emass = 0.0d;
		String name = "", cas = "", formula = "";
		String date = "2007.10.09";
		String year = "";
		String authors = null; //"Bruker Scientific"; should come from properties instead of defaults here
		//String contributor = "";
		String smiles = "", inchi = "", inchikey = "";
		String instrument = null; //"micrOTOF-Q"; should come from properties instead of defaults here
		String instrument_type = null; //"ESI-TOF";
		String compound_class = "CH$COMPOUND_CLASS: unknown";
		String ev = "";
		String ion_method = "ESI";
		String ion_mode = "pos";
		String ion = "";
		String ms = "MS/MS";
		String precursor = "";
		int numPeaks = 0;
		StringBuffer peaks = new StringBuffer();
		String link = "";
		List<String> synonyme = new ArrayList<String>();
		List<String> names = new ArrayList<String>();
		
		// list writer for matching compoundnames to accession numbers
		FileWriter listWriter = new FileWriter(listMol);
		
		while((line = br.readLine()) != null) {
			if(line.startsWith("Name:"))
				name = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("CAS:"))
				cas = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("Formula:")) {
				formula = line.substring(line.indexOf(":") + 1).trim();
				emass = computeEmass(formula);
			}
			if(line.startsWith("Synonym:")) {
				String syn = line.substring(line.indexOf(":") + 1).trim();
				if(syn.matches("C[0-9]{5}")) {	// found KEGG
					synonyme.add("KEGG " + syn);
				}
				else if(syn.contains("CID:")) {	// found Pubchem
					syn = syn.replaceAll(" ", "");
					synonyme.add("PUBCHEM " + syn);
				}
				else if(syn.contains(",")) {	// found names
					String[] split = syn.split(",");
					for (int i = 0; i < split.length; i++) {
						names.add(split[i]);
					}
				}
				else {							// found names
					names.add(syn);
				}
			}
			if(line.startsWith("Comment:")) {
				String temp = line.substring(line.indexOf(":") + 1).trim();
				if(!synonyme.contains("KEGG " + temp) && temp.matches("C[0-9]{5}"))	// KEGG ID in comment and not in Synonym
					synonyme.add("KEGG " + temp);
			}
			if(line.startsWith("Structure:")) {
				String mol = convertHexToString(line.substring(line.indexOf(":") + 1).trim());
				File mdl = new File(molPath, id + ".mol");
				System.out.println("write mol -> " + mdl);
				FileWriter fw = new FileWriter(mdl);
				fw.write(mol);
				fw.flush();
				fw.close();

				FileInputStream ins = new FileInputStream(mdl);
				MDLV2000Reader reader = new MDLV2000Reader(ins);
				IMolecule molecule=null;
				try {
					molecule = (IMolecule) reader.read(new ChemFile());
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				
				// Calculate SMILES 
				String molsmiles = "";
				SmilesGenerator generator = new SmilesGenerator();
				generator.setUseAromaticityFlag(true);
				molsmiles = generator.createSMILES(molecule);
				
				// get molecular Formula
				String molformula = MolecularFormulaManipulator.getString(MolecularFormulaManipulator.getMolecularFormula(molecule));
				
				// Debug: output FOrmula + SMILES
				System.out.println("Formula:" + molformula + " and SMILES: " + molsmiles);
								
				// 
				
				listWriter.append(name + "\t" + id + ".mol" + "\n");
			}
			if(line.startsWith("InChI:"))
				inchi = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("InChIKey:"))
				inchikey = line.substring(line.indexOf(":") + 1).trim();

			if(line.startsWith("Contributor:"))
				authors = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("InstType:"))
				instrument_type = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("InstName:")) {
				instrument = "Bruker " + line.substring(line.indexOf(":") + 1).trim();
			}
			if(line.startsWith("IoniMethod:"))
				ion_method = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("IonPol")) {
				ion_mode = line.substring(line.indexOf(":") + 1).trim();
				if(ion_mode.equalsIgnoreCase("pos")) {
					ion_mode = "ION_MODE POSITIVE";
					ion = "[M+H]+";
				}
				else if(ion_mode.equalsIgnoreCase("neg")) {
					ion_mode = "ION_MODE NEGATIVE";
					ion = "[M-H]-";
				}					
			}
			if(line.startsWith("MSMS:"))
				ms = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("PreIon:"))
				precursor = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith("ColEnergy:")) {
				ev = line.substring(line.indexOf(":") + 1).trim();
				if(ev.contains("."))
					ev = ev.substring(0, ev.indexOf("."));
			}
			if(line.startsWith("Date:")) {
				Date today = new Date();
				SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("yyyy.MM.dd");
				
				date = DATE_FORMAT.format(today) + " (Created " + line.substring(line.indexOf(":") + 1, line.indexOf(":") + 12).trim() + ")";
				DATE_FORMAT = new SimpleDateFormat("yyyy");
				year = DATE_FORMAT.format(today);
			}
			
			if(line.startsWith("Num Peaks:")) {
				numPeaks = Integer.parseInt(line.substring(line.indexOf(":") + 1).trim());
				
				// read in peaklist
				line = br.readLine();
				String[] split = line.split(" ");
				System.out.println("numPeaks -> " + numPeaks + "\t#split -> " + split.length);
				
				if((split.length == numPeaks*2) && (split.length % 2 == 0)) {		// all peaks in one line
					for (int i = 0; i < split.length; i=i+2) {
						peaks.append("  " + split[i] + " " + split[i+1] + "  " + split[i+1] + "\n");
					}
					peaks.append("//");
				}
				else if(split.length < numPeaks*2){		// peaks are spread over multiple lines
					// peaklist geht ueber mehrere zeilen
					// TODO
					br.mark(2000);
					String temp = "";
					do {
						temp = br.readLine();
						line += temp;
					}
					while(!temp.isEmpty() && !temp.equals("\n"));
					br.reset();
					
					split = line.split(" ");
					System.out.println(line);
					System.out.println(split.length);
					for (int i = 0; i < split.length; i=i+2) {
						peaks.append("  " + split[i] + " " + split[i+1] + "  " + split[i+1] + "\n");
					}
					peaks.append("//");
				}
				
				
				//if(split.length == 2*numPeaks) {
//				if(split.length % 2 == 0) {
//					if((split.length / 2) < numPeaks) {
//						System.err.println("reduced numpeaks from " + numPeaks + " to " + split.length / 2);
//						numPeaks = split.length / 2;
//					}
//					for (int i = 0; i < split.length; i=i+2) {
//						peaks.append("  " + split[i] + " " + split[i+1] + "  " + split[i+1] + "\n");
//					}
//					peaks.append("//");
//				}
			}
				
				
			/**
			 * write out current record
			 */
			if(line.isEmpty() || line.equals("\n")) {
				System.out.println("ID -> " + id);
				
				if(ev.isEmpty())
					ev = "10";
				if(instrument.isEmpty())
					instrument = "micrOTOF-Q";
				if(instrument_type.isEmpty())
					instrument_type = "ESI-TOF";
				if(ion_mode.isEmpty())
					ion_mode = "ION_MODE POSITIVE";
				if(ion_method.isEmpty())
					ion_method = "ESI";
				if(authors.isEmpty())
					authors = "Bruker";
				
				File temp = new File(recPath, id + ".txt");
				FileWriter fw = new FileWriter(temp);
				
				// write information
				fw.write("ACCESSION: " + id);
				fw.write("\n");
				fw.write("RECORD_TITLE: " + name + "; MS/MS; QqTOF; CE:" + ev + " V; " + ion);
				fw.write("\n");
				fw.write("DATE: " + date);
				fw.write("\n");
//				fw.write("AUTHORS: " + contributor);
				fw.write("AUTHORS: " + prop.getProperty("authors"));//, contributor));

				fw.write("\n");
				fw.write("COPYRIGHT: Copyright(C) " + year );
				fw.write("\n");
				fw.write("CH$NAME: " + name);
				fw.write("\n");
				if(names != null && names.size() > 0) {
					for (String s : names) {
						fw.write("CH$NAME: " + s.trim());
						fw.write("\n");
					}
				}
				fw.write("CH$COMPOUND_CLASS: N/A");
				fw.write("\n");
				fw.write("CH$FORMULA: " + formula);
				fw.write("\n");
				fw.write("CH$EXACT_MASS: " + String.format("%.4f", emass));
				fw.write("\n");
				fw.write("CH$SMILES: not available");
				fw.write("\n");
				fw.write("CH$IUPAC: " + inchi );
				fw.write("\n");
				if(synonyme != null && synonyme.size() > 0) {
					for (String s : synonyme) {
						fw.write("CH$LINK: " + s);
						fw.write("\n");
					}
				}
				if(!cas.isEmpty()) {
					fw.write("CH$LINK: CAS " + cas);
					fw.write("\n");
				}
					
				fw.write("AC$INSTRUMENT: " + instrument);
				fw.write("\n");
				fw.write("AC$INSTRUMENT_TYPE: " + instrument_type);
				fw.write("\n");
				fw.write("AC$MASS_SPECTROMETRY: " + ion_mode);
				fw.write("\n");
				fw.write("AC$MASS_SPECTROMETRY: COLLISION_ENERGY " + ev + " eV");
				fw.write("\n");
				
				if(!precursor.isEmpty()) {
					fw.write("MS$FOCUSED_ION: PRECURSOR_M/Z " + precursor);
					fw.write("\n");
					fw.write("MS$FOCUSED_ION: PRECURSOR_TYPE " + ion);
					fw.write("\n");
				}
				
				fw.write("PK$NUM_PEAK: " + numPeaks);
				fw.write("\n");
				fw.write("PK$PEAK: m/z int. rel.int.");
				fw.write("\n");
				fw.write(peaks.toString());
				
				// close FileWriter
				fw.flush();
				fw.close();
				
				// raise counter and create new ID
				counter++;
				id = getPrefix() + formatCounter(counter);
				
				// clear current variables
				peaks = new StringBuffer();
				link = "";
				name = "";
				cas = "";
				emass = 0;
				names = new ArrayList<String>();
				synonyme = new ArrayList<String>();
				authors = "";
				instrument = "";
				instrument_type = "";
				ion_method = "";
				ion_mode = "";
				ms = "";
				ev = "";
				date = "";
				year = "";
				numPeaks = 0;
			}
			
		}
		listWriter.flush();
		listWriter.close();
		br.close();
	}
	
	/**
	 * Compute emass.
	 *
	 * @param formula the formula
	 * @return the double
	 */
	public double computeEmass(String formula) {
		IMolecularFormula mf = MolecularFormulaManipulator.getMolecularFormula(formula, NoNotificationChemObjectBuilder.getInstance());
		return MolecularFormulaManipulator.getTotalExactMass(mf);	//MolecularFormulaManipulator.getNaturalExactMass(mf);
	}
	
	/**
	 * Format counter.
	 *
	 * @param i the i
	 * @return the string
	 */
	public String formatCounter(int i) {
		DecimalFormat format = new DecimalFormat("000000");
		
		if(this.prefix.length() == 2)
			format = new DecimalFormat("000000");
		if(this.prefix.length() == 3)
				format = new DecimalFormat("00000");
		
		return format.format(i);
	}
	

	/**
	 * Convert file.
	 *
	 * @param f the f
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public void convertFile(File f) throws IOException {
		if(!f.exists()) {
			System.err.println("File [" + f.getAbsolutePath() + "] not found!");
			return;
		}
		
		// generate proper date format
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd");
		Date d = new Date();
		String date = sdf.format(d);
		
		// generate moldata directory
		File dir = new File(outputPath, "moldata");
		boolean io = dir.mkdirs();
		if(dir.exists())
			io = true;
		if(!io) {
			System.err.println("Error creating directory [" + dir.getAbsolutePath() + "]!");
			System.err.println("Please check path and access rights - exiting!");
			System.exit(-1);
		}
		this.molPath = dir.getAbsolutePath();
		
		// generate new list file containing assignment of compoundnames to molfiles
		this.listMol = new File(dir, "list.tsv");
		
		// generate recdata directory
		dir = new File(outputPath, "recdata");
		io = dir.mkdir();
		if(dir.exists())
			io = true;
		if(!io) {
			System.err.println("Error creating directory [" + dir.getAbsolutePath() + "]!");
			System.err.println("Please check path and access rights - exiting!");
			System.exit(-1);
		}
		this.recPath = dir.getAbsolutePath();
		
		// read in library file
			//BufferedReader br = new BufferedReader(new FileReader(f));
		// use correct encoding
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(f.getAbsoluteFile()), "ISO-8859-1"));
		String line = "";
		
		int counter = 1;
		int numPeaks = 0;
		String id = this.getPrefix() + formatCounter(counter);
		double emass = 0.0d;
		String name = "", cas = "", pubchem="", nist = "", un = "", formula = "", mw = ""; 
		String year = "", link = "", comment = "";
		String smiles = "", inchi = "", inchikey = "";  
		String preion = "", prodion = "", trapdrive = "", skim = "", fragampl = "", isolwidth = "", targetgas = "", targetgaspres = "",
			reagention = "", reagentgaspres = "", peakwidth = "", reflector = "", psd = "", chargedecon = "", column = "", rettime = "",
			ssid = "", analid = "", analname = "", massrange = "";
		String authors = prop.getProperty("authors");

		String instrument = "micrOTOF-Q";
		String instrument_type = "ESI-TOF";
		String compound_class_unknown = "CH$COMPOUND_CLASS: unknown";
		String compound_class = "";
		String ev = "", ion = "", precursor = "";
		String ion_method = "ESI";
		String ion_mode = "unknown";
		String ms = "";

		StringBuffer peaks = new StringBuffer();
		List<String> synonyms = new ArrayList<String>();
		List<String> names = new ArrayList<String>();
		
		// list writer for matching compoundnames to accession numbers
		FileWriter listWriter = new FileWriter(listMol);
		Map<String, String> uniqueCompounds = new HashMap<String, String>();
		
		while((line = br.readLine()) != null) {
			if(line.startsWith(KEY_NAME)) {
				name = line.substring(line.indexOf(":") + 1).trim();
			}
			if(line.startsWith(KEY_CAS))
				cas = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_PUBCHEM))
				pubchem = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_NIST))
				nist = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_UN))
				un = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_MW))
				mw = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_INCHI))
				inchi = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_INCHIKEY))
				inchikey = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_FORMULA)) {
				formula = line.substring(line.indexOf(":") + 1).trim();
				emass = computeEmass(formula);
			}
			if(line.startsWith(KEY_SYNONYM)) {
				String syn = line.substring(line.indexOf(":") + 1).trim();
				if(syn.matches("C[0-9]{5}")) {	// found KEGG
					synonyms.add("KEGG " + syn);
				}
				else if(syn.contains("CID:")) {	// found Pubchem
					syn = syn.replaceAll(" ", "");
					synonyms.add("PUBCHEM " + syn);
				}
				// TODO: split synonyms correctly
//				else if(syn.contains(" ,")) {	// found names
//					String[] split = syn.split(" ,");
//					for (int i = 0; i < split.length; i++) {
//						names.add(split[i].trim());
//					}
//				}
				else {							// found names
					names.add(syn);
				}
			}
			if(line.startsWith(KEY_COMMENT)) {
				comment = line.substring(line.indexOf(":") + 1).trim();
				if(comment.contains("Compound class") || comment.contains("compound class")) {
					String temp = line.substring(line.indexOf(":") + 1).trim();
					temp = temp.replace("Compound class", "");
					temp = temp.replace("compound class", "");
					if(temp.contains(":"))
						temp = temp.replace(":", "");
					compound_class = temp.trim(); 	//line.substring(line.indexOf(":") + 1).trim();
				}
//				String temp = line.substring(line.indexOf(":") + 1).trim();
//				if(!synonyms.contains("KEGG " + temp) && temp.matches("C[0-9]{5}"))	// KEGG ID in comment and not in Synonym
//					synonyms.add("KEGG " + temp);
			}
			if(line.startsWith(KEY_STRUCTURE)) {
				String mol = convertHexToString(line.substring(line.indexOf(":") + 1).trim());
				File mdl = new File(molPath, id + ".mol");
				System.out.println("write mol -> " + mdl);
				FileWriter fw = new FileWriter(mdl);
				fw.write(mol);
				fw.flush();
				fw.close();
				
				InputStream is = new ByteArrayInputStream(mol.getBytes());
				MDLV2000Reader reader = new MDLV2000Reader(is);
				IChemFile chemFile = new ChemFile();
				IAtomContainer container = null;
				try {
					chemFile = (IChemFile) reader.read(chemFile);
					container = ChemFileManipulator.getAllAtomContainers(chemFile).get(0);
					SmilesGenerator sg = new SmilesGenerator();
					sg.setUseAromaticityFlag(true);
					smiles = sg.createSMILES(container);

					// get molecular Formula
					formula = MolecularFormulaManipulator.getString(MolecularFormulaManipulator.getMolecularFormula(container));

					emass = AtomContainerManipulator.getTotalExactMass(container);
					
//					// Generate factory - throws CDKException if native code does not load
//					InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
//					// Get InChIGenerator
//					InChIGenerator gen = factory.getInChIGenerator(container);
//
//					INCHI_RET ret = gen.getReturnStatus();
//					if (ret == INCHI_RET.WARNING) {
//						// InChI generated, but with warning message
//						System.out.println("InChI warning: " + gen.getMessage());
//					} else if (ret != INCHI_RET.OKAY) {
//						inchi = "";
//						// InChI generation failed
//						throw new CDKException("InChI failed: " + ret.toString() + " [" + gen.getMessage() + "]");
//					}
//
//					inchi = gen.getInchi();
												
				} catch (CDKException e) {
					System.err.println("Error loading molfile and generating SMILES!");
					smiles = "";
				}
				
				if(!uniqueCompounds.containsKey(name))
					uniqueCompounds.put(name, id + ".mol");
				//listWriter.append(name + "\t" + id + ".mol" + "\n");
			}
			if(line.startsWith(KEY_CONTRIBUTOR)) {
				authors = line.substring(line.indexOf(":") + 1).trim();
//				System.out.println(contributor);
//				byte[] arr = contributor.getBytes();
//				String s = new String(arr, "UTF8");
//				System.out.println(s);
			}
			
			if(line.startsWith(KEY_INSTTYPE)) {
				instrument_type = line.substring(line.indexOf(":") + 1).trim();
				// override from settings file
				if (prop.containsKey("instrument_type")) {
					instrument_type = prop.getProperty("instrument_type");
				}
			}
			if(line.startsWith(KEY_INSTNAME))
				instrument = "Bruker " + line.substring(line.indexOf(":") + 1).trim();			
			if(line.startsWith(KEY_IONIMETHOD))
				ion_method = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_IONPOLARITY)) {
				ion_mode = line.substring(line.indexOf(":") + 1).trim();
				if(ion_mode.equalsIgnoreCase("pos")) {
					ion_mode = "ION_MODE POSITIVE";
					ion = "[M+H]+";
				}
				else if(ion_mode.equalsIgnoreCase("neg")) {
					ion_mode = "ION_MODE NEGATIVE";
					ion = "[M-H]-";
				}
					
			}
			if(line.startsWith(KEY_MSMS)) {
				ms = line.substring(line.indexOf(":") + 1).trim();
				if(ms.equals("1"))
					ms = "MS";
				else if(ms.equals("2"))
					ms = "MS2";			// new
					//ms = "MS/MS";		// old
				else if(ms.equals("3"))
					ms = "MS3";
				else if(ms.equals("4"))
					ms = "MS4";
				else ms = "MS" + ms;
			}
			if(line.startsWith(KEY_PREION))
				precursor = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_PRODION))
				prodion = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_TRAPDRIVE))
				trapdrive = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_SKIM1))
				skim = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_FRAGAMPL))
				fragampl = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_ISOLWIDTH))
				isolwidth = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_TARGETGAS))
				targetgas = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_TARGETGASPRES))
				targetgaspres = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_REAGENTION))
				reagention = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_REAGENTGASPRES))
				reagentgaspres = line.substring(line.indexOf(":") + 1).trim();
			
			if(line.startsWith(KEY_COLENERGY)) {
				ev = line.substring(line.indexOf(":") + 1).trim();
				if(ev.contains("."))
					ev = ev.substring(0, ev.indexOf("."));
			}
			
			if(line.startsWith(KEY_PEAKWIDTH))
				peakwidth = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_REFLECTOR))
				reflector = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_PSD))
				psd = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_CHARGEDECONVOLVED))
				chargedecon = line.substring(line.indexOf(":") + 1).trim();
			
			if(line.startsWith(KEY_DATE)) {

				Date today = new Date();
				SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("yyyy.MM.dd");				
				
				date = line.substring(line.indexOf(":") + 1, line.indexOf(":") + 12).trim();
				date = DATE_FORMAT.format(today) + " (Created " + date.replaceAll("-", ".")  + ")";

				DATE_FORMAT = new SimpleDateFormat("yyyy");
				year = DATE_FORMAT.format(today);
			}
			
			if(line.startsWith(KEY_COLUMN))
				column = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_RETTIME))
				rettime = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_SSID))
				ssid = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_ANALID))
				analid = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_ANALNAME))
				analname = line.substring(line.indexOf(":") + 1).trim();
			if(line.startsWith(KEY_MASSRANGE)) {
				massrange = line.substring(line.indexOf(":") + 1).trim();
				massrange = massrange.replace(",", "-");
			}
			
			if(line.startsWith(KEY_NUMPEAKS)) {
				numPeaks = Integer.parseInt(line.substring(line.indexOf(":") + 1).trim());
				
				// read in peaklist
				line = br.readLine();
				String[] split = line.split(" ");
				System.out.println("numPeaks -> " + numPeaks + "\t#split -> " + split.length);
				
				if((split.length == numPeaks*2) && (split.length % 2 == 0)) {		// all peaks in one line
					for (int i = 0; i < split.length; i=i+2) {
						peaks.append("  " + split[i] + " " + split[i+1] + " " + split[i+1] + "\n");
					}
					peaks.append("//");
				}
				else if(split.length < numPeaks*2){		// peaks are spread over multiple lines
					br.mark(20000);
					String temp = "";
					do {
						temp = br.readLine();
						line += " " + temp;
					}
					while(temp != null && !temp.isEmpty() && !temp.equals("\n"));
					br.reset();
					
					split = line.split(" ");
					System.out.println(line);
					System.out.println(split.length);
					for (int i = 0; i < split.length; i=i+2) {
						peaks.append("  " + split[i] + " " + split[i+1] + " " + split[i+1] + "\n");
					}
					peaks.append("//");
				}
			}
				
				
			/**
			 * write out current record
			 */
			if(line.isEmpty() || line.equals("\n")) {
				System.out.println("ID -> " + id);
				
				if(ev.isEmpty())
					ev = "10";
				if(instrument.isEmpty())
					instrument = "micrOTOF-Q";
				if(instrument_type.isEmpty())
					instrument_type = "ESI-TOF";
				if(ion_mode.isEmpty())
					ion_mode = "ION_MODE POSITIVE";
				if(ion_method.isEmpty())
					ion_method = "ESI";
				if(authors.isEmpty())
					authors = "Bruker Library";
				
				File temp = new File(recPath, id + ".txt");
				FileWriter fw = new FileWriter(temp);
				
				// write information
				fw.write("ACCESSION: " + id);
				fw.write("\n");
				if(ion_method.contains("EI"))
					fw.write("RECORD_TITLE: " + name + "; " + instrument_type + "RT:" + rettime + " sec");
				else if(ion_method.contains("ESI"))
					fw.write("RECORD_TITLE: " + name + "; " + instrument_type + "; " + ms + "; CE:" + ev + " eV; " + ion);
				else fw.write("RECORD_TITLE: " + name + "; " + instrument_type + "; " + ms + "; CE:" + ev + " eV; " + ion);
				fw.write("\n");
				fw.write("DATE: " + date);
				fw.write("\n");
				fw.write("AUTHORS: " + authors);
				fw.write("\n");
				fw.write("COPYRIGHT: " + prop.getProperty("copyright"));
				fw.write("\n");
				fw.write("LICENSE: " + prop.getProperty("license"));

				fw.write("\n");
				if(!comment.isEmpty()) {
					fw.write("COMMENT: " + comment);
					fw.write("\n");
				}
				fw.write("CH$NAME: " + name);
				fw.write("\n");
				if(names != null && names.size() > 0) {
					for (String s : names) {
						fw.write("CH$NAME: " + s);
						fw.write("\n");
					}
				}
				if(!compound_class.isEmpty()) {
					fw.write("CH$COMPOUND_CLASS: " + compound_class);
					//fw.write("\n");
				}else fw.write(compound_class_unknown);
				fw.write("\n");
				fw.write("CH$FORMULA: " + formula);
				fw.write("\n");
				fw.write("CH$EXACT_MASS: " + String.format("%.4f", emass));
				fw.write("\n");
				if(smiles.isEmpty())
					fw.write("CH$SMILES: not available");
				else 
					fw.write("CH$SMILES: " + smiles);
				fw.write("\n");

				if(inchi.isEmpty())
					fw.write("CH$IUPAC: not available");
				else 
					fw.write("CH$IUPAC: " + inchi);
				fw.write("\n");
				
				if(!inchikey.isEmpty()) {
					fw.write("CH$LINK: INCHIKEY " + inchikey);
					fw.write("\n");
				}
				
				if(synonyms != null && synonyms.size() > 0) {
					for (String s : synonyms) {
						fw.write("CH$LINK: " + s);
						fw.write("\n");
					}
				}
				if(!cas.isEmpty()) {
					fw.write("CH$LINK: CAS " + cas);
					fw.write("\n");
				}
				if(!pubchem.isEmpty()) {
					fw.write("CH$LINK: PUBCHEM CID:" + pubchem);
					fw.write("\n");
				}
				if(!nist.isEmpty()) {
					fw.write("CH$LINK: NIST " + nist);
					fw.write("\n");
				}	
				if(!un.isEmpty()) {
					fw.write("CH$LINK: UN " + un);
					fw.write("\n");
				}
				fw.write("AC$INSTRUMENT: " + instrument);
				fw.write("\n");
				fw.write("AC$INSTRUMENT_TYPE: " + instrument_type);
				fw.write("\n");
				fw.write("AC$MASS_SPECTROMETRY: " + ion_mode);
				fw.write("\n");
				fw.write("AC$MASS_SPECTROMETRY: COLLISION_ENERGY " + ev + " eV");
				fw.write("\n");
				if(!ms.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: MS_TYPE " + ms);	// old
					//fw.write("AC$MASS_SPECTROMETRY: MS_TYPE " + ms);		// new
					fw.write("\n");
				}
				if(!ion_method.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: IONIZATION " + ion_method);	// old
					//fw.write("AC$MASS_SPECTROMETRY: ION_" + ion_mode);	// new
					fw.write("\n");
				}
				if(!trapdrive.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_TRAPDRIVE + " " + trapdrive);
					fw.write("\n");
				}
				if(!skim.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_SKIM1 + " " + skim);
					fw.write("\n");
				}
				if(!fragampl.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_FRAGAMPL + " " + fragampl);
					fw.write("\n");
				}
				if(!isolwidth.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_ISOLWIDTH + " " + isolwidth);
					fw.write("\n");
				}
				if(!targetgas.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_TARGETGAS + " " + targetgas);
					fw.write("\n");
				}
				if(!targetgaspres.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_TARGETGASPRES + " " + targetgaspres);
					fw.write("\n");
				}
				if(!reagention.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_REAGENTION + " " + reagention);
					fw.write("\n");
				}
				if(!reagentgaspres.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_REAGENTGASPRES + " " + reagentgaspres);
					fw.write("\n");
				}
				if(!peakwidth.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_PEAKWIDTH + " " + peakwidth);
					fw.write("\n");
				}
				if(!reflector.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_REFLECTOR + " " + reflector);
					fw.write("\n");
				}
				if(!psd.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_PSD + " " + psd);
					fw.write("\n");
				}
				if(!chargedecon.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_CHARGEDECONVOLVED + " " + chargedecon);
					fw.write("\n");
				}
				if(!column.isEmpty()) {
					fw.write("AC$CHROMATOGRAPHY: " + KEY_COLUMN + " " + column);
					fw.write("\n");
				}
				if(!rettime.isEmpty()) {
					fw.write("AC$CHROMATOGRAPHY: RETENTION_TIME " + rettime + " sec");
					fw.write("\n");
				}

				/** skip bruker specific data (filesystem) for public access */
//				if(!ssid.isEmpty()) {
//					fw.write("AC$ANALYTICAL_CONDITION: " + KEY_SSID + " " + ssid);
//					fw.write("\n");
//				}
//				if(!analid.isEmpty()) {
//					fw.write("AC$ANALYTICAL_CONDITION: " + KEY_ANALID + " " + analid);
//					fw.write("\n");
//				}
//				if(!analname.isEmpty()) {
//					fw.write("AC$ANALYTICAL_CONDITION: " + KEY_ANALNAME + " " + analname);
//					fw.write("\n");
//				}
				/** */

				if(!massrange.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: " + KEY_MASSRANGE + " " + massrange);
					fw.write("\n");
				}
				
				if(!precursor.isEmpty()) {
					fw.write("MS$FOCUSED_ION: PRECURSOR_M/Z " + precursor);
					fw.write("\n");
					fw.write("MS$FOCUSED_ION: PRECURSOR_TYPE " + ion);
					fw.write("\n");
				}
				if(!prodion.isEmpty()) {
					fw.write("MS$FOCUSED_ION: BASE_PEAK " + prodion);
					fw.write("\n");
				}

				fw.write("MS$DATA_PROCESSING: CONVERT from Bruker Library Editor (https://github.com/MassBank/MassBank2NIST)\n");

				fw.write("PK$NUM_PEAK: " + numPeaks);
				fw.write("\n");
				fw.write("PK$PEAK: m/z int. rel.int.");
				fw.write("\n");
				fw.write(peaks.toString());
				
				// close FileWriter
				fw.flush();
				fw.close();
				
				// raise counter and create new ID
				counter++;
				id = getPrefix() + formatCounter(counter);
				
				// clear current variables
				peaks = new StringBuffer();
				link = "";
				name = "";
				cas = "";
				emass = 0;
				names = new ArrayList<String>();
				synonyms = new ArrayList<String>();
				instrument = "";
				instrument_type = "";
				ion_method = "";
				ion_mode = "";
				ms = "";
				ev = "";
				date = "";
				year = "";
				numPeaks = 0;
				smiles = "";
				inchi = "";
			}
			
		}
		
		Set<String> keys = uniqueCompounds.keySet();
		for (Iterator<String> it = keys.iterator(); it.hasNext();) {
			String compound = it.next();
			String mol = uniqueCompounds.get(compound);
			listWriter.append(compound + "\t" + mol + "\n");
		}
		listWriter.flush();
		listWriter.close();
		br.close();
	}
	
	/**
	 * Sets the output path.
	 *
	 * @param outputPath the new output path
	 */
	public void setOutputPath(String outputPath) {
		this.outputPath = outputPath;
	}

	/**
	 * Gets the output path.
	 *
	 * @return the output path
	 */
	public String getOutputPath() {
		return outputPath;
	}

	/**
	 * Sets the prefix.
	 *
	 * @param prefix the new prefix
	 */
	public void setPrefix(String prefix) {
		this.prefix = prefix;
	}

	/**
	 * Gets the prefix.
	 *
	 * @return the prefix
	 */
	public String getPrefix() {
		return prefix;
	}
}