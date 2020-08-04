package de.ipbhalle.converter;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
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
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import net.sf.jniinchi.INCHI_RET;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import edu.ucdavis.fiehnlab.spectra.hash.core.Spectrum;
import edu.ucdavis.fiehnlab.spectra.hash.core.Splash;
import edu.ucdavis.fiehnlab.spectra.hash.core.SplashFactory;
import edu.ucdavis.fiehnlab.spectra.hash.core.types.Ion;
import edu.ucdavis.fiehnlab.spectra.hash.core.types.SpectraType;
import edu.ucdavis.fiehnlab.spectra.hash.core.types.SpectrumImpl;

// TODO: Auto-generated Javadoc
/**
 * The Class LibraryToMassBank.
 */
public class LibraryToMassBank {

	/** The prefix. */
	private String prefix = "XX";
	
	/** First Accession ID. */
	private int startswithid = 1;

	public int getStartswithid() {
		return startswithid;
	}

	public void setStartswithid(int startswithid) {
		this.startswithid = startswithid;
	}

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
		if (prop.containsKey("startswithid")) {
			lto.setStartswithid(Integer.parseInt(prop.getProperty("startswithid")));
		} else {
			lto.setStartswithid(1);
		}
		lto.prop = prop;
		
		File f = new File(libfile);
		try {
			try {
				lto.convertFile(f);
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.err.println("Done.");
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
	 * Compute emass.
	 *
	 * @param formula the formula
	 * @return the double
	 */
	public double computeEmass(String formula) {
		IMolecularFormula mf = MolecularFormulaManipulator.getMolecularFormula(formula, SilentChemObjectBuilder.getInstance());
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
	 * @throws CDKException 
	 */
	public void convertFile(File f) throws IOException, CDKException {
		if(!f.exists()) {
			System.err.println("File [" + f.getAbsolutePath() + "] not found!");
			return;
		}
		
		// generate proper date format
		String date = new SimpleDateFormat("yyyy.MM.dd").format(new Date());
		
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
		
		int counter = getStartswithid();
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
		String compound_class_unknown = "CH$COMPOUND_CLASS: Natural Product; N/A";
		String compound_class = "";
		String ev = "";
		String ion = "";
		String precursor = "";
		String ion_method = "ESI";
		String ion_mode = "unknown";
		String ms = "MS1";

		String lc_solvent_a = "";
		String lc_solvent_b = "";
		String lc_gradient = "";
		String lc_flow_rate = "";
		
		StringBuffer peaks = new StringBuffer();
		ArrayList<Ion> ionList = new ArrayList<Ion>();
						
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
				// Sanitize molfile
				// remove trailling empty lines: https://stackoverflow.com/a/37648209
				mol = mol.replaceAll("([\\n\\r]+\\s*)*$", "");
				
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
					System.err.println("Formula is " + formula);
					Isotopes.getInstance().configureAtoms(container);
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
												
				} catch (IllegalArgumentException e) {
					System.err.println("molfile contains PseudoAtom: " + e);
					// Generate factory - throws CDKException if native code does not load
					InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();;
					                                    
					// Get InChIToStructure
					InChIToStructure intostruct = factory.getInChIToStructure(
					inchi, DefaultChemObjectBuilder.getInstance()
					);

					INCHI_RET ret = intostruct.getReturnStatus();
					if (ret == INCHI_RET.WARNING) {
					// Structure generated, but with warning message
					System.out.println("InChI warning: " + intostruct.getMessage());
					} else if (ret != INCHI_RET.OKAY) {
					// Structure generation failed
					throw new CDKException("Structure generation failed failed: " + ret.toString()
					+ " [" + intostruct.getMessage() + "]");
					}

					container = intostruct.getAtomContainer();
					Isotopes.getInstance().configureAtoms(container);
					emass = AtomContainerManipulator.getTotalExactMass(container);					
					
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

			// obtain from settings file
			if (prop.containsKey("lc_solvent_a")) {
				lc_solvent_a = prop.getProperty("lc_solvent_a");
			}
			if (prop.containsKey("lc_solvent_b")) {
				lc_solvent_b = prop.getProperty("lc_solvent_b");
			}
			if (prop.containsKey("lc_gradient")) {
				lc_gradient = prop.getProperty("lc_gradient");
			}
			if (prop.containsKey("lc_flow_rate")) {
				lc_flow_rate = prop.getProperty("lc_flow_rate");
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
				line = line.trim();
				String[] split = line.split(" ");
				System.out.println("numPeaks -> " + numPeaks + "\t#split -> " + split.length);
				
				if(!((split.length == numPeaks*2) && (split.length % 2 == 0))) {		// peaks are spread over multiple lines
					br.mark(200000);
					String temp = "";
					do {
						temp = br.readLine();
						temp = temp.trim();
						line += " " + temp;
					}
					while(temp != null && !temp.isEmpty() && !temp.equals("\n"));
					br.reset();
				}
					
				split = line.split(" ");
					//System.out.println(line);
					//System.out.println(split.length);
				Double maxIntensity = new Double(0);
				for (int i = 0; i < split.length; i=i+2) {
					maxIntensity = Double.parseDouble(split[i+1]) > maxIntensity ? Double.parseDouble(split[i+1]) : maxIntensity;
				}
				
				for (int i = 0; i < split.length; i=i+2) {
					Integer relint = (int) Math.round( Double.parseDouble(split[i+1]) * 999 / maxIntensity);
					peaks.append("  " + split[i] + " " + split[i+1] + " " + relint + "\n");
					ionList.add(new Ion(Double.parseDouble(split[i]), Double.parseDouble(split[i+1])));
				}
				peaks.append("//\n");
			}
				
				
			/**
			 * write out current record
			 */
			if(line.isEmpty() || line.equals("\n")) {
				System.out.println("ID -> " + id);
				
				//if(ev.isEmpty()) // If collision energy is empty, then keep it empty.
				//	ev = "10";
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
				fw.write("ACCESSION: " + id + "\n");
				if(ion_method.contains("EI"))
					fw.write("RECORD_TITLE: " + name + "; " + instrument_type + "RT:" + rettime + " sec");
				else if(ion_method.contains("ESI"))
					fw.write("RECORD_TITLE: " + name + "; " + instrument_type + "; " + ms 
							+ ( ev.isEmpty() ? "" : ( "; CE:" + ev + " eV") + "; " + ion ) );
				else fw.write("RECORD_TITLE: " + name + "; " + instrument_type + "; " + ms + "; CE:" + ev + " eV; " + ion);
				fw.write("\n");
				fw.write("DATE: " + date);
				fw.write("\n");
				fw.write("AUTHORS: " + authors);
				fw.write("\n");
				fw.write("LICENSE: " + prop.getProperty("license"));
				fw.write("\n");
				fw.write("COPYRIGHT: " + prop.getProperty("copyright"));
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
				fw.write("CH$EXACT_MASS: " + String.format(Locale.US, "%.4f", emass));
				fw.write("\n");
				if(smiles.isEmpty())
					fw.write("CH$SMILES: N/A");
				else 
					fw.write("CH$SMILES: " + smiles);
				fw.write("\n");

				if(inchi.isEmpty())
					fw.write("CH$IUPAC: N/A");
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
				fw.write("AC$INSTRUMENT: " + instrument + "\n");
				fw.write("AC$INSTRUMENT_TYPE: " + instrument_type + "\n");
				fw.write("AC$MASS_SPECTROMETRY: MS_TYPE " + ms + "\n");
				fw.write("AC$MASS_SPECTROMETRY: " + ion_mode + "\n");
				if(!ev.isEmpty()) {
					fw.write("AC$MASS_SPECTROMETRY: COLLISION_ENERGY " + ev + " eV" + "\n");
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
					fw.write("AC$CHROMATOGRAPHY: COLUMN_NAME " + column + "\n");
				}
				if(!lc_gradient.isEmpty()) {
					fw.write("AC$CHROMATOGRAPHY: FLOW_GRADIENT " + lc_gradient + "\n");
				}
				if(!lc_flow_rate.isEmpty()) {
					fw.write("AC$CHROMATOGRAPHY: FLOW_RATE " + lc_flow_rate + "\n");
				}
				if(!lc_solvent_a.isEmpty()) {
					fw.write("AC$CHROMATOGRAPHY: SOLVENT A " + lc_solvent_a + "\n");
				}
				if(!lc_solvent_b.isEmpty()) {
					fw.write("AC$CHROMATOGRAPHY: SOLVENT B " + lc_solvent_b + "\n");
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

				Splash splashFactory = SplashFactory.create();
				Spectrum spectrum = new SpectrumImpl(ionList, SpectraType.MS);
				String splash = splashFactory.splashIt(spectrum);
				
				fw.write("PK$SPLASH: "+ splash + "\n");
				
				fw.write("PK$NUM_PEAK: " + numPeaks + "\n");
				fw.write("PK$PEAK: m/z int. rel.int." + "\n");
				fw.write(peaks.toString());
				
				// close FileWriter
				fw.flush();
				fw.close();
				
				// raise counter and create new ID
				counter++;
				id = getPrefix() + formatCounter(counter);
				
				// clear current variables
				peaks = new StringBuffer();
				ionList = new ArrayList<Ion>();
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
				ms = "MS1";
				ev = "";
				ion = "";
				precursor = "";
				date = new SimpleDateFormat("yyyy.MM.dd").format(new Date());
				year = "";
				numPeaks = 0;
				smiles = "";
				inchi = "";
				
				preion = ""; prodion = ""; trapdrive = ""; skim = ""; fragampl = ""; isolwidth = ""; targetgas = ""; targetgaspres = "";
				reagention = ""; reagentgaspres = ""; peakwidth = ""; reflector = ""; psd = ""; chargedecon = ""; column = ""; rettime = "";
				ssid = ""; analid = ""; analname = ""; massrange = "";

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