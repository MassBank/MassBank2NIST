/**
 * created by Michael Gerlich, Feb 27, 2012 - 11:11:20 AM
 */ 

package de.ipbhalle.converter;

/**
 * Main file which calls specialized classes in order to perform conversion 
 * between MassBank-to-Library and vice-versa.
 * 
 * @author mgerlich
 *
 */
public class ConverterMain {

	enum switches {mode, spectra, mol, output, prefix, library};
	enum modes {mb2lib, lib2mb};
	private final int NUM_SWITCHES = 4;
	private final String ARG_DELIMITER = "-";
	
	private String mode;
	private String spectraDir;
	private String outputDir;
	private String molDir;
	private String prefix;
	private String libFile;
	
	// -mode mb2lib -spectra /tmp/test/recdata/ -mol /tmp/test/moldata -output /tmp/test/
	
	public ConverterMain() {
		
	}
	
	private void setupArgs(String[] args) {
		if(args.length == NUM_SWITCHES) {
			System.out.println("Using arguments directly!");
			if(args[0].matches(modes.lib2mb.toString())) {
				this.libFile = args[1];
				this.outputDir = args[2];
				this.prefix = args[3];
				String[] newArgs = {libFile, outputDir, prefix};
				LibraryToMassBank.main(newArgs);
			}
			else if(args[0].matches(modes.mb2lib.toString())) {
				this.spectraDir = args[1];
				this.molDir = args[2];
				this.outputDir = args[3];
				String[] newArgs = {spectraDir, molDir, outputDir};
				MassBankToLibrary.main(newArgs);
			}
			else {
				modes[] _modes = modes.values();
				System.out.print("Wrong mode: either use ");
				for (int i = 0; i < _modes.length - 1; i++) {
					System.out.print("[" + _modes[i] + "] or ");
				}
				System.out.println("[" + _modes[_modes.length-1] + "]");
			}
		}
		else if(args.length == 2*NUM_SWITCHES) {
			System.out.println("Using switches and arguments!");
			for (int i = 0; i < args.length; i+=2) {
				if(args[i].startsWith(ARG_DELIMITER))
					args[i] = args[i].substring(1);
				if(switches.valueOf(args[i]) == switches.mode) {
					this.mode = args[i+1];
				}
				else if(switches.valueOf(args[i]) == switches.spectra) {
					this.spectraDir = args[i+1];
				}
				else if(switches.valueOf(args[i]) == switches.mol) {
					this.molDir = args[i+1];
				}
				else if(switches.valueOf(args[i]) == switches.output) {
					this.outputDir = args[i+1];
				}
				else if(switches.valueOf(args[i]) == switches.prefix) {
					this.prefix = args[i+1];
				}
				else if(switches.valueOf(args[i]) == switches.library) {
					this.libFile = args[i+1];
				}
				else {
					System.err.println("Unknown switch!");
					System.err.println(args[i] + " -> " + args[i+1]);
				}
			}
			
			// start program according to mode
			if(mode.equals(modes.mb2lib.toString())) {		// start MassBank-to-Library conversion
				String[] newArgs = {spectraDir, molDir, outputDir};
				MassBankToLibrary.main(newArgs);
			}
			else if(mode.equals(modes.lib2mb.toString())) {	// start Library-to-MassBank conversion
				if(prefix == null || prefix.isEmpty()) {
					System.err.println("Missing argument -prefix !");
					System.err.println("Please specify either 2- or 3-letter code as prefix.");
					System.exit(-1);
				}
				String[] newArgs = {libFile, outputDir, prefix};
				LibraryToMassBank.main(newArgs);
			}
			else {
				System.err.println("Unknown mode definded -> " + mode);
				
			}
		}
		else {
			System.err.println("Wrong number of arguments!");
			System.out.println("[MassBank-to-Library]: -"+switches.mode+" "+modes.mb2lib+"  -"+switches.spectra+" path -"+switches.mol+" path " +
					"-"+switches.output+" path");
			System.out.println("[Library-to-MassBank]: -"+switches.mode+" "+modes.lib2mb+"  -"+switches.library+" file -"+switches.output+" path " +
					"-"+switches.prefix+" prefix");
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// switch for MassBank-to-Library
		// -mode mb2lib -spectra path -mol path -output path
		//String[] test = {"mb2lib", "/tmp/test/recdata/", "/tmp/test/moldata", "/tmp/test/"};
		
		// switch for Library-TO-MassBank, 
		// -mode lib2mb -library path -output path -prefix prefix
		//String[] test = {"lib2mb", "/home/mgerlich/workspace_new/Edda/test_lib/Metabo_Q_Lib_PubChemID.library", "/tmp/test/", "BD"};
		// -mode lib2mb -library /home/mgerlich/workspace_new/Edda/test_lib/Metabo_Q_Lib_PubChemID.library -output /tmp/test/ -prefix BD
		
		ConverterMain con = new ConverterMain();
		con.setupArgs(args);
	}


}
