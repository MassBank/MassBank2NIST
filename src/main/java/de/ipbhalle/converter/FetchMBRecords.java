/**
 * created by Michael Gerlich, Jan 31, 2011 - 1:04:57 PM
 */ 

package de.ipbhalle.converter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Map;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.methods.PostMethod;
import org.openscience.cdk.interfaces.IAtomContainer;
import de.ipbhalle.metfusion.utilities.MassBank.MassBankUtilities;
import massbank.MassBankCommon;

public class FetchMBRecords {

	private static String serverUrl = "http://www.massbank.jp/"; 
		//"http://msbi.ipb-halle.de/MassBank/";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MassBankCommon mbCommon = new MassBankCommon();
		
		/**
         * build up parameter string for MassBank search
         */
        String typeName = MassBankCommon.CGI_TBL[MassBankCommon.CGI_TBL_NUM_TYPE][MassBankCommon.CGI_TBL_TYPE_RCDIDX];
        
        // http://www.massbank.jp/jsp/Result.jsp?type=rcdidx&idxtype=inst&srchkey=GC-EI-TOF-MS&sortKey=name&sortAction=1&pageNo=1&exec=
        String param = "idxtype=inst&srchkey=EI-MS&sortKey=name&sortAction=1&pageNo=1&exec=";
        	//"idxtype=inst&srchkey=GC-EI-TOF-MS&sortKey=name&sortAction=1&pageNo=1&exec=";
		/**
		 * 
		 */
		
		// retrieve result list
		ArrayList<String> result = mbCommon.execMultiDispatcher(serverUrl, typeName, param);
		System.out.println(result.size());
		
		for (int i = 0; i < result.size(); i++) {
			System.out.println(i+1 + "/" + result.size());
			
			String s = result.get(i);
            String[] split = s.split("\t");
            if(split.length == 6) {
                // name; instrument
                // id
                // ionization mode
                // sum formula
                // score
                // site
                String name = split[0].substring(0, split[0].indexOf(";"));
                String id = split[1].trim();
                double score = Double.parseDouble(split[4].substring(split[4].indexOf(".")));
                String site = split[5];

            	File f = new File("EI+GC", id + ".txt");
        		if(!f.exists()) {
        			String reqStr = serverUrl;
        			//String reqStr = "http://www.massbank.jp/";
        			reqStr += "jsp/" + MassBankCommon.DISPATCHER_NAME;
        			
        			HttpClient client = new HttpClient();
        			PostMethod method = new PostMethod( reqStr );
        			method.addParameter("type", "disp");		// display record
        			method.addParameter("id", id);				// specify record
        			method.addParameter("site", site);
        			
        			try {
        				client.executeMethod(method);
        				//String result1 = method.getResponseBodyAsString();
        				String result1 = "";
        				
        				// getResponseBodyAsStream()
        				InputStream is = method.getResponseBodyAsStream();
        				StringBuilder sb = new StringBuilder();
        				String line = "";
        				if (is != null) {
        					try {
        						BufferedReader reader = new BufferedReader(
        								new InputStreamReader(is, "UTF-8"));
        						while ((line = reader.readLine()) != null) {
        							if (line.equals("") || line.equals("\n"))
        								sb.append(line);
        							else
        								sb.append(line).append("\n"); // .append("\n");
        						}
        					} finally {
        						is.close();
        					}
        					result1 = sb.toString();
        				}
        				method.releaseConnection();
        				
        				if(result1.contains("ACCESSION")) {
        					String record = result1.substring(result1.indexOf("ACCESSION"), result1.indexOf("</pre>"));
        					FileWriter fw = new FileWriter(f);
        					fw.write(record);
        					fw.flush();
        					fw.close();
        				}
        			}
        			catch(IOException e) {
        				
        			}
        		}

                // create AtomContainer via SMILES
        		MassBankUtilities mbu = new MassBankUtilities();
                //Map<String, String> links = MassBankUtilities.retrieveLinks(id, site);
        		Map<String, String> links = mbu.retrieveLinks(id, site);
                String smiles = links.get("smiles");
                IAtomContainer container = null;
                boolean fetch = false;
                // first look if container is present, then download if not
                //container = MassBankUtilities.getContainer(id, basePath);
                if(container == null) {
                    //fetch = MassBankUtilities.fetchMol(name, id, site, "EI+GC/");
                	fetch = mbu.fetchMol(name, id, site, "EI+GC/");
                    if(fetch) {
                        System.out.println("container via fetch");
                        //container = MassBankUtilities.getMolFromAny(id, basePath, smiles);
                        //container = MassBankUtilities.getContainer(id, "EI+GC");
                        container = mbu.getContainer(id, "EI+GC");
                    }
                    else {
                        System.out.println("container via smiles");
                        //container = MassBankUtilities.getMolFromSmiles(smiles);
                        container = mbu.getMolFromSmiles(smiles);

                        if(container != null) {
                            // write out molfile
                            File mol = new File("EI+GC", id + ".mol");
                            //MassBankUtilities.writeContainer(mol, container);
                            mbu.writeContainer(mol, container);
                        }
                    }
                }
            }
		}
		
	}

}
