package de.ipbhalle.wrapper;

import java.util.List;

import org.openscience.cdk.interfaces.IAtomContainer;

public class Container {

	private String id;
	private String moldata;
	private List<IAtomContainer> containers;
	private boolean isomorph;
	
	public Container(String id, String moldata, List<IAtomContainer> containers, boolean isomorph) {
		this.setId(id);
		this.setMoldata(moldata);
		this.setContainers(containers);
		this.setIsomorph(isomorph);
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getId() {
		return id;
	}

	public void setMoldata(String moldata) {
		this.moldata = moldata;
	}

	public String getMoldata() {
		return moldata;
	}

	public void setContainers(List<IAtomContainer> containers) {
		this.containers = containers;
	}

	public List<IAtomContainer> getContainers() {
		return containers;
	}

	public void setIsomorph(boolean isomorph) {
		this.isomorph = isomorph;
	}

	public boolean isIsomorph() {
		return isomorph;
	}
}
