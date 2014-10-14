package beast.geo;


import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.Set;

import beast.continuous.SphericalDiffusionModel;

public class Vertex {
	double lat1; 
	double long1;
	double [] cart;
	
	static NumberFormat formatter = new DecimalFormat("#0.00"); 
	
	Vertex(double lat1, double long1) {
		this.lat1 = lat1;
		this.long1 = long1;
		adjacentTraingles = new HashSet<Triangle>();
		cart = SphericalDiffusionModel.spherical2Cartesian(lat1, long1);
	}
	
	Set<Triangle> adjacentTraingles;

	@Override
	public String toString() {
		return "(" + formatter.format(lat1) + "," + formatter.format(long1) +")";
	}
}