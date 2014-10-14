package beast.geo;


import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.Set;

import beast.continuous.SphericalDiffusionModel;

public class Vertex {
	static NumberFormat formatter = new DecimalFormat("#0.00"); 

	double lat1; 
	double long1;
	double [] cart;
	
	Set<GraphNode> adjacentGNodes;
	
	Vertex(double lat1, double long1) {
		this.lat1 = lat1;
		this.long1 = long1;
		adjacentGNodes = new HashSet<GraphNode>();
		cart = SphericalDiffusionModel.spherical2Cartesian(lat1, long1);
	}
	
	
	boolean hasLatLongInsideBBox(double minLat, double minLong,
			double maxLat, double maxLong) {
		if (lat1 >= minLat && lat1 <= maxLat && long1 >= minLong && long1 <= maxLong) {
			return true;
		}
		return false;
	}

	@Override
	public String toString() {
		return "(" + formatter.format(lat1) + "," + formatter.format(long1) +")";
	}
}