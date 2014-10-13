package beast.geo;

import java.util.ArrayList;
import java.util.List;

import beast.continuous.SphericalDiffusionModel;

public class Vertex {
	double lat1; 
	double long1;
	double [] cart;
	
	Vertex(double lat1, double long1) {
		this.lat1 = lat1;
		this.long1 = long1;
		traingles = new ArrayList<Triangle>();
		cart = SphericalDiffusionModel.spherical2Cartesian(lat1, long1);
	}
	
	List<Triangle> traingles;

}