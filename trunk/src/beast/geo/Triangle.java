package beast.geo;

import java.util.List;

import beast.continuous.SphericalDiffusionModel;


/** latitudes in -90,90 longitudes in -180,180 **/
public class Triangle {
	Vertex v1, v2, v3;
//	double lat1, long1;
//	double lat2, long2;
//	double lat3, long3;
//	double[] cart1;
//	double[] cart2;
//	double[] cart3;
	Triangle [] neightbours;
	
//	Triangle(double lat1, double long1, double lat2, double long2, double lat3, double long3) {
//		this.lat1 = lat1;
//		this.long1 = long1;
//		this.lat2 = lat2;
//		this.long2 = long2;
//		this.lat3 = lat3;
//		this.long3 = long3;
//		
//		cart1 = SphericalDiffusionModel.spherical2Cartesian(lat1, long1);
//		cart2 = SphericalDiffusionModel.spherical2Cartesian(lat2, long2);
//		cart3 = SphericalDiffusionModel.spherical2Cartesian(lat3, long3);
//	}

	Triangle(Vertex v1, Vertex v2, Vertex v3) {
		this.v1 = v1;
		this.v2 = v2;
		this.v3 = v3;
	}

	/** return centre of triangle in latitude/longitude **/
	public double [] getCenter() {
		double [] mean = new double[3];
		for (int i = 0; i < 3; i++) {
			mean[i] = (v1.cart[i] + v2.cart[i] + v3.cart[i]) / 3.0;
		}
		normalise(mean);
		double [] center = SphericalDiffusionModel.cartesian2Sperical(mean);
		return center;
	}

	/** return halfway point between two points of the triangle in latitude/longitude 
	 * @param leaveout: identify which corner (1,2,3) to leave out
	 * **/
	public Vertex getHalfway(int leaveout) {
		double [] mean = new double[3];
		double [] c1;
		double [] c2;
		switch(leaveout) {
		case 1:
				c1 = v2.cart; c2 = v3.cart; break;
		case 2:
				c1 = v1.cart; c2 = v3.cart; break;
		case 3:
				c1 = v1.cart; c2 = v2.cart; break;
		default:
				throw new RuntimeException("leaveout should be one of 1, 2, 3, not " + leaveout);
		}
		
		for (int i = 0; i < 3; i++) {
			mean[i] = (c1[i] + c2[i]) / 2.0;
		}
		normalise(mean);
		double [] center = SphericalDiffusionModel.cartesian2Sperical(mean);
		return new Vertex(center[0], center[1]);
	}
	
	public void split(List<Triangle> newTriangles) {
		Vertex pos12 = getHalfway(3);
		Vertex pos13 = getHalfway(2);
		Vertex pos23 = getHalfway(1);
		Triangle t;
		t = new Triangle(v1, pos12, pos13);
		newTriangles.add(t);
		t = new Triangle(v2, pos13, pos23);
		newTriangles.add(t);
		t = new Triangle(v3, pos13, pos23);
		newTriangles.add(t);
		t = new Triangle(pos12, pos13, pos23);
		newTriangles.add(t);
	}

	private void normalise(double[] position) {
		double len = Math.sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
		position[0] /= len;
		position[1] /= len;
		position[2] /= len;
	}

	/** assumes that the triangle is much smaller than the bounding box, that is the bounding box is not inside the triangle **/
	public boolean hasPointsInside(double minLat, double minLong,
			double maxLat, double maxLong) {
		if (v1.lat1 >= minLat && v1.lat1 <= maxLat && v1.long1 >= minLong && v1.long1 <= maxLong) {
			return true;
		}
		if (v2.lat1 >= minLat && v2.lat1 <= maxLat && v2.long1 >= minLong && v2.long1 <= maxLong) {
			return true;
		}
		if (v3.lat1 >= minLat && v3.lat1 <= maxLat && v3.long1 >= minLong && v3.long1 <= maxLong) {
			return true;
		}
		return false;
	}
}
