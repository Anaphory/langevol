package beast.geo;

import java.util.List;

import beast.continuous.SphericalDiffusionModel;


/** latitudes in -90,90 longitudes in -180,180 **/
public class Triangle {
	double lat1, long1;
	double lat2, long2;
	double lat3, long3;
	double[] cart1;
	double[] cart2;
	double[] cart3;
	
	Triangle(double lat1, double long1, double lat2, double long2, double lat3, double long3) {
		this.lat1 = lat1;
		this.long1 = long1;
		this.lat2 = lat2;
		this.long2 = long2;
		this.lat3 = lat3;
		this.long3 = long3;
		
		cart1 = SphericalDiffusionModel.spherical2Cartesian(lat1, long1);
		cart2 = SphericalDiffusionModel.spherical2Cartesian(lat2, long2);
		cart3 = SphericalDiffusionModel.spherical2Cartesian(lat3, long3);
	}

	/** return centre of triangle in latitude/longitude **/
	public double [] getCenter() {
		double [] mean = new double[3];
		for (int i = 0; i < 3; i++) {
			mean[i] = (cart1[i] + cart2[i] + cart3[i]) / 3.0;
		}
		normalise(mean);
		double [] center = SphericalDiffusionModel.cartesian2Sperical(mean);
		return center;
	}

	/** return halfway point between two points of the triangle in latitude/longitude 
	 * @param leaveout: identify which corner (1,2,3) to leave out
	 * **/
	public double [] getHalfway(int leaveout) {
		double [] mean = new double[3];
		double [] c1;
		double [] c2;
		switch(leaveout) {
		case 1:
				c1 = cart2; c2 = cart3; break;
		case 2:
				c1 = cart1; c2 = cart3; break;
		case 3:
				c1 = cart1; c2 = cart2; break;
		default:
				throw new RuntimeException("leaveout should be one of 1, 2, 3, not " + leaveout);
		}
		
		for (int i = 0; i < 3; i++) {
			mean[i] = (c1[i] + c2[i]) / 2.0;
		}
		normalise(mean);
		double [] center = SphericalDiffusionModel.cartesian2Sperical(mean);
		return center;
	}
	
	public void split(List<Triangle> newTriangles) {
		double[]pos12 = getHalfway(3);
		double[]pos13 = getHalfway(2);
		double[]pos23 = getHalfway(1);
		Triangle t;
		t = new Triangle(lat1, long1, pos12[0], pos12[1], pos13[0], pos13[1]);
		newTriangles.add(t);
		t = new Triangle(lat2, long2, pos12[0], pos12[1], pos23[0], pos23[1]);
		newTriangles.add(t);
		t = new Triangle(lat3, long3, pos13[0], pos13[1], pos23[0], pos23[1]);
		newTriangles.add(t);
		t = new Triangle(pos12[0], pos12[1], pos13[0], pos13[1], pos23[0], pos23[1]);
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
		if (lat1 >= minLat && lat1 <= maxLat && long1 >= minLong && long1 <= maxLong) {
			return true;
		}
		if (lat2 >= minLat && lat2 <= maxLat && long2 >= minLong && long2 <= maxLong) {
			return true;
		}
		if (lat3 >= minLat && lat3 <= maxLat && long3 >= minLong && long3 <= maxLong) {
			return true;
		}
		return false;
	}
}
