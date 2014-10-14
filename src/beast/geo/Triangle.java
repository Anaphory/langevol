package beast.geo;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.continuous.SphericalDiffusionModel;


/** latitudes in -90,90 longitudes in -180,180 **/
public class Triangle extends GraphNode {
	
	Vertex v1, v2, v3;

	Triangle(Vertex v1, Vertex v2, Vertex v3) {
		this.v1 = v1;
		this.v2 = v2;
		this.v3 = v3;
		id = -1;
		v1.adjacentTraingles.add(this);
		v2.adjacentTraingles.add(this);
		v3.adjacentTraingles.add(this);
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
	
	public void split(List<GraphNode> newTriangles) {
		Vertex pos12 = getHalfway(3);
		Vertex pos13 = getHalfway(2);
		Vertex pos23 = getHalfway(1);
		Triangle t;
		t = new Triangle(v1, pos12, pos13);
		newTriangles.add(t);
		t = new Triangle(v2, pos12, pos23);
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
	
	void calcNeighbours() {
		Set<Triangle> neighbourset = new HashSet<Triangle>();
		
		Set<Triangle> s = new HashSet<Triangle>();
		s.addAll(v1.adjacentTraingles);
		s.retainAll(v2.adjacentTraingles);
		neighbourset.addAll(s);

		s.clear();
		s.addAll(v2.adjacentTraingles);
		s.retainAll(v3.adjacentTraingles);
		neighbourset.addAll(s);
		
		s.clear();
		s.addAll(v3.adjacentTraingles);
		s.retainAll(v1.adjacentTraingles);
		neighbourset.addAll(s);

		neighbourset.remove(this);
		neighbours = neighbourset.toArray(new Triangle[]{});
		
		// remove triangles outside bounding box
		// these can be identified since they have id = -1
		boolean b = false;
		for (int i = 0; i < neighbours.length; i++) {
			if (neighbours[i].id < 0) {
				neighbourset.remove(neighbours[i]);
				b = true;
			}
		}
		if (b) {
			neighbours = neighbourset.toArray(new Triangle[]{});
		}
		
		distance = new double[neighbours.length];
		Arrays.fill(distance, 1.0);
	}

	
	public void clear() {
		v1.adjacentTraingles = null;
		v2.adjacentTraingles = null;
		v3.adjacentTraingles = null;
	}
	
	@Override
	public String toString() {
		return id + "";//"<" + v1.toString() + ", " + v2.toString() + ", " + v3.toString() +">";
	}
}
