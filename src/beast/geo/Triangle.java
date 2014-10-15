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
		v1.adjacentGNodes.add(this);
		v2.adjacentGNodes.add(this);
		v3.adjacentGNodes.add(this);
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
	
	public void split(List<GraphNode> newTriangles) {
		Vertex pos12 = getHalfway(v1, v2);
		Vertex pos13 = getHalfway(v1, v3);
		Vertex pos23 = getHalfway(v2, v3);
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


	/** assumes that the triangle is much smaller than the bounding box, that is the bounding box is not inside the triangle **/
	public boolean hasPointsInside(double minLat, double minLong,
			double maxLat, double maxLong) {
		return
				v1.hasLatLongInsideBBox(minLat, minLong, maxLat, maxLong) ||
				v2.hasLatLongInsideBBox(minLat, minLong, maxLat, maxLong) ||
				v3.hasLatLongInsideBBox(minLat, minLong, maxLat, maxLong);
//		
//		if (v1.lat1 >= minLat && v1.lat1 <= maxLat && v1.long1 >= minLong && v1.long1 <= maxLong) {
//			return true;
//		}
//		if (v2.lat1 >= minLat && v2.lat1 <= maxLat && v2.long1 >= minLong && v2.long1 <= maxLong) {
//			return true;
//		}
//		if (v3.lat1 >= minLat && v3.lat1 <= maxLat && v3.long1 >= minLong && v3.long1 <= maxLong) {
//			return true;
//		}
//		return false;
	}
	
	void calcNeighbours() {
		Set<GraphNode> neighbourset = new HashSet<GraphNode>();
		
		Set<GraphNode> s = new HashSet<GraphNode>();
		s.addAll(v1.adjacentGNodes);
		s.retainAll(v2.adjacentGNodes);
		neighbourset.addAll(s);

		s.clear();
		s.addAll(v2.adjacentGNodes);
		s.retainAll(v3.adjacentGNodes);
		neighbourset.addAll(s);
		
		s.clear();
		s.addAll(v3.adjacentGNodes);
		s.retainAll(v1.adjacentGNodes);
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
		v1.adjacentGNodes = null;
		v2.adjacentGNodes = null;
		v3.adjacentGNodes = null;
	}

	@Override
	public void addVertices(Set<Vertex> vertices) {
		vertices.add(v1);
		vertices.add(v2);
		vertices.add(v3);
	}
	
}
