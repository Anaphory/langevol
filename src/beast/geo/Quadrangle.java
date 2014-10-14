package beast.geo;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import beast.continuous.SphericalDiffusionModel;

public class Quadrangle extends GraphNode {
	
	Vertex v1, v2, v3, v4;

	/** create Quadrangle with corners defined by four vertices **/
	Quadrangle(Vertex v1, Vertex v2, Vertex v3, Vertex v4) {
		this.v1 = v1;
		this.v2 = v2;
		this.v3 = v3;
		this.v4 = v4;
		id = -1;
		v1.adjacentGNodes.add(this);
		v2.adjacentGNodes.add(this);
		v3.adjacentGNodes.add(this);
		v4.adjacentGNodes.add(this);
	}


	/** create Quadrangle with corners on bounding box defined by two latitudes and two longitudes **/
	public Quadrangle(double minLat, double minLong, double maxLat, double maxLong) {
		this(new Vertex(minLat, minLong),
				new Vertex(minLat, maxLong),
				new Vertex(maxLat, minLong),
				new Vertex(maxLat, maxLong));
	}


	@Override
	double[] getCenter() {
		double [] mean = new double[3];
		for (int i = 0; i < 3; i++) {
			mean[i] = (v1.cart[i] + v2.cart[i] + v3.cart[i] + v4.cart[i]) / 4.0;
		}
		normalise(mean);
		double [] center = SphericalDiffusionModel.cartesian2Sperical(mean);
		return center;
	}

	public void clear() {
		v1.adjacentGNodes = null;
		v2.adjacentGNodes = null;
		v3.adjacentGNodes = null;
		v4.adjacentGNodes = null;
	}

	public boolean hasPointsInside(double minLat, double minLong,
			double maxLat, double maxLong) {
		return
				v1.hasLatLongInsideBBox(minLat, minLong, maxLat, maxLong) ||
				v2.hasLatLongInsideBBox(minLat, minLong, maxLat, maxLong) ||
				v3.hasLatLongInsideBBox(minLat, minLong, maxLat, maxLong) ||
				v4.hasLatLongInsideBBox(minLat, minLong, maxLat, maxLong);
	}

	
	void calcNeighbours() {
		Set<GraphNode> neighbourset = new HashSet<GraphNode>();
		
		Set<GraphNode> s = new HashSet<GraphNode>();
		addNeighbors(v1, v2, s);
		neighbourset.addAll(s);
		addNeighbors(v2, v3, s);
		neighbourset.addAll(s);
		addNeighbors(v3, v4, s);
		neighbourset.addAll(s);
		addNeighbors(v4, v1, s);
		neighbourset.addAll(s);

		neighbourset.remove(this);
		neighbours = neighbourset.toArray(new Quadrangle[]{});
		
		// remove graphnodes outside bounding box
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


	private void addNeighbors(Vertex v12, Vertex v22, Set<GraphNode> s) {
		s.clear();
		s.addAll(v2.adjacentGNodes);
		s.retainAll(v3.adjacentGNodes);
	}

}
