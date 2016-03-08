package beast.geo;


import java.util.HashSet;
import java.util.List;
import java.util.Set;

import sphericalGeo.SphericalDiffusionModel;

public class Quadrangle extends GraphNode {
	
	Vertex v1, v2, v3, v4;

	/** create Quadrangle with corners defined by four vertices **/
	Quadrangle(Vertex v1, Vertex v2, Vertex v3, Vertex v4) {
		this.v1 = v1;
		this.v2 = v2;
		this.v3 = v3;
		this.v4 = v4;
		id = -1;
		v2.adjacentGNodes.add(this);
		v1.adjacentGNodes.add(this);
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
	public double[] getCenter() {
		double [] mean = new double[3];
		for (int i = 0; i < 3; i++) {
			mean[i] = (v1.cart[i] + v2.cart[i] + v3.cart[i] + v4.cart[i]) / 4.0;
		}
		normalise(mean);
		double [] center = SphericalDiffusionModel.cartesian2Sperical(mean, true);
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

	@Override
	void calcNeighbours(boolean allNeighbors, boolean useGreatCircleDistance) {
		Set<GraphNode> neighbourset = new HashSet<GraphNode>();
		
		Set<GraphNode> s = new HashSet<GraphNode>();
		addNeighbors(v1, v2, s, allNeighbors);
		neighbourset.addAll(s);
		addNeighbors(v2, v3, s, allNeighbors);
		neighbourset.addAll(s);
		addNeighbors(v3, v4, s, allNeighbors);
		neighbourset.addAll(s);
		addNeighbors(v4, v1, s, allNeighbors);
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
			neighbours = neighbourset.toArray(new GraphNode[]{});
		}
		setUpDistances(useGreatCircleDistance);
	}
		


	private void addNeighbors(Vertex v2, Vertex v3, Set<GraphNode> s, boolean allNeighbors) {
		s.clear();
		s.addAll(v2.adjacentGNodes);
		if (!allNeighbors)
			s.retainAll(v3.adjacentGNodes);
	}


	@Override
	public void addVertices(Set<Vertex> vertices) {
		vertices.add(v1);
		vertices.add(v2);
		vertices.add(v3);
		vertices.add(v4);
	}

	public void split(List<GraphNode> newQuadrangles) {
		Vertex v12 = getHalfway(v1, v2);
		Vertex v23 = getHalfway(v2, v3);
		Vertex v34 = getHalfway(v3, v4);
		Vertex v41 = getHalfway(v4, v1);
		Vertex center = getHalfway(v1, v2, v3, v4);
		newQuadrangles.add(new Quadrangle(v1, v12, center, v41));
		newQuadrangles.add(new Quadrangle(v2, v23, center, v12));
		newQuadrangles.add(new Quadrangle(v3, v34, center, v23));
		newQuadrangles.add(new Quadrangle(v4, v41, center, v34));
	}

}
