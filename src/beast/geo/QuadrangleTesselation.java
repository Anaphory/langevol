package beast.geo;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;

@Description("Tesselates (part of) a sphere with equal sized quadrangle")
public class QuadrangleTesselation extends SphereTesselation {

	@Override
	public void initAndValidate() throws Exception {
		nodes = new ArrayList<GraphNode>();


		System.err.println("#nodes = " + nodes.size() + " before filtering");
		filterNodesInBoundingBox();

		List<Vertex> vertices = amalgamateNeighbors();

		// renumber remaining quadrangles
		renumber();

		// set up adjacency graph -- requires vertices to have adjacentGraphNode
		// to be set up
		for (GraphNode t : nodes) {
			((Triangle) t).calcNeighbours();
		}

		// save memory
		for (Vertex v : vertices) {
			v.adjacentGNodes = null;
		}

		// log some stats
		System.err.println("#nodes= " + nodes.size());
	}

}
