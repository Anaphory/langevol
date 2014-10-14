package beast.geo;

public class GraphNode {
	int id;

	/** adjacent triangles **/
	GraphNode [] neighbours;

	/** distance to neighbor **/
	double [] distance;

	double [] getCenter() {return null;}
}
