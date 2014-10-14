package beast.geo;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

import beast.core.BEASTObject;
import beast.core.Description;

@Description("Graph connecting nodes on a surface")
public class Graph extends BEASTObject {
	List<GraphNode> nodes;

	@Override
	public void initAndValidate() throws Exception {
	}

	class DistantGNode {
		
		double distance;
		GraphNode gnode;

		DistantGNode(double distance, GraphNode gnode) {
			this.distance = distance;
			this.gnode = gnode;
		}
		
		@Override
		public String toString() {
			return "(" + gnode.id +":" + Vertex.formatter.format(distance)+")";
		}
	}

	class DistanceGNodeComparator implements Comparator<DistantGNode> {

		@Override
		public int compare(DistantGNode t1, DistantGNode t2) {
			if (t1.distance > t2.distance) {
				return 1;
			} else if (t1.distance < t2.distance) {
				return -1;
			}
			return 0;
		}
	}
	
	public int shortestPath(GraphNode t1, double w1, GraphNode t2, double w2, List<GraphNode> path1, List<GraphNode> path2) {
		final DistanceGNodeComparator comparator = new DistanceGNodeComparator();

		boolean [] done1 = new boolean [nodes.size()];
		double [] dist1 = new double[nodes.size()];
		int[] prev1 = new int[nodes.size()];
		boolean [] done2 = new boolean [nodes.size()]; 
		double [] dist2 = new double[nodes.size()];
		int[] prev2 = new int[nodes.size()];
		
		done1[t1.id] = true;
		prev1[t1.id] = t1.id;
		done2[t2.id] = true;
		prev2[t2.id] = t2.id;
		PriorityQueue<DistantGNode> queue1 = new PriorityQueue<DistantGNode>(comparator);
		queue1.add(new DistantGNode(0, t1));

		PriorityQueue<DistantGNode> queue2 = new PriorityQueue<DistantGNode>(comparator);
		queue2.add(new DistantGNode(0, t2));
		
		while (queue1.size() > 0) {
			DistantGNode ct1 = queue1.peek();
			DistantGNode ct2 = queue2.peek();
			
			if ((ct1.distance / w1) < (ct2.distance / w2)) {
				int connection = doStep(queue1, dist1, prev1, done1, done2);
				if (connection >= 0) {
					buildPaths(connection, t1.id, prev1, path1);
					buildPaths(connection, t2.id, prev2, path2);
					return connection;					
				}
			} else {
				int connection = doStep(queue2, dist2, prev2, done2, done1);
				if (connection >= 0) {
					buildPaths(connection, t1.id, prev1, path1);
					buildPaths(connection, t2.id, prev2, path2);
					return connection;					
				}
			}
		}
		return -1;
	}

	public int shortestPath(GraphNode t1, double w1, GraphNode t2, double w2, GraphNode t3, double w3, List<GraphNode> path1, List<GraphNode> path2, List<GraphNode> path3) {
		final DistanceGNodeComparator comparator = new DistanceGNodeComparator();

		boolean [] done1 = new boolean [nodes.size()];
		double [] dist1 = new double[nodes.size()];
		int[] prev1 = new int[nodes.size()];
		PriorityQueue<DistantGNode> queue1 = new PriorityQueue<DistantGNode>(comparator);
		
		boolean [] done2 = new boolean [nodes.size()]; 
		double [] dist2 = new double[nodes.size()];
		int[] prev2 = new int[nodes.size()];
		PriorityQueue<DistantGNode> queue2 = new PriorityQueue<DistantGNode>(comparator);

		boolean [] done3 = new boolean [nodes.size()]; 
		double [] dist3 = new double[nodes.size()];
		int[] prev3 = new int[nodes.size()];
		PriorityQueue<DistantGNode> queue3 = new PriorityQueue<DistantGNode>(comparator);

		done1[t1.id] = true;
		done2[t2.id] = true;
		done3[t3.id] = true;

		prev1[t1.id] = t1.id;
		prev2[t2.id] = t2.id;
		prev3[t3.id] = t3.id;
		
		queue1.add(new DistantGNode(0, t1));
		queue2.add(new DistantGNode(0, t2));
		queue3.add(new DistantGNode(0, t3));

		while (queue1.size() > 0) {
			DistantGNode ct1 = queue1.peek();
			DistantGNode ct2 = queue2.peek();
			DistantGNode ct3 = queue3.peek();
			
			if ((ct1.distance / w1) <= (ct2.distance / w2) && (ct1.distance / w1) <= (ct3.distance / w3)) {
				int connection = doStep(queue1, dist1, prev1, done1, done2, done3);
				if (connection >= 0) {
					buildPaths(connection, t1.id, prev1, path1);
					buildPaths(connection, t2.id, prev2, path2);
					buildPaths(connection, t3.id, prev3, path3);
					return connection;					
				}
			} else if ((ct2.distance / w2) <= (ct1.distance / w1) && (ct2.distance / w2) <= (ct3.distance / w3)) {
				int connection = doStep(queue2, dist2, prev2, done2, done1, done3);
				if (connection >= 0) {
					buildPaths(connection, t1.id, prev1, path1);
					buildPaths(connection, t2.id, prev2, path2);
					buildPaths(connection, t3.id, prev3, path3);
					return connection;					
				}
			} else {
				int connection = doStep(queue3, dist3, prev3, done3, done1, done2);
				if (connection >= 0) {
					buildPaths(connection, t1.id, prev1, path1);
					buildPaths(connection, t2.id, prev2, path2);
					buildPaths(connection, t3.id, prev3, path3);
					return connection;					
				}
			}
		}
		return -1;
	}

	/** calc all minimal distances from t1 **/
	public double [] distances(GraphNode t1) {
		boolean [] done1 = new boolean [nodes.size()];
		double [] dist1 = new double[nodes.size()];
		int[] prev1 = new int[nodes.size()];
		List<GraphNode> queue = new ArrayList<GraphNode>();
		
		done1[t1.id] = true;
		prev1[t1.id] = t1.id;
		queue.add(t1);
		while (queue.size() > 0) {
			GraphNode gnode = queue.remove(queue.size() - 1);
			double dist = dist1[gnode.id];
			for (int i = 0; i < gnode.neighbours.length; i++) {
				GraphNode t = gnode.neighbours[i];
				double d = gnode.distance[i];
				if (!done1[t.id] || dist + d < dist1[t.id]) {
					dist1[t.id] = dist + d;
					prev1[t.id] = gnode.id;
					done1[t.id] = true;
					queue.add(t);
				}
			}
		}
		return dist1;
	}

	
	int doStep(PriorityQueue<DistantGNode> queue1, double [] dist1, int [] prev1, boolean [] done1, boolean [] done2) {
		GraphNode gnode = queue1.poll().gnode;
		double dist = dist1[gnode.id];
		for (int i = 0; i < gnode.neighbours.length; i++) {
			GraphNode t = gnode.neighbours[i];
			double d = gnode.distance[i];
			if (!done1[t.id] || dist + d < dist1[t.id]) {
				dist1[t.id] = dist + d;
				prev1[t.id] = gnode.id;
				done1[t.id] = true;
				if (done2[t.id]) {
					// we have a connection
					return t.id;
				}
				queue1.add(new DistantGNode(dist + d, t));
			}
		}
		return -1;
	}
	
	int doStep(PriorityQueue<DistantGNode> queue1, double [] dist1, int [] prev1, boolean [] done1, boolean [] done2, boolean [] done3) {
		GraphNode gnode = queue1.poll().gnode;
		double dist = dist1[gnode.id];
		for (int i = 0; i < gnode.neighbours.length; i++) {
			GraphNode t = gnode.neighbours[i];
			double d = gnode.distance[i];
			if (!done1[t.id] || dist + d < dist1[t.id]) {
				dist1[t.id] = dist + d;
				prev1[t.id] = gnode.id;
				done1[t.id] = true;
				if (done2[t.id] && done3[t.id]) {
					// we have a connection
					return t.id;
				}
				queue1.add(new DistantGNode(dist + d, t));
			}
		}
		return -1;
	}

	private void buildPaths(int end, int id1, int[] prev1, List<GraphNode> path1) {
		int i = end;
		path1.add(nodes.get(i));
		do {
			i = prev1[i];
			path1.add(0, nodes.get(i));
		} while (i != id1);
	}
}
