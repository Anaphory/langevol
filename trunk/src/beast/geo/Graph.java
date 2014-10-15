package beast.geo;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.imageio.ImageIO;

import beast.core.BEASTObject;
import beast.core.Description;

@Description("Graph connecting nodes on the earth's surface")
public class Graph extends BEASTObject {
	List<GraphNode> nodes;
	
	DistanceMatrix distances;

	LatLongMap latLongMap;
	
	@Override
	public void initAndValidate() throws Exception {
	}
	
	/** used to quickly find a GraphNode closest to a given latitude/longitude pair 
	 * in a graph **/
	class LatLongMap {
		double maxLat = 90, minLat = -90;
		double maxLong= 180, minLong= -180;

		double deltaLat = Double.POSITIVE_INFINITY;
		double deltaLong = Double.POSITIVE_INFINITY;

		GraphNode [][] map;
		
		LatLongMap(List<GraphNode> nodes) {
			final int X = 20; 
			final int MULTI_X = 2*X/3; 
			
			Set<Vertex> vertices = new HashSet<Vertex>();
			for (GraphNode t : nodes) {
				t.addVertices(vertices);
			}
			
			// determine map boundary
			Vertex dummy = (Vertex) vertices.toArray()[0];
			maxLat = dummy.lat1;
			minLat = dummy.lat1;
			maxLong = dummy.long1;
			minLong = dummy.long1;

			for (Vertex v : vertices) {
				maxLat = Math.max(maxLat, v.lat1);
				minLat = Math.min(minLat, v.lat1);
				maxLong = Math.max(maxLong, v.long1);
				minLong = Math.min(minLong, v.long1);
			}

			// determine step size
			for (GraphNode node : nodes) {
				double [] center = node.getCenter();
				for (GraphNode nb : node.neighbours) {
					double [] nbcenter = nb.getCenter();
					double dLat = Math.abs(center[0] - nbcenter[0]);
					double dLong = Math.abs(center[1] - nbcenter[1]);
					if (dLat < dLong) {
						deltaLong = Math.min(deltaLong, dLong);
					} else {
						deltaLat = Math.min(deltaLat, dLat);
					}
				}
			}
			
			// set up the map
			int latSteps =  X * (int)((maxLat - minLat + deltaLat * 0.9999) / deltaLat);
			int longSteps = X * (int)((maxLong - minLong + deltaLong * 0.9999) / deltaLong);
			
//			latSteps = 90*10;
//			longSteps = 180*10;
//			deltaLat = 1;
//			deltaLong = 1;
//			minLong = 0;
//			minLat = 0;
//			maxLong = 180;
//			maxLat = 90;

			deltaLat = (maxLat - minLat) / latSteps;
			deltaLong = (maxLong - minLong) / longSteps;
			map = new GraphNode[latSteps][longSteps];

			double [][] centers = new double[nodes.size()][];
			for (GraphNode node : nodes) {
				if (centers[node.id] == null) { 
					centers[node.id] = node.getCenter();
				} else {
					System.err.println("duplicate id found " + node.id + " " + centers[node.id]);
				}
			}

			for (GraphNode node : nodes) {
				double [] center = centers[node.id];
				int iLat = (int)((center[0] - minLat + deltaLat/2) / deltaLat);
				int iLong = (int)((center[1] - minLong + deltaLong/2) / deltaLong);
				setMap(iLat, iLong, node);
				
				// check whether surrounding map locations are also closest to node
				// in a MULTI_XxMULTI_X grid with (iLat,iLong) at the center
				for (int i = -MULTI_X; i <= MULTI_X; i++) {
					int y = iLat + i;
					if (y >= 0 && y < latSteps) {
						double lat0 = minLat + y * deltaLat + deltaLat/2;
						for (int j = - MULTI_X; j <= MULTI_X; j++) {
							int x = iLong + j;
							if (x >= 0 && x < longSteps) {
								double long0 = minLong + x * deltaLong + deltaLong/2;
								double nodeDist = (center[0] - lat0) * (center[0] - lat0) + (center[1] - long0) * (center[1] - long0);
								
								boolean foundCloser = false;
								double currentDist = nodeDist;
								for (GraphNode nb : node.neighbours) {
									double [] nbcenter = centers[nb.id];
									double nbDist = (nbcenter[0] - lat0) * (nbcenter[0] - lat0) + (nbcenter[1] - long0) * (nbcenter[1] - long0);
									if (nbDist < currentDist) {
										if (map[y][x] != null) {
											double [] occupyingcenter = centers[map[y][x].id];
											double ocDist = (occupyingcenter[0] - lat0) * (occupyingcenter[0] - lat0) + (occupyingcenter[1] - long0) * (occupyingcenter[1] - long0);
											if (ocDist > currentDist) {
												if (i == 0 && j == 0) {
													int iLat2 = (int)((nbcenter[0] - minLat) / deltaLat);
													int iLong2 = (int)((nbcenter[1] - minLong) / deltaLong);
													System.err.println(iLat2 + " " + iLong2);
												}
												setMap(y, x, nb);
											}
										} else {
											setMap(y, x, nb);
										}
										foundCloser = true;
										currentDist = nbDist;
									}
								}
								if (!foundCloser) {
									if (map[y][x] != null) {
										double [] occupyingcenter = centers[map[y][x].id];
										double ocDist = (occupyingcenter[0] - lat0) * (occupyingcenter[0] - lat0) + (occupyingcenter[1] - long0) * (occupyingcenter[1] - long0);
										if (ocDist > currentDist) {
											setMap(y, x, node);
										}
									} else {
										setMap(y, x, node);
									}
								}
								
							}
						}
					}
				}
			}
			
			// sanity check
			for (GraphNode node : nodes) {
				GraphNode closest = getClosestNode(centers[node.id][0], centers[node.id][1]);
				if (node.id != closest.id) {
					throw new RuntimeException("node is not closest to itself");
				}
			}
			// sanity check 2
			int [] bm = new int[map.length * map[0].length];
			int k = 0;
			for (int i = 0; i < map.length; i++) {
				for (int j = 0; j < map[0].length; j++) {
					GraphNode node = map[i][j];
					if (node != null) {
						int id = node.id;
						bm[k++] = Color.HSBtoRGB(((id % 19)/19.0f), 1.0f, 1.0f);
					} else {
						bm[k++] = 0;
					}
				}
			}
			BufferedImage bMap = new BufferedImage(longSteps, latSteps, BufferedImage.TYPE_INT_RGB);
			bMap.setRGB(0, 0, longSteps, latSteps, bm, 0, longSteps);
			try {
				ImageIO.write(bMap, "png", new File("/tmp/graphmap.png"));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		private void setMap(int iLat, int iLong, GraphNode node) {
//			if (iLat == 370 && iLong == 1396) {
//				int k = 3;
//				k++;
//			}
			map[iLat][iLong] = node;
		}

		GraphNode getClosestNode(double latitude, double longitude) {
			int iLat = (int)((latitude - minLat + deltaLat/2) / deltaLat);
			int iLong = (int)((longitude - minLong + deltaLong/2) / deltaLong);
			return map[iLat][iLong];
		}
	} // class LatLongMap
	

	void setUpLatLongMap() {
		latLongMap = new LatLongMap(nodes);
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
	} // class DistantGNode

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
	} // class DistanceGNodeComparator
	
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
	
	/** calc all pairwise distances **/
	public double[][] distances() {
		long start = System.currentTimeMillis();
		
		double [][] distances = new double[nodes.size()][];
		for (int i = 0; i < nodes.size(); i++) {
			distances[i] = distances(nodes.get(i));
			if (i % 10 == 9) {
				long end = System.currentTimeMillis();
				System.err.println((end-start)/1000 + " sec " + ((nodes.size()-i) * (end-start)/(1000 * (i+1))) +" sec to go");
			}
		}
		return distances;
	}
	
	
	/** threaded version of distances() method **/
	public DistanceMatrix distances(int threads) {
		double [][] distances = new double[nodes.size()][];
		final CountDownLatch m_nCountDown = new CountDownLatch(threads);

		/** calculates distance for nodes from to too **/
		class DistanceRunner implements Runnable {
			int from, too;
			double[][] distances;
			long start;
			
			public DistanceRunner(double[][] distances, int from, int too, long start) {
				this.distances = distances;
				this.from = from;
				this.too = too;
				this.start = start;
			}

			@Override
			public void run() {
				for (int i = from; i < too; i++) {
					if (i < distances.length) {
						distances[i] = distances(nodes.get(i));
					}
					if (i % 10 == 9) {
						long end = System.currentTimeMillis();
						System.err.println("Thread " + (from/(too-from)) + "  Done " + (i - from + 1) + " to go " + (too - i) + " " + (end-start)/1000 + " sec " + ((too-i) * (end-start)/(1000 * (i-from+1))) +" sec to go");
					}
				}
				System.err.println("Done " + from + " to " + too);
	  		    m_nCountDown.countDown();
			}
		}
		
		/** manage the threads **/
		ExecutorService g_exec = Executors.newFixedThreadPool(threads);
		long start = System.currentTimeMillis();
		
		int stepsize = (nodes.size() + threads - 1)/ threads;
		for (int i = 0; i < threads; i++) {
    		DistanceRunner coreRunnable = new DistanceRunner(distances, i * stepsize, (i+1) * stepsize, start);
    		g_exec.execute(coreRunnable);
    	}
		try {
			m_nCountDown.await();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		long end = System.currentTimeMillis();
		System.err.println(" Done in "  + (end-start)/1000 + " sec ");
		
		return new DistanceMatrix(distances);
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

	public double getDistance(double[] start, double[] stop) {
		if (latLongMap == null) {
			throw new RuntimeException("call setUpLatLongMap() before calling getDistance()");
		}
		GraphNode startNode = latLongMap.getClosestNode(start[0], start[1]);
		GraphNode endNode = latLongMap.getClosestNode(stop[0], stop[1]);
		double distance = distances.distances[startNode.id][endNode.id];
		return distance;
	}
}
