package beast.geo;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;

import beast.continuous.SphericalDiffusionModel;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.util.Randomizer;

@Description("Tesselates (part of) a sphere with equal sized quadrangle on a grid")
public class GridTesselation extends SphereTesselation {

		public Input<File> mapFileInput = new Input<File>("map","geographic world map in Mercator projection -- red coloured items will have distance defined by 'reddistance'. "
				+ "If not specified, all distances are 1.0",
				new File("/home/remco/data/geo/aboriginal25.bmp"));
		public Input<Double> reddistanceInput = new Input<Double>("reddistance","distance fpr red coloured items", 1.0);
		
		GraphNode [][] grid;
		
		/** flag indicating we use native code **/
		boolean m_bNative = false;
		
		public GridTesselation() {
			bboxInput.setRule(Validate.REQUIRED);
		}
	
		@Override
		public void initAndValidate() throws Exception {
			
			parseBBox();
			Quadrangle q = new Quadrangle(minLat, minLong, maxLat, maxLong);
			double [] center = q.getCenter();
			
			int depth = depthInput.get();
			if (depth <= 0) {
				throw new RuntimeException("depth must be a positive number");
			}
			
			// create vertices
			double long0 = minLong - center[1];
			double longStep = (maxLong - minLong)/depth;
			double lat0 = minLat - center[0];
			double latStep = (maxLat - minLat)/depth;
			Vertex [][] vertex = new Vertex[depth][depth];
			for (int i = 0; i < depth; i++) {
				double lat1 = lat0 + latStep * i;
				for (int j = 0; j < depth; j++) {
					double long1 = long0 + longStep * j;
					double [] point = SphericalDiffusionModel.reverseMap(lat1, long1, center[0], center[1]);
					vertex[i][j] = new Vertex(point[0], point[1]);
				}
			}
			
			grid = new GraphNode[depth-1][depth-1];
			nodes = new ArrayList<GraphNode>();
			for (int i = 0; i < depth - 1; i++) {
				for (int j = 0; j < depth - 1; j++) {
					Quadrangle q1 = new Quadrangle(vertex[i][j],vertex[i+1][j], vertex[i+1][j+1],vertex[i][j+1]);
					nodes.add(q1);
					grid[i][j] = q1;
				}
			}
			
			
			boolean all = allNeighborsInput.get();
			for (int i = 0; i < depth - 1; i++) {
				for (int j = 0; j < depth - 1; j++) {
					List<GraphNode> neighbors = new ArrayList<GraphNode>();
					if (i > 0) {
						neighbors.add(grid[i-1][j]);
					}
					if (all && i > 0 && j < depth-2) {
						neighbors.add(grid[i-1][j + 1]);
					}
					if (j < depth-2) {
						neighbors.add(grid[i][j+1]);
					}
					if (all && i < depth - 2 && j < depth-2) {
						neighbors.add(grid[i+1][j + 1]);
					}
					if (i < depth-2) {
						neighbors.add(grid[i+1][j]);
					}
					if (all && i < depth - 2 && j > 0) {
						neighbors.add(grid[i+1][j - 1]);
					}
					if (j > 0) {
						neighbors.add(grid[i][j-1]);
					}
					if (all && i > 0  && j > 0) {
						neighbors.add(grid[i-1][j - 1]);
					}
					grid[i][j].neighbours = neighbors.toArray(new GraphNode[]{});
					grid[i][j].setUpDistances(useGreatCircleInput.get());
				}
			}
			
			System.err.println("#nodes = " + nodes.size());
	
			// renumber remaining quadrangles
			renumber();
	
			// adjust distances to map
			final BufferedImage image = ImageIO.read(mapFileInput.get());
			int w = image.getWidth();
			int h = image.getHeight();

			double redDistance = reddistanceInput.get();
			for (GraphNode t : nodes) {
				double [] c = t.getCenter();
				int x =(int)( w * (c[1]+180) / 360.0);
				int y =(int)( h * (c[0]+90) / 180.0);
				int color = image.getRGB(x, y) & 0xFFFFFF;
				if (color == 0x00FF00) {
					t.scaleDistance(redDistance);
				}
			}

			setUpLatLongMap();
			
			// log some stats
			System.err.println("#nodes= " + nodes.size());
		}
	
	
		/** return distance between node (x1, y1) and neighbouring node nr (x2, y2) **/
		double getDistance(int x1, int y1, int x2, int y2) {
			GraphNode v1 = grid[x1][y1];
			int depth = grid.length;
			if (x2 < 0 || x2 >= depth || y2 < 0 || y2 >= depth) {
				return Double.POSITIVE_INFINITY;
			}
			if (Math.abs(x1 - x2) > 1 || Math.abs(y1 - y2) > 1 ) {
				throw new RuntimeException("node (x1, y1) must be neighbouring node (x2, y2)");
			}
			GraphNode v2 = grid[x2][y2];
			int k = 0;
			for (GraphNode v : v1.neighbours) {
				if (v2.id == v.id) {
					return v1.distance[k];
				}
				k++;
			}
			throw new RuntimeException("node (x2, y2) is not in neighbor set of (x1, y1)");
		}
		
		public enum DISTANCE_METHOD {by_graph, parallel, by_native}; 
		
		public double[] distances(GraphNode t1, DISTANCE_METHOD method) {
			switch(method) { 
			case by_graph : 
				return super.distances(t1);
			case parallel : 
				return parallelDistances(t1);
			case by_native: 
				return nativeDistances(t1);
			}
			return null;
		}

		private double[] parallelDistances(GraphNode t1) {
			if (!allNeighborsInput.get()) {
				throw new RuntimeException("parallelDistances() assumes allNeighbours = true");
			}
			int depth = grid.length;
			double [] distances = new double[(depth+2)*(depth+2)];
			Arrays.fill(distances, 1e100);
			int id = t1.id;
			distances[(id/depth + 1) * (depth+2) + (id%depth) + 1] = 0;
			double [] pairwise = getPairWise();

			double [] distances2 = new double[(depth+2)*(depth+2)];
			System.arraycopy(distances, 0, distances2, 0, distances.length);
			
			for (int step = 0; step < 50; step++) {
				for (int i = 0; i < depth; i++) {
					for (int j = 0; j < depth; j++) {
						int k = ((i+1) * (depth + 2) + j +1);
						double d = distances[k];
						d = Math.min(d, distances[k - depth - 2] + pairwise[k*8]);
						d = Math.min(d, distances[k - depth - 1] + pairwise[k*8+1]);
						d = Math.min(d, distances[k + 1] + pairwise[k*8+2]);
						d = Math.min(d, distances[k + depth + 3] + pairwise[k*8+3]);
						d = Math.min(d, distances[k + depth + 2] + pairwise[k*8+4]);
						d = Math.min(d, distances[k + depth + 1] + pairwise[k*8+5]);
						d = Math.min(d, distances[k - 1] + pairwise[k*8+6]);
						d = Math.min(d, distances[k - depth - 3] + pairwise[k*8+7]);
						distances2[k] = d;
					}
				}
				double [] tmp = distances; 
				distances = distances2;
				distances2 = tmp;
			}
			
			
			return distances;
		}
		
		double [] getPairWise() {
			int depth = grid.length;
			double [] pairwise = new double[(depth+2)*(depth+2) * 8];
			Arrays.fill(pairwise, 1e100);
			for (int i = 0; i < depth; i++) {
				for (int j = 0; j < depth; j++) {
					int k = ((i+1) * (depth + 2) + j +1)* 8;
					pairwise[k++] = getDistance(i, j, i-1, j);
					pairwise[k++] = getDistance(i, j, i-1, j + 1);
					pairwise[k++] = getDistance(i, j, i, j+1);
					pairwise[k++] = getDistance(i, j, i+1, j + 1);
					pairwise[k++] = getDistance(i, j, i+1, j);
					pairwise[k++] = getDistance(i, j, i+1, j - 1);
					pairwise[k++] = getDistance(i, j, i, j-1);
					pairwise[k++] = getDistance(i, j, i-1, j - 1);
				}
			}
			return pairwise;
		}

		native void setDevice(int device);
		native double [] doNativeDistances(double [] distances, double [] pairwise, int depth, int iterations);
		native float [] doNativeDistancesSinglePrecission(float [] distances, float[] pairwise, int depth, int iterations);

		private double[] nativeDistances(GraphNode t1) {
			if (!allNeighborsInput.get()) {
				throw new RuntimeException("parallelDistances() assumes allNeighbours = true");
			}
			int depth = grid.length;
			double [] distances = new double[(depth+2)*(depth+2)];
			Arrays.fill(distances, 1e100);
			int id = t1.id;
			distances[(id/depth + 1) * (depth+2) + (id%depth) + 1] = 0;
			double [] pairwise = getPairWise();

			distances = doNativeDistances(distances, pairwise, depth, 50);
			
			return distances;
		}

		
		void loadLibrary() {
			try {
				System.loadLibrary("Tesselate");
				m_bNative = true;
				//setDevice(1);
				System.err.println("Going native!!!");
			} catch (Throwable e) {
				System.err.println(e.getMessage());
				System.err.println("Going slow!!!");
			}
		}
		
		public static void main(String[] args) throws Exception {
			final GridTesselation tessel = new GridTesselation();
			tessel.loadLibrary();
			tessel.initByName("depth", 50, "bbox", "5 120 47 154", "reddistance", "1.0", 
					"allNeighbors", true, "useGreatCircle", true);
			
			int n = tessel.nodes.size();
			Randomizer.setSeed(122L);
			final GraphNode t1 = tessel.getLowerLeftCorner();//tessel.nodes.get(Randomizer.nextInt(n));
			//double [] dist = new double [] {0.0};

			
			final double [] dist = tessel.distances(t1, DISTANCE_METHOD.by_graph);
			final double [] dist2 = tessel.distances(t1, DISTANCE_METHOD.by_native);
			double max = 0;
			int depth = tessel.grid.length;
			for (int i = 0; i < depth; i++) {
				for (int j = 0; j < depth; j++) {
					max = Math.max(max, Math.abs(dist[i * depth + j] - dist2[(i+1)*(depth+2) + j + 1]));
				}
			}
			System.err.println("max diff = " + max);

		}
}
