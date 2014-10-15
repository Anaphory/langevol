package beast.geo;


import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

import beast.continuous.SphericalDiffusionModel;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.util.Randomizer;

@Description("Tesselates (part of) a sphere with equal sized quadrangle")
public class QuadrangleTesselation extends SphereTesselation {
	public Input<File> mapFileInput = new Input<File>("map","geographic world map in Mercator projection -- red coloured items will have distance defined by 'reddistance'. "
			+ "If not specified, all distances are 1.0",
			new File("/home/remco/data/geo/aboriginal25.bmp"));
	public Input<Boolean> removeOveWaterInput = new Input<Boolean>("removeWater","remove items over water as defined in map", true);
	public Input<Double> reddistanceInput = new Input<Double>("reddistance","distance fpr red coloured items", 1.0);
	
	QuadrangleTesselation() {
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
		
		nodes = new ArrayList<GraphNode>();
		for (int i = 0; i < depth - 1; i++) {
			for (int j = 0; j < depth - 1; j++) {
				Quadrangle q1 = new Quadrangle(vertex[i][j],vertex[i+1][j], vertex[i+1][j+1],vertex[i][j+1]);
				nodes.add(q1);
			}
		}
		

		System.err.println("#nodes = " + nodes.size() + " before filtering");
//		filterNodesInBoundingBox();


		// filter nodes that are in the sea only
		List<GraphNode> filteredTriangles = new ArrayList<>();
		final BufferedImage image = ImageIO.read(mapFileInput.get());
		int w = image.getWidth();
		int h = image.getHeight();
		if (removeOveWaterInput.get()) {
			for (GraphNode t : nodes) {
				double [] c = t.getCenter();
				int x =(int)( w * (c[1]+180) / 360.0);
				int y =(int)( h * (c[0]+90) / 180.0);
				int color = image.getRGB(x, y) & 0xFFFFFF;
				if (color != 0xFF0000) {
					filteredTriangles.add(t);
				}
			}
			nodes = filteredTriangles;
		}
		
		// collect vertices
		List<Vertex> vertices = new ArrayList<Vertex>();
		for (GraphNode t : nodes) {
			vertices.add(((Quadrangle)t).v1);
			vertices.add(((Quadrangle)t).v2);
			vertices.add(((Quadrangle)t).v3);
			vertices.add(((Quadrangle)t).v4);
		}
		amalgamateNeighborsInVertices(vertices);

		// renumber remaining quadrangles
		renumber();

		// set up adjacency graph -- requires vertices to have adjacentGraphNode
		// to be set up
		for (GraphNode t : nodes) {
			t.calcNeighbours();
		}
		
		// only use those components that are connected to center
		GraphNode centerNode = nodes.get(0);
		double minDist = Math.abs(center[0] - ((Quadrangle) centerNode).v1.lat1) + Math.abs(center[1] - ((Quadrangle) centerNode).v1.long1);
		for (GraphNode node : nodes) {
			double dist = Math.abs(center[0] - ((Quadrangle) node).v1.lat1) + Math.abs(center[1] - ((Quadrangle) node).v1.long1);
			if (dist < minDist) {
				minDist = dist;
				centerNode = node; 
			}
		}
		List<GraphNode> queue = new ArrayList<GraphNode>();
		filteredTriangles = new ArrayList<>();
		queue.add(centerNode);
		boolean [] done = new boolean[vertices.size()];
		done[centerNode.id] = true;
		filteredTriangles.add(centerNode);
		while (!queue.isEmpty()) {
			GraphNode node = queue.remove(queue.size() - 1);
			done[node.id] = true;
			for (GraphNode neighbour : node.neighbours) {
				if (!done[neighbour.id]) {
					done[neighbour.id] = true;
					queue.add(neighbour);
					filteredTriangles.add(neighbour);
				}
			}
		}
		nodes = filteredTriangles;

		// renumber remaining quadrangles
		renumber();

		// save memory
		for (Vertex v : vertices) {
			v.adjacentGNodes = null;
		}
		
		// adjust distances to map
		double redDistance = reddistanceInput.get();
		for (GraphNode t : nodes) {
			double [] c = t.getCenter();
			int x =(int)( w * (c[1]+180) / 360.0);
			int y =(int)( h * (c[0]+90) / 180.0);
			int color = image.getRGB(x, y) & 0xFFFFFF;
			if (color == 0x00FF00) {
				Arrays.fill(t.distance, redDistance);
			}
		}

		// recalc bounding box
		maxLat = vertices.get(0).lat1;
		minLat = vertices.get(0).lat1;
		maxLong = vertices.get(0).long1;
		maxLong = vertices.get(0).long1;

		for (Vertex v : vertices) {
			maxLat = Math.max(maxLat, v.lat1);
			minLat = Math.min(minLat, v.lat1);
			maxLong = Math.max(maxLong, v.long1);
			minLong = Math.min(minLong, v.long1);
		}
		
		setUpLatLongMap();
		
		// log some stats
		System.err.println("#nodes= " + nodes.size());
	}


	public static void main(String[] args) throws Exception {
		JFrame frame = new JFrame();
		final QuadrangleTesselation tessel = new QuadrangleTesselation();
		//tessel.initByName("depth", 8, "bbox", "10 112 40 154");
		tessel.initByName("depth", 50, "bbox", "5 120 47 154", "reddistance", "0.3");
		
		final List<GraphNode> path1 = new ArrayList<GraphNode>();
		final List<GraphNode> path2 = new ArrayList<GraphNode>();
		final List<GraphNode> path3 = new ArrayList<GraphNode>();
		int n = tessel.nodes.size();
		Randomizer.setSeed(122L);
		final GraphNode t1 = tessel.nodes.get(Randomizer.nextInt(n));
		final GraphNode t2 = tessel.nodes.get(Randomizer.nextInt(n));
		final GraphNode t3 = tessel.nodes.get(Randomizer.nextInt(n));
		tessel.shortestPath(t1, 1.0, t2, 1.0, t3, 1.0, path1, path2, path3);
//		System.err.println("path1 : +" + path1.size() + " path2 : +" + path2.size() + " path3 : +" + path3.size());
		
		final double [] dist = tessel.distances(t1);
//		System.err.println(Arrays.toString(dist));
		
//		DistanceMatrix distances = tessel.distances(4);
//		FileOutputStream fos = new FileOutputStream("/tmp/distances.ser");
//		ObjectOutputStream out = new ObjectOutputStream(fos);
//		out.writeObject(distances);
//		out.close();

		
		
		final BufferedImage image = ImageIO.read(new File("World98b.png"));
		JPanel panel = new JPanel() {
			private static final long serialVersionUID = 1L;
			final boolean scale = true;
			int w;
			double max = - 1;

			
			protected void paintComponent(java.awt.Graphics g) {
				if (max < 0) {			
					for (double d : dist) {
						max = Math.max(max, d);
					}
					System.err.println("max distance = " + max);
				}
				
				this.w = getWidth()/2;
				g.setColor(Color.white);
				g.fillRect(0, 0, getWidth(), getHeight());
				double w = 0, h = 0;
				if (!scale) {
					g.drawImage(image, 0, 0, getWidth(), getHeight(), 0, 0, image.getWidth(), image.getHeight(), null);
					w = getWidth()/360.0;
					h = getHeight()/180.0;
				} else {
					w = getWidth()/(tessel.maxLong - tessel.minLong);
					h = getHeight()/(tessel.maxLat - tessel.minLat);
					g.drawImage(image, 0, 0, getWidth(), getHeight(), 
							(int)(image.getWidth() * (180+tessel.minLong) / 360.0), 
							(int)(image.getHeight() * (90+tessel.minLat) / 180.0), 
							(int)(image.getWidth() * (180+tessel.maxLong) / 360.0),
							(int)(image.getHeight() * (90+tessel.maxLat) / 180.0), null);
				}
				
				
				g.setColor(Color.red);
				g.drawRect((int)((tessel.minLong + 180) * w), (int)((tessel.minLat + 90) * h), (int)((tessel.maxLong -tessel.minLong) * w), (int)((tessel.maxLat - tessel.minLat) * h));
				g.setColor(Color.blue);
				int k = 0;
				for (GraphNode gn : tessel.nodes) {
					Quadrangle t = (Quadrangle) gn;
					if (!scale) {
						drawLine(g, (int)((t.v1.long1 + 180) * w), (int)((t.v1.lat1 + 90) * h), (int)((t.v2.long1 + 180) * w), (int)((t.v2.lat1 + 90) * h));
						drawLine(g, (int)((t.v2.long1 + 180) * w), (int)((t.v2.lat1 + 90) * h), (int)((t.v3.long1 + 180) * w), (int)((t.v3.lat1 + 90) * h));
						drawLine(g, (int)((t.v3.long1 + 180) * w), (int)((t.v3.lat1 + 90) * h), (int)((t.v4.long1 + 180) * w), (int)((t.v4.lat1 + 90) * h));
						drawLine(g, (int)((t.v4.long1 + 180) * w), (int)((t.v4.lat1 + 90) * h), (int)((t.v1.long1 + 180) * w), (int)((t.v1.lat1 + 90) * h));
					} else {
						g.setColor(Color.blue);
						drawLine(g, (int)((t.v1.long1 - tessel.minLong) * w), (int)((t.v1.lat1 - tessel.minLat) * h), (int)((t.v2.long1 - tessel.minLong) * w), (int)((t.v2.lat1 - tessel.minLat) * h));
						drawLine(g, (int)((t.v2.long1 - tessel.minLong) * w), (int)((t.v2.lat1 - tessel.minLat) * h), (int)((t.v3.long1 - tessel.minLong) * w), (int)((t.v3.lat1 - tessel.minLat) * h));
						drawLine(g, (int)((t.v3.long1 - tessel.minLong) * w), (int)((t.v3.lat1 - tessel.minLat) * h), (int)((t.v4.long1 - tessel.minLong) * w), (int)((t.v4.lat1 - tessel.minLat) * h));
						drawLine(g, (int)((t.v4.long1 - tessel.minLong) * w), (int)((t.v4.lat1 - tessel.minLat) * h), (int)((t.v1.long1 - tessel.minLong) * w), (int)((t.v1.lat1 - tessel.minLat) * h));
						
						g.setColor(new Color(Color.HSBtoRGB((float)(dist[k]/max),0.95f, 0.75f)));
						int ow = (int)((t.v3.long1 - tessel.minLong) * w) - (int)((t.v1.long1 - tessel.minLong) * w);
						int oh = (int)((t.v3.lat1 - tessel.minLat) * w) - (int)((t.v1.lat1 - tessel.minLat) * w);
						g.fillOval((int)((t.v1.long1 - tessel.minLong) * w), (int)((t.v1.lat1 - tessel.minLat) * h), ow, oh);
						
					}
					//double [] center = t.getCenter();
					//g.drawString(t.neightbours.length +"", (int)((center[1] - tessel.minLong) * w), (int)((center[0] - tessel.minLat) * h));
					//g.drawString(d[gn.id] +"", (int)((center[1] - tessel.minLong) * w), (int)((center[0] - tessel.minLat) * h));
					k++;
				}
				double [] center = t1.getCenter();
				g.setColor(Color.WHITE);
				int W = 15;
				g.fillOval((int)((center[1] - tessel.minLong) * w)-W/2, (int)((center[0] - tessel.minLat) * h)-W/2, W, W);
				g.setColor(Color.BLACK);
				g.drawOval((int)((center[1] - tessel.minLong) * w)-W/2, (int)((center[0] - tessel.minLat) * h)-W/2, W, W);

				drawPath(g, Color.red, path1, w, h);
				drawPath(g, Color.black, path2, w, h);
				drawPath(g, Color.green, path3, w, h);
			}

			private void drawPath(Graphics g, Color color, List<GraphNode> path1, double w, double h) {
				((Graphics2D)g).setStroke(new BasicStroke(4f));
				g.setColor(color);
				for (int i = 0; i < path1.size() - 1; i++) {
					double [] center1 = path1.get(i).getCenter();
					double [] center2 = path1.get(i+1).getCenter();
					g.drawLine((int)((center1[1] - tessel.minLong) * w), (int)((center1[0] - tessel.minLat) * h),
							(int)((center2[1] - tessel.minLong) * w), (int)((center2[0] - tessel.minLat) * h));
					
				}				
			}

			private void drawLine(Graphics g, int x1, int y1, int x2, int y2) {
				if ((x1 < w && x2 > w) || (x1 > w && x2 < w )) {
					//return;
				}
				
				g.drawLine(x1,y1,x2,y2);
				
			};
		};
		frame.add(panel);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(1024,768);
		frame.setVisible(true);
	}

}
