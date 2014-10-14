package beast.geo;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;

@Description("Tesselates a sphere with equal sized triangles")
public class SphereTesselation extends Graph {
	public Input<Integer> depthInput = new Input<Integer>("depth", "number of splits of base tesselation - higher means more traingles (4^depth)", 3);

	public Input<String> bboxInput = new Input<String>("bbox", "bounding box as space separated list (in min-latitude min-longitude max-latitude max-longitude)"
			+ " e.g. ");


	double maxLat = 90, minLat = -90;
	double maxLong= 180, minLong= -180;
	
	@Override
	public void initAndValidate() throws Exception {
		nodes = new ArrayList<GraphNode>();
		double phi = Math.atan(0.5) * 180 / Math.PI;
		
		// north pole
		Vertex np = new Vertex(-90,0);
		Vertex sp = new Vertex(90,0);
		
		Vertex n1 = new Vertex(-phi, -180);
		Vertex n2 = new Vertex(-phi, -108);
		Vertex n3 = new Vertex(-phi,  -36);
		Vertex n4 = new Vertex(-phi,   36);
		Vertex n5 = new Vertex(-phi,  108);
		
		Vertex s1 = new Vertex(phi, -144);
		Vertex s2 = new Vertex(phi,  -72);
		Vertex s3 = new Vertex(phi,    0);
		Vertex s4 = new Vertex(phi,   72);
		Vertex s5 = new Vertex(phi,  144);
		
		nodes.add(new Triangle(np, n1, n2));
		nodes.add(new Triangle(np, n2, n3));
		nodes.add(new Triangle(np, n3, n4));
		nodes.add(new Triangle(np, n4, n5));
		nodes.add(new Triangle(np, n5, n1));

		nodes.add(new Triangle(s1, n1, n2));
		nodes.add(new Triangle(s2, n2, n3));
		nodes.add(new Triangle(s3, n3, n4));
		nodes.add(new Triangle(s4, n4, n5));
		nodes.add(new Triangle(s5, n5, n1));

		nodes.add(new Triangle(sp, s1, s2));
		nodes.add(new Triangle(sp, s2, s3));
		nodes.add(new Triangle(sp, s3, s4));
		nodes.add(new Triangle(sp, s4, s5));
		nodes.add(new Triangle(sp, s5, s1));

		nodes.add(new Triangle(n2, s1, s2));
		nodes.add(new Triangle(n3, s2, s3));
		nodes.add(new Triangle(n4, s3, s4));
		nodes.add(new Triangle(n5, s4, s5));
		nodes.add(new Triangle(n1, s5, s1));
		
		for (int i = 0; i < depthInput.get(); i++) {
			List<GraphNode> newTriangles = new ArrayList<>();
			for (GraphNode t : nodes) {
				((Triangle)t).split(newTriangles);
			}
			nodes = newTriangles;
		}

		
		System.err.println("#triangels = " + nodes.size() + " before filtering");
		// filtering out traingle outside bounding box (if any)
		if (bboxInput.get() != null) {
			String str = bboxInput.get().trim();
			String [] strs = str.split("\\s+");
			if (strs.length != 4) {
				throw new RuntimeException("bbox input must contain 4 numbers");
			}
			minLat = Double.parseDouble(strs[0]);
			minLong = Double.parseDouble(strs[1]);
			maxLat = Double.parseDouble(strs[2]);
			maxLong = Double.parseDouble(strs[3]);
			if (minLat >= maxLat) {
				throw new RuntimeException("bbox input first latitude must be smaller than second latitude");
			}
			if (minLong >= maxLong) {
				throw new RuntimeException("bbox input first longitude must be smaller than second longitude");
			}
			
			List<GraphNode> filteredTriangles = new ArrayList<>();
			for (GraphNode t : nodes) {
				if (((Triangle)t).hasPointsInside(minLat, minLong, maxLat, maxLong)) {
					filteredTriangles.add(t);
				}
			}
			nodes = filteredTriangles;
		}
		
		// collect vertices
		List<Vertex> vertices = new ArrayList<Vertex>();
		for (GraphNode t : nodes) {
			vertices.add(((Triangle)t).v1);
			vertices.add(((Triangle)t).v2);
			vertices.add(((Triangle)t).v3);
		}
		
		// sort vertices by latitude/longitude in order to find duplicates
		class VertexComparator implements Comparator<Vertex> {
			final double EPSILON = 1e-8;
			@Override
			public int compare(Vertex v1, Vertex v2) {
				if (Math.abs(v1.lat1 - v2.lat1) > EPSILON) {
					if (v1.lat1 > v2.lat1) {
						return 1;
					} else {
						return -1;
					}
				}
				if (Math.abs(v1.long1 - v2.long1) > EPSILON) {
					if (v1.long1 > v2.long1) {
						return 1;
					} else {
						return -1;
					}
				}
				return 0;
			}
		};
		
		VertexComparator comparator = new VertexComparator();
		vertices.sort(comparator);
		
		// join adjacentTraingles of duplicate vertices
		for (int i = 0; i < vertices.size() - 1; i++) {
			Vertex v1 = vertices.get(i);
			Vertex v2 = vertices.get(i+1);
			if (comparator.compare(v1, v2) == 0) { 
				v1.adjacentTraingles.addAll(v2.adjacentTraingles);
				v2.adjacentTraingles.addAll(v1.adjacentTraingles);
			}
		}
		
		// renumber remaining triangles
		int i = 0;
		for (GraphNode t : nodes) {
			t.id = i++;
		}

		// set up adjacency graph -- requires vertices to have adjacentTraingles to be set up
		for (GraphNode t : nodes) {
			((Triangle)t).calcNeighbours();
		}
		
		// save memory
		for (Vertex v : vertices) {
			v.adjacentTraingles = null;
		}

		
		// log some stats
		int k = 0;
		for (GraphNode t : nodes) {
			if (t.neighbours.length != 3) {
				k++;
			}
		}
		System.err.println("#triangels = " + nodes.size());	
		System.err.println("#triangels with less than 3 neighbors = " + k);	
		
	}


	public static void main(String[] args) throws Exception {
		JFrame frame = new JFrame();
		final SphereTesselation tessel = new SphereTesselation();
//		tessel.initByName("depth", 8, "bbox", "10 112 40 154");
		tessel.initByName("depth",4);
		
		final List<GraphNode> path1 = new ArrayList<GraphNode>();
		final List<GraphNode> path2 = new ArrayList<GraphNode>();
		final List<GraphNode> path3 = new ArrayList<GraphNode>();
		int n = tessel.nodes.size();
		Randomizer.setSeed(129L);
		final GraphNode t1 = tessel.nodes.get(Randomizer.nextInt(n));
		final GraphNode t2 = tessel.nodes.get(Randomizer.nextInt(n));
		final GraphNode t3 = tessel.nodes.get(Randomizer.nextInt(n));
		tessel.shortestPath(t1, 1.0, t2, 1.0, t3, 1.0, path1, path2, path3);
		System.err.println("path1 : +" + path1.size() + " path2 : +" + path2.size() + " path3 : +" + path3.size());
		
		final double [] d = tessel.distances(t1);
		System.err.println(Arrays.toString(d));
		
		final BufferedImage image = ImageIO.read(new File("World98b.png"));
		JPanel panel = new JPanel() {
			private static final long serialVersionUID = 1L;
			int w;

			protected void paintComponent(java.awt.Graphics g) {
				this.w = getWidth()/2;
				g.setColor(Color.white);
				g.fillRect(0, 0, getWidth(), getHeight());
//				g.drawImage(image, 0, 0, getWidth(), getHeight(), 0, 0, image.getWidth(), image.getHeight(), null);
//				double w = getWidth()/360.0;
//				double h = getHeight()/180.0;
				double w = getWidth()/(tessel.maxLong - tessel.minLong);
				double h = getHeight()/(tessel.maxLat - tessel.minLat);
				g.drawImage(image, 0, 0, getWidth(), getHeight(), 
						(int)(image.getWidth() * (180+tessel.minLong) / 360.0), 
						(int)(image.getHeight() * (90+tessel.minLat) / 180.0), 
						(int)(image.getWidth() * (180+tessel.maxLong) / 360.0),
						(int)(image.getHeight() * (90+tessel.maxLat) / 180.0), null);
				
				
				g.setColor(Color.red);
				g.drawRect((int)((tessel.minLong + 180) * w), (int)((tessel.minLat + 90) * h), (int)((tessel.maxLong -tessel.minLong) * w), (int)((tessel.maxLat - tessel.minLat) * h));
				g.setColor(Color.blue);
				for (GraphNode gn : tessel.nodes) {
					Triangle t = (Triangle) gn;
//					drawLine(g, (int)((t.v1.long1 + 180) * w), (int)((t.v1.lat1 + 90) * h), (int)((t.v2.long1 + 180) * w), (int)((t.v2.lat1 + 90) * h));
//					drawLine(g, (int)((t.v1.long1 + 180) * w), (int)((t.v1.lat1 + 90) * h), (int)((t.v3.long1 + 180) * w), (int)((t.v3.lat1 + 90) * h));
//					drawLine(g, (int)((t.v2.long1 + 180) * w), (int)((t.v2.lat1 + 90) * h), (int)((t.v3.long1 + 180) * w), (int)((t.v3.lat1 + 90) * h));
					drawLine(g, (int)((t.v1.long1 - tessel.minLong) * w), (int)((t.v1.lat1 - tessel.minLat) * h), (int)((t.v2.long1 - tessel.minLong) * w), (int)((t.v2.lat1 - tessel.minLat) * h));
					drawLine(g, (int)((t.v1.long1 - tessel.minLong) * w), (int)((t.v1.lat1 - tessel.minLat) * h), (int)((t.v3.long1 - tessel.minLong) * w), (int)((t.v3.lat1 - tessel.minLat) * h));
					drawLine(g, (int)((t.v2.long1 - tessel.minLong) * w), (int)((t.v2.lat1 - tessel.minLat) * h), (int)((t.v3.long1 - tessel.minLong) * w), (int)((t.v3.lat1 - tessel.minLat) * h));
					double [] center = t.getCenter();
					//g.drawString(t.neightbours.length +"", (int)((center[1] - tessel.minLong) * w), (int)((center[0] - tessel.minLat) * h));
					g.drawString(d[gn.id] +"", (int)((center[1] - tessel.minLong) * w), (int)((center[0] - tessel.minLat) * h));
				}
//				drawPath(g, Color.red, path1, w, h);
//				drawPath(g, Color.black, path2, w, h);
//				drawPath(g, Color.green, path3, w, h);
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
					return;
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
