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
import beast.evolution.alignment.distance.GreatCircleDistance;
import beast.util.Randomizer;

@Description("Tesselates (part of) a sphere with equal sized quadrangle on a grid")
public class GridTesselation extends SphereTesselation {

		public Input<File> mapFileInput = new Input<File>("map","geographic world map in Mercator projection -- red coloured items will have distance defined by 'reddistance'. "
				+ "If not specified, all distances are 1.0",
				new File("/home/remco/data/geo/aboriginal25.bmp"));
		public Input<Double> reddistanceInput = new Input<Double>("reddistance","distance fpr red coloured items", 1.0);
		
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
			
			GraphNode [][] grid = new GraphNode[depth-1][depth-1];
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
					grid[i][j].setUpDistances();
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
	
	
	
	

}
