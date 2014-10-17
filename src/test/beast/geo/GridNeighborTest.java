package test.beast.geo;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import beast.evolution.alignment.distance.GreatCircleDistance;
import beast.util.Randomizer;
import beast.geo.GraphNode;
import beast.geo.GridTesselation;

public class GridNeighborTest {

	public static void main(String[] args) throws Exception {
		test1(false);
		test1(true);
	}
	
	public static void test1(boolean allNeighbors) throws Exception {
		final GridTesselation tessel = new GridTesselation();
		tessel.initByName("depth", 51, "bbox", "5 120 47 154", "reddistance", "1.0", "allNeighbors", allNeighbors);
		
		int n = tessel.nodes.size();
		Randomizer.setSeed(122L);
		final GraphNode t1 = tessel.getLowerLeftCorner();//tessel.nodes.get(Randomizer.nextInt(n));
		//double [] dist = new double [] {0.0};

		
		final double [] dist = tessel.distances(t1);
		System.err.println(Arrays.toString(dist));
		final double [] dist2 = new double[dist.length];
		double [] center = t1.getCenter();
		for (int i = 0; i < dist2.length; i++) {
			GraphNode node = tessel.nodes.get(i);
			double [] ocenter = node.getCenter();
			dist2[node.id] = GreatCircleDistance.pairwiseDistance(center, ocenter);
		}
		System.err.println(Arrays.toString(dist2));

		
		BufferedImage bMap = new BufferedImage(1024, 1024, BufferedImage.TYPE_INT_RGB);
		Graphics g =  bMap.getGraphics();
		g.setClip(new Rectangle(0, 0, 1024, 1024));
		g.setColor(Color.white);
		g.fillRect(0, 0, 1024, 1024);
		g.setColor(Color.black);
		paint(g, tessel, 1024, 1024, dist, dist2);
		
		try {
			if (allNeighbors) {
	 			ImageIO.write(bMap, "png", new File("/tmp/gc_vs_approximate.png"));
				System.err.println("Result written to /tmp/gc_vs_approximate.png");
			} else {
				ImageIO.write(bMap, "png", new File("/tmp/gc_vs_approximate_all_neighbors.png"));
				System.err.println("Result written to /tmp/gc_vs_approximate_all_neighbors.png");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	static 	void paint(Graphics g, GridTesselation tessel, int width, int height, double [] dist, double [] dist2) {
		double w = 0, h = 0;
		g.setColor(Color.white);
		g.fillRect(0, 0, width, height);
		
			//g.drawImage(image, 0, 0, getWidth(), getHeight(), 0, 0, image.getWidth(), image.getHeight(), null);
			w = width/360.0;
			h = height/180.0;
			w = width/(tessel.maxLong - tessel.minLong);
			h = height/(tessel.maxLat - tessel.minLat);
			//g.drawImage(image, 0, 0, getWidth(), getHeight(), 
			//		(int)(image.getWidth() * (180+tessel.minLong) / 360.0), 
			//		(int)(image.getHeight() * (90+tessel.minLat) / 180.0), 
			//		(int)(image.getWidth() * (180+tessel.maxLong) / 360.0),

			double max = -1;
		if (max < 0) {			
//			for (double d : dist) {
//				max = Math.max(max, d);
//			}
			for (double d : dist2) {
				max = Math.max(max, d);
			}
			System.err.println("max distance = " + max);
		}
		
		w = width/2;
//		g.setColor(Color.white);
//		g.fillRect(0, 0, getWidth(), getHeight());
		
		int h1 = height;
		int w1 = width;
		//double w = getWidth()/(tessel.maxLong - tessel.minLong);
		//double h = getHeight()/(tessel.maxLat - tessel.minLat);

		g.setColor(Color.blue);
		double sumsse = 0.0;
		double c = 0.75;//1.0 / 1.33; 
		if (tessel.allNeighborsInput.get()) {
			c = 1/1.06;
		}
		for (int i = 0; i < dist.length; i++) {
			g.drawOval((int)(dist2[i]*w1/max)-1, h1-(int)(dist[i]*c*h1/max)-1, 3, 3);
			double [] center = tessel.nodes.get(i).getCenter();
			//g.drawLine((int)(dist2[i]*w1/max), h1-(int)(dist[i]*c*h1/max), (int)((center[1] - tessel.minLong) * w), (int)((center[0] - tessel.minLat) * h));
			
			sumsse += Math.abs(dist[i]*c- dist2[i]);// * (dist[i]*c- dist2[i]);
		}
		System.err.println(c + " " + (sumsse/dist.length));
		
		g.setColor(Color.black);
		g.drawLine(0,h1,w1,0);
		int k = 20;
		for (int i = 0; i < k; i++) {
			g.drawLine(0,i*h1/k,w1,i*h1/k);
			g.drawLine(i*w1/k,0, i*w1/k,h1);
		}
		
	}


}
