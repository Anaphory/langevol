package beast.geo;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;

@Description("Tesselates a sphere with equal sized triangles")
public class SphereTesselation extends BEASTObject {
	public Input<Integer> depthInput = new Input<Integer>("depth", "number of splits of base tesselation - higher means more traingles (4^depth)", 3);

	public Input<String> bboxInput = new Input<String>("bbox", "bounding box as space separated list (in min-latitude min-longitude max-latitude max-longitude)"
			+ " e.g. ");

	List<Triangle> triangles;

	double maxLat = 90, minLat = -90;
	double maxLong= 180, minLong= -180;
	
	@Override
	public void initAndValidate() throws Exception {
		triangles = new ArrayList<Triangle>();
		double phi = Math.atan(0.5) * 180 / Math.PI;
		for (int i = -180; i < 180; i+=72) {
			// to the north pole
			Triangle t = new Triangle(new Vertex(90, i + 36), new Vertex(phi, i), new Vertex(phi, i + 72));
			triangles.add(t);
			// to the south pole
			t = new Triangle(new Vertex(-90,i + 72), new Vertex(-phi, i + 36), new Vertex(-phi, i + 108));
			triangles.add(t);
			// on the equator
			t = new Triangle(new Vertex(phi, i + 72), new Vertex(-phi, i + 36), new Vertex(-phi, i + 108));
			triangles.add(t);
			t = new Triangle(new Vertex(-phi, i + 36), new Vertex(phi, i), new Vertex(phi, i + 72.0));
			triangles.add(t);
		}
		for (int i = 0; i < depthInput.get(); i++) {
			List<Triangle> newTriangles = new ArrayList<>();
			for (Triangle t : triangles) {
				t.split(newTriangles);
			}
			triangles = newTriangles;
		}

		
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
			
			List<Triangle> filteredTriangles = new ArrayList<>();
			for (Triangle t : triangles) {
				if (t.hasPointsInside(minLat, minLong, maxLat, maxLong)) {
					filteredTriangles.add(t);
				}
			}
			triangles = filteredTriangles;
		}

	}
	
	
	public static void main(String[] args) throws Exception {
		JFrame frame = new JFrame();
		final SphereTesselation tessel = new SphereTesselation();
		tessel.initByName("depth", 5, "bbox", "10 112 40 154");
		final BufferedImage image = ImageIO.read(new File("World98b.png"));
		JPanel panel = new JPanel() {
			private static final long serialVersionUID = 1L;
			int w;

			protected void paintComponent(java.awt.Graphics g) {
				w = getWidth()/2;
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
				for (Triangle t : tessel.triangles) {
//					drawLine(g, (int)((t.long1 + 180) * w), (int)((t.lat1 + 90) * h), (int)((t.long2 + 180) * w), (int)((t.lat2 + 90) * h));
//					drawLine(g, (int)((t.long1 + 180) * w), (int)((t.lat1 + 90) * h), (int)((t.long3 + 180) * w), (int)((t.lat3 + 90) * h));
//					drawLine(g, (int)((t.long2 + 180) * w), (int)((t.lat2 + 90) * h), (int)((t.long3 + 180) * w), (int)((t.lat3 + 90) * h));
					drawLine(g, (int)((t.v1.long1 - tessel.minLong) * w), (int)((t.v1.lat1 - tessel.minLat) * h), (int)((t.v2.long1 - tessel.minLong) * w), (int)((t.v2.lat1 - tessel.minLat) * h));
					drawLine(g, (int)((t.v1.long1 - tessel.minLong) * w), (int)((t.v1.lat1 - tessel.minLat) * h), (int)((t.v3.long1 - tessel.minLong) * w), (int)((t.v3.lat1 - tessel.minLat) * h));
					drawLine(g, (int)((t.v2.long1 - tessel.minLong) * w), (int)((t.v2.lat1 - tessel.minLat) * h), (int)((t.v3.long1 - tessel.minLong) * w), (int)((t.v3.lat1 - tessel.minLat) * h));
				}
			}

			private void drawLine(Graphics g, int x1, int y1, int x2, int y2) {
				if (x1 < w && x2 > w) {
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
