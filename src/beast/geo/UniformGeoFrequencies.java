package beast.geo;


import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.substitutionmodel.Frequencies;

@Description("Uniform distribution over graph nodes in a graph inside a region specified by a KML file")
public class UniformGeoFrequencies extends Frequencies {
	public Input<Graph> graphInput = new Input<Graph>("graph", "graph containing all locations", Validate.REQUIRED);
	public Input<File> kmlFileInput = new Input<File>("kml", "kml file with polygons over admissable locations. If not specified, all locations are admissable.");

	

	Graph graph;
	IntegerParameter location;

	int taxonNr;
	boolean[] isAdmissable;

	boolean onroot = false;


	public UniformGeoFrequencies() {
		frequenciesInput.setRule(Validate.OPTIONAL);
		dataInput.setRule(Validate.OPTIONAL);
		
	}
	@Override
	public void initAndValidate() {
		graph = graphInput.get();

		if (kmlFileInput.get()==null) {
			// every graph node is admissable
			isAdmissable = new boolean[graph.getSize()];
			Arrays.fill(isAdmissable, true);
		} else {
			
				List<List<Double>> coordinates;
				try {
					coordinates = parseKML();
				} catch (SAXException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					throw new RuntimeException(e.getMessage());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					throw new RuntimeException(e.getMessage());
				} catch (ParserConfigurationException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					throw new RuntimeException(e.getMessage());
				}
				try {
					calcAdmissableNodes(coordinates);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					throw new RuntimeException(e.getMessage());
				}
		
		}
		calcFreqs();
	}

	
	private void calcFreqs() {
		int n = isAdmissable.length;
		freqs = new double[n];
		int sum = 0;
		for (int i = 0; i < n; i++) {
			if (isAdmissable[i]) {
				sum++;
			}
		}
		
		for (int i = 0; i < n; i++) {
			if (isAdmissable[i]) {
				freqs[i] = 1.0/sum;
			}
		}
		
		
	}
	private void calcAdmissableNodes(List<List<Double>> coordinates) throws IOException {
		boolean debug = true;//Boolean.valueOf(System.getProperty("beast.debug"));
		
		isAdmissable = new boolean[graph.getSize()];
		
		double minLong = 360, maxLong = -360, maxLat = -90, minLat = 180;
		for (List<Double> coords : coordinates) {
			for (int i = 0; i < coords.size(); i += 2) {
				double latitude = coords.get(i);				
				double longitude = coords.get(i+1);
				minLong = Math.min(minLong, longitude);
				maxLong = Math.max(maxLong, longitude);
				minLat = Math.min(minLat, latitude);
				maxLat = Math.max(maxLat, latitude);
			}
		}
		
		
		// create bitmap to draw in
		int width = 1024;
		int height = 1024;
		double w = width/(maxLong - minLong);
		double h = height/(maxLat - minLat);
		
		// draw polygons in an image
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics g = image.getGraphics();

		if (debug) {
			final BufferedImage worldimage = ImageIO.read(new File("World98b.png"));
			System.err.println((int)(worldimage.getWidth() * (180+minLong) / 360.0));
					System.err.println((int)(worldimage.getHeight() * (90+minLat) / 180.0)); 
							System.err.println((int)(worldimage.getWidth() * (180+maxLong) / 360.0));
			System.err.println((int)(worldimage.getHeight() * (90+maxLat) / 180.0));

			
			g.drawImage(worldimage, 0, 0, width, height, 
					(int)(worldimage.getWidth() * (180+minLong) / 360.0), 
					(int)(worldimage.getHeight() * (minLat) / 180.0), 
					(int)(worldimage.getWidth() * (180+maxLong) / 360.0),
					(int)(worldimage.getHeight() * (maxLat) / 180.0), null);
			ImageIO.write(image, "png", new File("/tmp/kmlrange" + getID( )+".png"));
		}
		
		for (List<Double> coords : coordinates) {
			int nPoints = coords.size()/2;
			int [] xPoints = new int[nPoints];
			int [] yPoints = new int[nPoints];
			for (int i = 0; i < coords.size(); i += 2) {
				xPoints[i/2] = (int)((coords.get(i+1) - minLong) * w);
				yPoints[i/2] = (int)((coords.get(i) - minLat) * h);				
			}
			g.setColor(Color.blue);
			g.fillPolygon(xPoints, yPoints, nPoints);
		}		
		
		
		// determine whether a graph node is admissable
		for (GraphNode node : graph.nodes){
			double [] center = node.getCenter();
			int x = (int)((center[1] - minLong) * w);
			int y = (int)((center[0] + 90 - minLat) * h);
			if (x>=0 && x < width && y >=0 && y < height) {
				int color = image.getRGB(x, y) & 0xFFFFFF;
				if (color == 0x0000FF) {
					isAdmissable[node.id] = true;
					if (debug) {
						g.setColor(Color.red);
						g.fillOval(x-5, y-5, 11, 11);
					}
				} else {
					if (debug) {
						g.setColor(Color.black);
						g.fillOval(x-3, y-3, 7, 7);
					}
				}
			}
		}
		
		if (debug) {
			ImageIO.write(image, "png", new File("/tmp/kmlprior" + getID() +".png"));
		}
		
	}


	private List<List<Double>> parseKML() throws SAXException, IOException, ParserConfigurationException {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		factory.setValidating(false);
		org.w3c.dom.Document doc = factory.newDocumentBuilder().parse(kmlFileInput.get());
		doc.normalize();

		List<List<Double>> coordinates = new ArrayList<List<Double>>();

		// grab 'coordinates' elements out of the KML file
		NodeList oCoordinates = doc.getElementsByTagName("coordinates");
		for (int iNode = 0; iNode < oCoordinates.getLength(); iNode++) {
			Node oCoordinate = oCoordinates.item(iNode);
			String sCoordinates = oCoordinate.getTextContent();
			List<Double> polygon = new ArrayList<>();
			String[] sStrs = sCoordinates.split("\\s+");
			for (String sStr : sStrs) {
				if (sStr.contains(",")) {
					String[] sCoords = sStr.split(",");
					polygon.add(90 - Double.parseDouble(sCoords[1].trim()));
					polygon.add(Double.parseDouble(sCoords[0].trim()));
				}
			}
			coordinates.add(polygon);
		}
		return coordinates;
	}
	
    @Override
    protected boolean requiresRecalculation() {
    	return false;
    }
    
    @Override
    public double[] getFreqs() {
    	return freqs;
    }
    
	public static void main(String[] args) throws Exception {
		System.setProperty("beast.debug", "true");
		
		Taxon taxon = new Taxon();
		taxon.setID("corvus");
		
		TaxonSet taxonSet = new TaxonSet();
		taxonSet.initByName("taxon", taxon);

		final QuadrangleTesselation tessel = new QuadrangleTesselation();
		tessel.initByName("depth", 51, "bbox", "-90 -180 89 179", 
				"reddistance", "1.0", "allNeighbors", false, "removeWater", false);
		
		UniformGeoFrequencies distr = new UniformGeoFrequencies();
		distr.initByName(//"taxon", taxon, "taxonset", taxonSet, 
				"graph", tessel,
				"kml", "/home/remco/data/beast/corvus/geo/kml/Dendrocitta_vagabunda_22705836.kml");
				//"location", new IntegerParameter());
	}

}
