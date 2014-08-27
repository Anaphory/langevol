package beast.evolution.geo;

import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("Set of lines on a mercator projection, upper left corner = 0,0 lower right corned = w,h")
public class LineSet extends BEASTObject {
	public Input<Double> wInput = new Input<Double>("w", "width of the drawing area, corresponds to 360 degrees", Validate.REQUIRED);
	public Input<Double> hInput = new Input<Double>("h", "height of the drawing area, correspondes to 180 degrees", Validate.REQUIRED);
	public Input<List<Line>> linesInput = new Input<List<Line>>("line", "set of lines on (0,0)-(w,h) plane", new ArrayList<Line>());

	double w;
	double h;
	
	double [][] latitudes;
	double [][] longitudes;
	
	@Override
	public void initAndValidate() throws Exception {
		w = wInput.get();
		h = hInput.get();
		
		double halfw = w / 2.0;
		double halfh = h / 2.0;
		
		List<Line> lines = linesInput.get();
		latitudes = new double [lines.size()][2];
		longitudes = new double [lines.size()][2];
		for (int i = 0; i < lines.size(); i++) {
			double x1 = lines.get(i).xInput.get();
			double y1 = lines.get(i).yInput.get();
			double x2 = x1 + lines.get(i).wInput.get();
			double y2 = y1 + lines.get(i).hInput.get();
			
			
			longitudes[i][0] = 360.0 * (x1-halfw)/w; 
			longitudes[i][1] = 360.0 * (x2-halfw)/w; 
			latitudes[i][0] = 180.0 * (halfh - y1)/h; 
			latitudes[i][1] = 180.0 * (halfh - y2)/h; 
		}
	}
	
	public double getLatStart(int i)  {return latitudes[i][0];}
	public double getLatStop(int i)   {return latitudes[i][1];}
	public double getLongStart(int i) {return longitudes[i][0];}
	public double getLongStop(int i)  {return longitudes[i][1];}

	public double[][] getLatitudes() { return latitudes;}
	public double[][] getLongitudes() {return longitudes;}

}
