package beast.geo;

import com.sun.org.glassfish.gmbal.Description;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.ContinuousSubstitutionModel;


@Description("Approximates diffusion model based on min distance between two points")
public class DistanceDiffusionModel extends ContinuousSubstitutionModel {
	public Input<Graph> graphInput = new Input<Graph>("graph","graph representing distances between points", Validate.REQUIRED);
	public Input<RealParameter> precisionInput = new Input<RealParameter>("precision", "precision governing diffusion process", Validate.REQUIRED);
	
	Graph graph;
	RealParameter precision;
	
	@Override
	public void initAndValidate() throws Exception {
		graph = graphInput.get();
		precision = precisionInput.get();
		super.initAndValidate();
	}
	
	
	@Override
	public double getLogLikelihood(double[] start, double[] stop, double time) {
		double distance = graph.getDistance(start, stop);

		double inverseVariance = precision.getValue(0) / time;
        double logP = Math.log(distance) + 0.5 * Math.log(inverseVariance) -0.5 * distance * distance * inverseVariance;
        return logP;
	}

}
