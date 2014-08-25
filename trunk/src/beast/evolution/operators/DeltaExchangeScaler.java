package beast.evolution.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

@Description("Simultaneously delta exchange and scale")
public class DeltaExchangeScaler extends Operator {
    public final Input<Double> deltaInput = new Input<Double>("delta", "Magnitude of change for two randomly picked values.", 1.0);

    public Input<RealParameter> upInput = new Input<RealParameter>("deltaexchange",
            "delta exchange is applied to this state node", Validate.REQUIRED);
    public Input<RealParameter> downInput = new Input<RealParameter>("scale",
            "scale operation is applied to this state node", Validate.REQUIRED);

    RealParameter realparameter;
    RealParameter scaleNode;
    double delta = 1;
    
    @Override
	public void initAndValidate() throws Exception {
    	realparameter = upInput.get();
    	scaleNode = downInput.get();
    	delta = deltaInput.get();
	}

	@Override
	public double proposal() {
        final int dim = realparameter.getDimension();

        final int dim1 = Randomizer.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = Randomizer.nextInt(dim);
        }

        // operate on real parameter
        double scalar1 = realparameter.getValue(dim1);
        double scalar2 = realparameter.getValue(dim2);

        // exchange a random delta
        final double d = Randomizer.nextDouble() * delta;
        scalar1 -= d;
        scalar2 += d;

        if (scalar1 < realparameter.getLower() || scalar1 > realparameter.getUpper() ||
                scalar2 < realparameter.getLower() || scalar2 > realparameter.getUpper()) {
            return Double.NEGATIVE_INFINITY;
        }

        double value = scaleNode.getValue();
        double newValue = (dim1 < dim2 ? value * Math.exp(d) : value * Math.exp(-d));
        if (newValue < scaleNode.getLower() || newValue > scaleNode.getUpper()) {
        	return Double.NEGATIVE_INFINITY;
        }
        
        realparameter.setValue(dim1, scalar1);
        realparameter.setValue(dim2, scalar2);
        scaleNode.setValue(newValue);
        
        return 0;
	}

}
