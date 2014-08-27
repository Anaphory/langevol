package beast.evolution.geo;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.ClusterTree;
import beast.util.Randomizer;

public class ClusterTraitTree extends ClusterTree {
	public Input<RealParameter> locationInput = new Input<RealParameter>("trait", "location trait where last two elements represent latitude, longitude", Validate.REQUIRED);

	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		
		RealParameter parameter = locationInput.get();

		Double [] values = new Double[parameter.getDimension()];
		for (int i = 0; i < parameter.getDimension(); i++) {
			values[i] = parameter.getValue(i);
		}
		int dim = parameter.getMinorDimension1();
		initInternalNodes(getRoot(), values, dim);

		
		RealParameter tmp = new RealParameter(values);	        
		tmp.setBounds(parameter.getLower(), parameter.getUpper());
		parameter.assignFromWithoutID(tmp);
	}

	/** set trait value as mean of its children **/
	private void initInternalNodes(Node node, Double[] values, int dim) {
		double jitter = 0.0;// jitterInput.get();
		if (!node.isLeaf()) {
			for (Node child : node.getChildren()) {
				initInternalNodes(child, values, dim);
			}
			for (int i = 0; i < dim; i++) {
				double value = 0;
				for (Node child : node.getChildren()) {
					value += values[child.getNr() * dim + i];
				}
				values[node.getNr() * dim + i] = value / dim + Randomizer.nextDouble() * jitter - jitter / 2.0;
			}
		}
	
	}

}
