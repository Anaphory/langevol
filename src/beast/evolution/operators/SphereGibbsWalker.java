package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.ContinuousSubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

public class SphereGibbsWalker extends Operator {
	public Input<RealParameter> locationInput = new Input<RealParameter>("location", "latitude/longitude pairs representing location", Validate.REQUIRED);
	public Input<TreeInterface> treeInput = new Input<>("tree", "tree for whcih to infer geography", Validate.REQUIRED);
    //public Input<Operator> operatorInput = new Input<Operator>("operator" ,"optional tree operator -- locations of filthy nodes will get a new locaiton");
    public Input<ContinuousSubstitutionModel> modelInput = new Input<>("model" ,"substitution model that determines the distribution over locations", Validate.REQUIRED);

    public Input<Double> minLatInput = new Input<Double>("minLat","minimum latitude for grid", -90.0);
    public Input<Double> maxLatInput = new Input<Double>("maxLat","maximum latitude for grid", 90.0);
    public Input<Double> latStepInput = new Input<Double>("latStep","stepsize in latitude direction for grid", 1.0);
    public Input<Double> minLongInput = new Input<Double>("minLong","minimum longitude for grid", -180.0);
    public Input<Double> maxLongInput = new Input<Double>("maxLong","maximum longitude for grid", 180.0);
    public Input<Double> longStepInput = new Input<Double>("longStep","stepsize in longitude direction for grid", 1.0);
    
    
    double minLat, maxLat, minLong, maxLong, latStep, longStep;
	int nY;
	int nX;
    double [] logP;
    
    ContinuousSubstitutionModel model;
    TreeInterface tree;
    RealParameter location;
    
    @Override
	public void initAndValidate() throws Exception {
    	minLat = minLatInput.get();
    	maxLat = maxLatInput.get();
    	latStep = latStepInput.get();
    	minLong = minLongInput.get();
    	maxLong = maxLongInput.get();
    	longStep = longStepInput.get();
    	nY = (int)((maxLat - minLat)/latStep);
    	nX = (int)((maxLong - minLong)/longStep);
    	logP = new double[nX * nY];
    	
    	model = modelInput.get();
    	tree = treeInput.get();
    	location = locationInput.get();
	}

	@Override
	public double proposal() {
		int iNode = tree.getLeafNodeCount() + Randomizer.nextInt(tree.getInternalNodeCount());
		Node node = tree.getNode(iNode);
		int parent = node.getParent().getNr();
		int child1 = node.getLeft().getNr();
		int child2 = node.getRight().getNr();

		double timeNode = node.getLength();
		double timeChild1 = node.getLeft().getLength();
		double timeChild2 = node.getRight().getLength();
		double [] parentLoc = new double[2];
		double [] child1Loc = new double[2];
		double [] child2Loc = new double[2];
		double [] nodeLoc = new double[2];
		parentLoc[0] = location.getValue(parent * 2);
		parentLoc[1] = location.getValue(parent * 2 + 1);
		child1Loc[0] = location.getValue(child1 * 2);
		child1Loc[1] = location.getValue(child1 * 2 + 1);
		child2Loc[0] = location.getValue(child2 * 2);
		child2Loc[1] = location.getValue(child2 * 2 + 1);
		
		int k = 0;
		for (double lat = minLat; lat < maxLat; lat += latStep) {
			nodeLoc[0] = lat;
			for (double _long = minLong; _long < maxLong; _long += longStep) {
				nodeLoc[1] = _long;
				logP[k] = model.getLogLikelihood(parentLoc, nodeLoc, timeNode) +
						model.getLogLikelihood(nodeLoc, child1Loc, timeChild1) +
						model.getLogLikelihood(nodeLoc, child2Loc, timeChild2);
				k++;
			}
		}
		double max = 0;
		for (double d : logP) {
			max = Math.max(max, d);
		}
		for (int i = 0; i < logP.length; i++) {
			logP[i] = Math.exp(logP[i]-max);
		}
		k = Randomizer.randomChoicePDF(logP);
		double newLat = minLat + (k / nX) * latStep;
		double newLong = minLong + (k % nX) * longStep;
		location.setValue(iNode * 2, newLat);
		location.setValue(iNode * 2 + 1, newLong);
		return 0;
	}

}
