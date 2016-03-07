package beast.evolution.operators;

import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.ContinuousSubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

@Description("Gibbs like operator to walk on a sphere")
public class SphereGibbsWalker extends MoveAboutMeanOperator {
    //public Input<Operator> operatorInput = new Input<Operator>("operator" ,"optional tree operator -- locations of filthy nodes will get a new locaiton");
    public Input<ContinuousSubstitutionModel> modelInput = new Input<>("model" ,"substitution model that determines the distribution over locations", Validate.REQUIRED);

    public Input<Double> minLatInput = new Input<Double>("minLat","minimum latitude for grid", -90.0);
    public Input<Double> maxLatInput = new Input<Double>("maxLat","maximum latitude for grid", 90.0);
    public Input<Double> latStepInput = new Input<Double>("latStep","stepsize in latitude direction for grid", 2.0);
    public Input<Double> minLongInput = new Input<Double>("minLong","minimum longitude for grid", -180.0);
    public Input<Double> maxLongInput = new Input<Double>("maxLong","maximum longitude for grid", 180.0);
    //public Input<Double> longStepInput = new Input<Double>("longStep","stepsize in longitude direction for grid", 2.0);

    
    double minLat, maxLat, minLong, maxLong, latStep, longStep;
	int nY;
	int nX;
    double [] logP;
    
    ContinuousSubstitutionModel model;

    
    @Override
	public void initAndValidate() {
    	minLat = minLatInput.get();
    	maxLat = maxLatInput.get();
    	latStep = latStepInput.get();
    	minLong = minLongInput.get();
    	maxLong = maxLongInput.get();
    	// ensure there are as many steps in long direcation as in lat direction
    	longStep = latStep * (maxLong - minLong) / (maxLat - minLat);
    	nY = (int)((maxLat - minLat)/latStep);
    	nX = (int)((maxLong - minLong)/longStep);
    	logP = new double[nX * nY];
    	
    	model = modelInput.get();
    	tree = treeInput.get();
    	location = locationInput.get();
    	
    	operator = operatorInput.get();
	}

	@Override
	public double proposal() {
		if (operator != null) {
			// do operator, then do random walk on all filthy nodes
			double logHR = operator.proposal();
			if (Randomizer.nextBoolean()) {
				for (Node node : tree.getNodesAsArray()) {
					if (node.isDirty() == Tree.IS_FILTHY) {
	                    // Why not isLeaf?? (JH) (Q2R)
						if (node.getNr() >= tree.getLeafNodeCount()) {
							logHR += doproposal(node.getNr());
						}
					}
				}
			} else {
				logHR = traverse(tree.getRoot(), new boolean[1]);
			}
			return logHR;
		}
		
		final int iNode = tree.getLeafNodeCount() + Randomizer.nextInt(tree.getInternalNodeCount());
		if (Randomizer.nextBoolean()) {
	        Node node = tree.getNode(iNode);
	        while (!node.isRoot()) {
	        	doproposal(node.getNr());
	        	node = node.getParent();
	        }
		} else {
        	doproposal(iNode);
		}
        return 0;
	}
    

	double [] child1Loc = new double[2];
	double [] child2Loc = new double[2];
	double [] nodeLoc = new double[2];
	double [] parentLoc = new double[2];
	double timeChild1;
	double timeChild2;
	double timeNode;
	double offset0 ;
	double offset1;
    final static double THRESHOLD = -100;
    
	public double doproposal(int iNode) {
		//int iNode = tree.getLeafNodeCount() + Randomizer.nextInt(tree.getInternalNodeCount());
		Node node = tree.getNode(iNode);

		
		int child1 = node.getLeft().getNr();
		int child2 = node.getRight().getNr();
		timeChild1 = node.getLeft().getLength();
		timeChild2 = node.getRight().getLength();
		child1Loc[0] = location.getValue(child1 * 2);
		child1Loc[1] = location.getValue(child1 * 2 + 1);
		child2Loc[0] = location.getValue(child2 * 2);
		child2Loc[1] = location.getValue(child2 * 2 + 1);
		
		offset0 = location.getValue(iNode * 2);
		offset0 = offset0 - Math.round(offset0 / latStep) * latStep;
		offset1 = location.getValue(iNode  * 2 + 1);
		offset1 = offset1 - Math.round(offset1/ longStep) * longStep;
		double max = Double.NEGATIVE_INFINITY;

		Arrays.fill(logP, 0.0);
		
		if (node.isRoot()) {
			if (true) {return Double.NEGATIVE_INFINITY;}
			int k = 0;
			for (double lat = minLat; lat < maxLat; lat += latStep) {
				nodeLoc[0] = lat + offset0;
				for (double _long = minLong; _long < maxLong; _long += longStep) {
					nodeLoc[1] = _long + offset1;
					double d = 
					model.getLogLikelihood(nodeLoc, child1Loc, timeChild1) +
							model.getLogLikelihood(nodeLoc, child2Loc, timeChild2);
					logP[k] = d;
					max = Math.max(max, d);
					k++;
				}
			}
		} else {
			int parent = node.getParent().getNr();

			timeNode = node.getLength();
			parentLoc[0] = location.getValue(parent * 2);
			parentLoc[1] = location.getValue(parent * 2 + 1);
			
			
			double center0 = (child1Loc[0] * timeChild1 + child2Loc[0] * timeChild2 + parentLoc[0] * timeNode) / (timeChild1 + timeChild2 + timeNode);
			double center1 = (child1Loc[1] * timeChild1 + child2Loc[1] * timeChild2 + parentLoc[1] * timeNode) / (timeChild1 + timeChild2 + timeNode);
			int pos0 = (int) Math.round((center0 - minLat - offset0) / latStep);
			int pos1 = (int) Math.round((center1 - minLong - offset1) / longStep);
			double d  = calcLogP(pos0, pos1);
			logP[((pos0 + nY) % nY)  * nX + (pos1 + nX) % nX] = d;

			max = d;
			
			boolean progress = true;
			for (int i = 2; progress && i < nX; i += 2) {
				progress = false;
				for (int j = 0; j < i; j++) {
					int _lat = pos0 + i/2;
					int _long = pos1 - i/2 + j;
					d = calcLogP(_lat, _long);
					logP[((_lat + nY) % nY)  * nX + (_long + nX) % nX] = d;
					max = Math.max(max, d);
					if (d- max > THRESHOLD) {progress = true;} 

					_lat = pos0 - i/2;
					_long = pos1 + i/2 - j;
					d = calcLogP(_lat, _long);
					logP[((_lat + nY) % nY)  * nX + (_long + nX) % nX] = d;
					max = Math.max(max, d);
					if (d- max > THRESHOLD) {progress = true;} 

					_lat = pos0 - i/2 + j;
					_long = pos1 - i/2;
					d = calcLogP(_lat, _long);
					logP[((_lat + nY) % nY)  * nX + (_long + nX) % nX] = d;
					max = Math.max(max, d);
					if (d- max > THRESHOLD) {progress = true;} 

					_lat = pos0 + i/2 - j;
					_long = pos1 + i/2;
					d = calcLogP(_lat, _long);
					logP[((_lat + nY) % nY)  * nX + (_long + nX) % nX] = d;
					max = Math.max(max, d);
					if (d- max > THRESHOLD) {progress = true;} 
				
				}
			}
			
//			int k = 0;
//			for (double lat = minLat; lat < maxLat; lat += latStep) {
//				nodeLoc[0] = lat + offset0;
//				for (double _long = minLong; _long < maxLong; _long += longStep) {
//					nodeLoc[1] = _long + offset1;
//					logP[k] = model.getLogLikelihood(parentLoc, nodeLoc, timeNode) +
//							model.getLogLikelihood(nodeLoc, child1Loc, timeChild1) +
//							model.getLogLikelihood(nodeLoc, child2Loc, timeChild2);
//					k++;
//				}
//			}
		}
//		if (max == Double.NEGATIVE_INFINITY) {
//			// numeric issues, give up
//			return Double.NEGATIVE_INFINITY;
//		}
		for (int i = 0; i < logP.length; i++) {
			if (logP[i] < 0) {
				logP[i] = Math.exp(logP[i]-max);
			}
		}
		int k = Randomizer.randomChoicePDF(logP);
		int pos0 = (k / nX);
		double newLat = minLat + pos0 * latStep + offset0;
		int pos1 = (k % nX);
		double newLong = minLong + pos1 * longStep + offset1;
		location.setValue(iNode * 2, newLat);
		location.setValue(iNode * 2 + 1, newLong);
		return 0;
	}

	private double calcLogP(int _lat, int _long) {
		nodeLoc[0] = minLat + _lat * latStep + offset0;
		nodeLoc[1] = minLong + _long * longStep + offset1;
		return model.getLogLikelihood(parentLoc, nodeLoc, timeNode) +
				model.getLogLikelihood(nodeLoc, child1Loc, timeChild1) +
				model.getLogLikelihood(nodeLoc, child2Loc, timeChild2);
	}

}
