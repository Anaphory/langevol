package beast.geo;

import java.util.Arrays;

import beast.continuous.DiscretePFApproxMultivariateTraitLikelihood;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;


public class LocationAllGibbsOperator extends Operator {
    public Input<IntegerParameter> parameterInput =
            new Input<IntegerParameter>("parameter", "the parameter to operate a random walk on.", Validate.REQUIRED);
	
	public Input<DiscretePFApproxMultivariateTraitLikelihood> likelihoodInput = new Input<DiscretePFApproxMultivariateTraitLikelihood>("likelihood2", "DiscretePFApproxMultivariateTraitLikelihood for logging locations", Validate.REQUIRED);

	DiscretePFApproxMultivariateTraitLikelihood likelihood;
    IntegerParameter param;
    GraphDistanceBasedDiffusionModel model;
    Graph graph;
    Tree tree;
    BranchRateModel clockModel;
	int n;

	final static int MAX_NR_PARTIALS = 20;
	double [][] partials;
	double [][][] matrices;
	int [][] graphNodeIDs;
	
	
	double [][] transitionprobabilities;
    IntegerParameter location;

	@Override
	public void initAndValidate() {
		//super.initAndValidate();
		param = parameterInput.get();
		likelihood = likelihoodInput.get();
		tree = (Tree) likelihood.treeInput.get();
        model = likelihood.modelInput.get();
        graph =  model.graphInput.get();
		n = graph.getSize();
		
        clockModel = likelihood.branchRateModelInput.get();
        location = model.locationInput.get();
        if (location != param) {
        	throw new RuntimeException("location != param");
        }

		
//		partials = new double[tree.getNodeCount()][];
//		graphNodeIDs = new int[tree.getNodeCount()][];
//		matrices = new double[tree.getNodeCount()][][];
//		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
//			graphNodeIDs[i] = new int[0];
//			matrices[i] = new double[1][];
//		}
//		for (int i = tree.getLeafNodeCount(); i < partials.length; i++) {
//			partials[i] = new double[n];
//			graphNodeIDs[i] = new int[MAX_NR_PARTIALS];
//			matrices[i] = new double[MAX_NR_PARTIALS][n];
//		}
//		int n = graph.getSize();
//		transitionprobabilities = new double[n][n*n];
		
	}
	
	@Override
	public double proposal() {
		try {
			likelihood.calculateLogP();
			int nrOfLeafs = tree.getLeafNodeCount();
			for (int i = nrOfLeafs; i < nrOfLeafs*2-1; i++) {
				// TODO: Figure out how to fix this method call
				// Filocation.setValue(i, likelihood.getPostion(i));
			}
			return Double.POSITIVE_INFINITY;
		} catch (Exception e) {
			e.printStackTrace();
			return Double.NEGATIVE_INFINITY;
		}
		
		
//		traverse(tree.getRoot());
//		sample(tree.getRoot());
		
//		return Double.POSITIVE_INFINITY;
	}

	private void sample(Node node) {
		int iNode = node.getNr();
		
        if (!node.isLeaf()) {

        	if (node.isRoot()) {
        		int loc = Randomizer.randomChoicePDF(partials[iNode]);
        		location.setValue(iNode, loc);
        	} else {
        		int parent = location.getValue(node.getParent().getNr());
        		int [] myGgraphNodeIDs = graphNodeIDs[iNode];
        		double [][] myMatrices = matrices[iNode];
        		double [] myPartials = partials[iNode];
        		double [] probs = new double [myGgraphNodeIDs.length];
        		for (int i = 0; i < myGgraphNodeIDs.length; i++) {
        			double a = myMatrices[i][parent];
        			double b = myPartials[myGgraphNodeIDs[i]];
        			if (a > 0 && b > 0) {
        				probs[i] = Math.log(a) + Math.log(b);
        			}
        		}
        		model.logPtoP(probs);
        		int locID = Randomizer.randomChoicePDF(probs);
        		int loc = myGgraphNodeIDs[locID];
        		location.setValue(iNode, loc);
        	}
        	
            sample(node.getLeft());
            sample(node.getRight());
        }
	}

	void traverse(Node node) {
	        if (!node.isLeaf()) {
	            final Node child1 = node.getLeft();
	            traverse(child1);

	            final Node child2 = node.getRight();
	            traverse(child2);
                calculatePartials(child1, child2, node);
	        }
	} // traverse

	private void calculatePartials(Node left, Node right, Node node) {
		double branchRateLeft = clockModel.getRateForBranch(left); 
		double timeLeft = branchRateLeft * (node.getHeight() - left.getHeight());
		double branchRateRight = clockModel.getRateForBranch(right); 
		double timeRight = branchRateRight * (node.getHeight() - right.getHeight());
		
		int iNode = node.getNr();
		int iLeft = left.getNr();
		int iRight = right.getNr();
		
		int [] graphNodeIDsLeft = graphNodeIDs[iLeft];
		double [][] pLeft = matrices[iLeft];
		if (left.isLeaf()) {
			GraphNode gnode = graph.nodes.get(location.getValue(iLeft));
	        double [] logP = model.getLogTransitionProbabilities(gnode, timeLeft);
	        model.logPtoP(logP);
	        pLeft[0] = logP;
		} else {
			for (int i = 0; i < graphNodeIDsLeft.length; i++) {
				int iGraphNode = graphNodeIDsLeft[i];
				GraphNode gnode = graph.nodes.get(iGraphNode);
		        pLeft[i] = model.getLogTransitionProbabilities(gnode, timeLeft);
				GraphDistanceBasedDiffusionModel.logPtoP(pLeft[i]);
			}
		}


		
		int [] graphNodeIDsRight = graphNodeIDs[iRight];
		double [][] pRight = matrices[iRight];
		if (right.isLeaf()) {
			GraphNode gnode = graph.nodes.get(location.getValue(iRight));
	        double [] logP = model.getLogTransitionProbabilities(gnode, timeRight);
	        model.logPtoP(logP);
	        pRight[0] = logP;
		} else {
			for (int i = 0; i < graphNodeIDsRight.length; i++) {
				int iGraphNode = graphNodeIDsRight[i];
				GraphNode gnode = graph.nodes.get(iGraphNode);
		        pRight[i] = model.getLogTransitionProbabilities(gnode, timeRight);
				GraphDistanceBasedDiffusionModel.logPtoP(pRight[i]);
			}
		}

		int leftGraphNode= location.getValue(iLeft);
		int rightGraphNode= location.getValue(iRight);
		
		// combine
		double [] nodePartials = partials[iNode];
		double [] leftPartials = partials[iLeft];
		double [] rightPartials = partials[iRight];
		for (int i = 0; i < n; i++) {
			double sum1 = 0.0;
			if (left.isLeaf()) {
				sum1 = pLeft[0][i];
			} else {
				for (int j = 0; j < graphNodeIDsLeft.length; j++) {
					sum1 += pLeft[j][i] * leftPartials[graphNodeIDsLeft[j]];
				}
			}
			double sum2 = 0.0;
			if (right.isLeaf()) {
				sum2 = pRight[0][i];
			} else {
				for (int j = 0; j < graphNodeIDsRight.length; j++) {
					sum2 += pRight[j][i] * rightPartials[graphNodeIDsRight[j]];
				}
			}
			if (sum1 > 0 && sum2 > 0) {
				nodePartials[i] = Math.log(sum1) + Math.log(sum2);
			}
//			if (Double.isNaN(nodePartials[i])) {
//				nodePartials[i] = 0.0;
//			}
		}

		model.logPtoP(nodePartials);

		double [] pParent = null;
		if (!node.isRoot()) {
			Node parent = node.getParent();
			int iParent = parent.getNr();
			GraphNode gnode = graph.nodes.get(location.getValue(iParent));
			double branchRateParent = clockModel.getRateForBranch(parent);
			double timeParent = branchRateParent * (parent.getHeight() - node.getHeight());
	        pParent = model.getLogTransitionProbabilities(gnode, timeParent);
	        model.logPtoP(pParent);
	        for (int i = 0; i < n; i++) {
	        	pParent[i] *= nodePartials[i];
	        }
		} else {
			pParent = nodePartials;
		}
		
		
		// find top MAX_NR_PARTIALS partials 
		int [] myGraphNodeIDs = graphNodeIDs[iNode];
		int k = myGraphNodeIDs.length;
		double [] myPartials = new double[k];
		
		double min = Double.POSITIVE_INFINITY;
		int iMin = -1;
		
		double [] tmp = new double[MAX_NR_PARTIALS];
		for (int i = 0; i < k; i++) {
			myGraphNodeIDs[i] = i;
			myPartials[i] = nodePartials[i];
			tmp[i] = pParent[i];
			if (pParent[i] < min) {
				min = pParent[i];
				iMin = i;
			}
		}
		for (int i = myGraphNodeIDs.length; i < n; i++) {
//			if (nodePartials[i] > min) {
//				myPartials[iMin] = nodePartials[i];
//				myGraphNodeIDs[iMin] = i;
//				min = Double.POSITIVE_INFINITY;
//				iMin = -1;
//				for (int j = 0; j < k; j++) {
//					if (myPartials[j] < min) {
//						min = myPartials[j];
//						iMin = j;
//					}
//				}
//			}
			if (pParent[i] > min) {
				myPartials[iMin] = nodePartials[i];
				tmp[iMin] = pParent[i];
				myGraphNodeIDs[iMin] = i;
				min = Double.POSITIVE_INFINITY;
				iMin = -1;
				for (int j = 0; j < k; j++) {
					if (tmp[j] < min) {
						min = tmp[j];
						iMin = j;
					}
				}
			}
		}
		
		
		
		double max = 0.0;
		for (double d : myPartials) {
			max = Math.max(max, d);
		}
		//System.out.println(max);
		for (int i = 0; i < n; i++) {
			nodePartials[i] /= max;
		}
	}

}
