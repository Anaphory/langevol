package beast.geo;

import beast.evolution.tree.Node;
import beast.util.Randomizer;


public class LocationAllGibbsOperator extends LocationGibbsOperator {


	final static int MAX_NR_PARTIALS = 20;
	double [][] partials;
	double [][][] matrices;
	int [][] graphNodeIDs;
	
	
	double [][] transitionprobabilities;
	int n;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		n = graph.getSize();
		partials = new double[tree.getNodeCount()][];
		graphNodeIDs = new int[tree.getNodeCount()][];
		matrices = new double[tree.getNodeCount()][][];
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			graphNodeIDs[i] = new int[0];
			//matrices[i] = new double[0][0];
		}
		for (int i = tree.getLeafNodeCount(); i < partials.length; i++) {
			partials[i] = new double[n];
			graphNodeIDs[i] = new int[MAX_NR_PARTIALS];
			matrices[i] = new double[MAX_NR_PARTIALS][n];
		}
		int n = graph.getSize();
		transitionprobabilities = new double[n][n*n];
		
	}
	
	@Override
	public double proposal() {
		traverse(tree.getRoot());
		sample(tree.getRoot());
		
		return Double.POSITIVE_INFINITY;
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
        			probs[i] = myMatrices[i][parent] * myPartials[myGgraphNodeIDs[i]]; 
        		}
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
		double timeRight = branchRateLeft * (node.getHeight() - right.getHeight());
		
		int iNode = node.getNr();
		int iLeft = left.getNr();
		int iRight = right.getNr();
		
		int [] graphNodeIDsLeft = graphNodeIDs[iLeft];
		double [][] pLeft = matrices[iLeft]; 
		for (int i = 0; i < graphNodeIDsLeft.length; i++) {
			int iGraphNode = graphNodeIDsLeft[i];
			GraphNode gnode = graph.nodes.get(iGraphNode);
	        pLeft[i] = model.getLogTransitionProbabilities(gnode, timeLeft);
			GraphDistanceBasedDiffusionModel.logPtoP(pLeft[i]);
		}

		int [] graphNodeIDsRight = graphNodeIDs[iRight];
		double [][] pRight = matrices[iRight]; 
		for (int i = 0; i < graphNodeIDsRight.length; i++) {
			int iGraphNode = graphNodeIDsRight[i];
			GraphNode gnode = graph.nodes.get(iGraphNode);
	        pRight[i] = model.getLogTransitionProbabilities(gnode, timeRight);
			GraphDistanceBasedDiffusionModel.logPtoP(pRight[i]);
		}

		
		// combine
		double [] nodePartials = partials[iNode];
		double [] leftPartials = partials[iLeft];
		double [] rightPartials = partials[iRight];
		for (int i = 0; i < n; i++) {
			double sum1 = 0.0;
			if (left.isLeaf()) {
				double [] matrix = new double[1];
				model.getTransitionProbabilities(left, left.getHeight(), node.getHeight(), branchRateLeft, matrix);
				sum1 = matrix[0];
			} else {
				for (int j = 0; j < graphNodeIDsLeft.length; j++) {
					sum1 += pLeft[j][i] * leftPartials[j];
				}
			}
			double sum2 = 0.0;
			if (right.isLeaf()) {
				double [] matrix = new double[1];
				model.getTransitionProbabilities(right, right.getHeight(), node.getHeight(), branchRateRight, matrix);
				sum2 = matrix[0];
			} else {
				for (int j = 0; j < graphNodeIDsRight.length; j++) {
					sum2 += pRight[j][i] * rightPartials[j];
				}
			}
			nodePartials[i] = sum1 * sum2;
		}
		
		// find top MAX_NR_PARTIALS partials 
		int [] myGraphNodeIDs = graphNodeIDs[iNode];
		int k = myGraphNodeIDs.length;
		double [] myPartials = new double[k];
		
		double min = Double.POSITIVE_INFINITY;
		int iMin = -1;
		for (int i = 0; i < k; i++) {
			myGraphNodeIDs[i] = i;
			myPartials[i] = nodePartials[i];
			if (nodePartials[i] < min) {
				min = nodePartials[i];
				iMin = i;
			}
		}
		for (int i = myGraphNodeIDs.length; i < n; i++) {
			if (nodePartials[i] > min) {
				myPartials[iMin] = nodePartials[i];
				myGraphNodeIDs[iMin] = i;
				min = Double.POSITIVE_INFINITY;
				iMin = -1;
				for (int j = 0; j < k; j++) {
					if (myPartials[i] < min) {
						min = myPartials[i];
						iMin = i;
					}
				}
			}
		}
		
	}

}
