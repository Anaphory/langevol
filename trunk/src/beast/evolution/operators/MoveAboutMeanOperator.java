package beast.evolution.operators;


import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.AlignmentFromTraitMap;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeTraitMap;
import beast.util.Randomizer;

public class MoveAboutMeanOperator extends Operator {
	public Input<RealParameter> locationInput = new Input<RealParameter>("location", "latitude/longitude pairs representing location", Validate.REQUIRED);
	public Input<TreeInterface> treeInput = new Input<>("tree", "tree for whcih to infer geography", Validate.REQUIRED);
	public Input<AlignmentFromTraitMap> traitMapInput = new Input<AlignmentFromTraitMap>("trait", "map that defines the trait", Validate.REQUIRED);
    public Input<BranchRateModel.Base> branchRateModelInput = new Input<BranchRateModel.Base>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.", new StrictClockModel());
    public Input<Double> windowSizeInput =
            new Input<Double>("windowSize", "the size of the standard deviation for perturbing points", 5.0);

    
    public Input<Operator> operatorInput = new Input<Operator>("operator" ,"optional tree operator -- locations of filthy " +
            "nodes will get a new location");

    TreeInterface tree;
    RealParameter location;
    BranchRateModel.Base clockModel;
	Operator operator;
	
    double [][] position;
	double [] branchLengths;
	double [] sumLengths;
	double windowSize;
	
    @Override
	public void initAndValidate() throws Exception {
    	tree = treeInput.get();
    	clockModel = branchRateModelInput.get();
    	operator = operatorInput.get();
    	location = locationInput.get();
    	if (location.getMinorDimension1() != 2) {
    		throw new RuntimeException("Trait should have minor dimension = 2");
    	}
    	
		// initialise leaf positions
		position = new double[tree.getNodeCount()][2];
		AlignmentFromTraitMap data = traitMapInput.get();
		TreeTraitMap traitMap = data.getTraitMap(); 
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			position[i] = traitMap.getTrait(tree, nodes[i]);
		}
		branchLengths = new double[tree.getNodeCount()];
		sumLengths = new double[tree.getNodeCount()];
	}

	@Override
	public double proposal() {
		double logHR = 0;
		calcBranchLengths(tree, branchLengths, sumLengths, clockModel);
		caclMeanPositions(tree, branchLengths, sumLengths, position);
		
		if (operator != null) {
			// do operator, then do random walk on all filthy nodes
			logHR = operator.proposal();
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

		return logHR;
	}
	
	double traverse(Node node, boolean [] updated) {
		double logHR = 0.0;
		if (node.isLeaf()) {
			return 0.0;
		} else {
			boolean needsUpdate = false;
			for (Node child : node.getChildren()) {
				updated[0] = false;
				logHR += traverse(child, updated);
				if (updated[0]) {
					needsUpdate = true;
				}
			}
			if (needsUpdate || node.isDirty() == Tree.IS_FILTHY) {
				logHR += doproposal(node.getNr());
				updated[0] = true;
			}
		}
		return logHR;
	}

	double doproposal(final int nodeNr) {
		
		double orig0 = location.getValue(nodeNr * 2 + 0);
		double d0 = position[nodeNr][0] - orig0;
		double orig1 = location.getValue(nodeNr * 2 + 1);
		double d1 = position[nodeNr][1] - orig1;
		
		double logHR = - d0 * d0 / (2.0 * windowSize) - d1 * d1 / (2.0 * windowSize);
		
		d0 = Randomizer.nextGaussian() * windowSize;
		d1 = Randomizer.nextGaussian() * windowSize;
		position[nodeNr][0] += d0;
		position[nodeNr][1] += d1;

		
		logHR += + d0 * d0 / (2.0 * windowSize) + d1 * d1 / (2.0 * windowSize);

		// assign new values to location
		location.setValue(nodeNr * 2 + 0, position[nodeNr][0]);
		location.setValue(nodeNr * 2 + 1, position[nodeNr][1]);
		
		return logHR;
	}

	public static void calcBranchLengths(TreeInterface tree, double [] branchLengths, double [] sumLengths, BranchRateModel clockModel) {
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				branchLengths[node.getNr()] = node.getLength() * clockModel.getRateForBranch(node);
				if (!node.isLeaf()) {
					Node child1 = node.getLeft();
					Node child2 = node.getRight();
					sumLengths[node.getNr()] = branchLengths[node.getNr()] +
							child1.getLength() * clockModel.getRateForBranch(child1) +
							child2.getLength() * clockModel.getRateForBranch(child2);
				}
			}
		}
	}
	
	public static void caclMeanPositions(TreeInterface tree, double [] branchLengths, double [] sumLengths, double [][] position) {
		initByMean(tree.getRoot(),position, branchLengths);
		resetMean2(tree.getRoot(), position, branchLengths, sumLengths);
	}
	
	public static void initByMean(Node node, double [][] position, double [] branchLengths) {
		if (!node.isLeaf()) {
			initByMean(node.getLeft(), position, branchLengths);
			initByMean(node.getRight(), position, branchLengths);
			
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			setHalfWayPosition(nodeNr, child1, child2, position, branchLengths);
		}
	}
	
	/** top down recalculation **/
	public static void resetMean2(Node node, double [][] position, double [] branchLengths, double [] sumLengths) {
		if (!node.isLeaf()) {
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (node.isRoot()) {
				setHalfWayPosition(nodeNr, child1, child2, position, branchLengths);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent, position, branchLengths, sumLengths);
			}

			resetMean2(node.getLeft(), position, branchLengths, sumLengths);
			resetMean2(node.getRight(), position, branchLengths, sumLengths);
		}
	}

	public static void setHalfWayPosition(int nodeNr, int child1, int child2, double [][] position, double [] branchLengths) {
		// start in weighted middle of the children
		position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
		position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
		
		// optimise for subst model
		double [] newPos = new double[2];
		newPos[0] = position[nodeNr][0]; 
		newPos[1] = position[nodeNr][1]; 

	}

	public static void setHalfWayPosition(int nodeNr, int child1, int child2, int parent, double [][] position, double [] branchLengths, double [] sumLengths) {
		// start in weighted middle of the children and parent location
		position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2] + position[parent][0] * branchLengths[nodeNr]) / sumLengths[nodeNr];
		position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2] + position[parent][1] * branchLengths[nodeNr]) / sumLengths[nodeNr];

		// optimise for subst model
		double [] newPos = new double[2];
		newPos[0] = position[nodeNr][0]; 
		newPos[1] = position[nodeNr][1]; 
	}
}
