package beast.continuous;


import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.AlignmentFromTraitMap;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.ContinuousSubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeTraitMap;
import beast.util.Randomizer;

@Description("Approximate likelihood by MAP approximation of internal states")
public class ApproxMultivariateTraitLikelihood extends GenericTreeLikelihood {
	public Input<Boolean> scaleByBranchLengthInput = new Input<Boolean>("scale", "scale by branch lengths for initial position", true);
	
	ContinuousSubstitutionModel substModel;
	TreeInterface tree;
	BranchRateModel clockModel;
	
	double [][] position;
	double [] branchLengths;
	double [] sumLengths;
	double [] parentweight;

	boolean needsUpdate = true;
	boolean scaleByBranchLength;
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		clockModel = branchRateModelInput.get();
		SiteModel siteModel = (SiteModel) siteModelInput.get();
		substModel = (ContinuousSubstitutionModel) siteModel.substModelInput.get();
		tree = treeInput.get();
		scaleByBranchLength = scaleByBranchLengthInput.get();
		
		// initialise leaf positions
		position = new double[tree.getNodeCount()][2];
		AlignmentFromTraitMap data = (AlignmentFromTraitMap) dataInput.get();
		TreeTraitMap traitMap = data.getTraitMap(); 
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			position[i] = traitMap.getTrait(tree, nodes[i]);
		}
		branchLengths = new double[tree.getNodeCount()];
		sumLengths = new double[tree.getNodeCount()];
		parentweight = new double[tree.getNodeCount()];
		
	}
	
	@Override
    public double getCurrentLogP() {
        double logP = Double.NaN;
		try {
			logP = calculateLogP();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return logP;
    }

	
	@Override
	public double calculateLogP() throws Exception {
		logP = 0.0;
		calcBranchLengths();
		caclPositions();
		logP = calcLogP();
		needsUpdate = false;
		return logP;
	}

	void calcBranchLengths() {
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

	/** traverse tree **/
	double calcLogP() {
		double logP = 0;
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				logP += substModel.getLogLikelihood(
						position[node.getParent().getNr()], 
						position[node.getNr()], 
						branchLengths[node.getNr()]);
			}
		}
		return logP;
	}
	
	void caclPositions() {
		final double EPSILON = 1e-8;
		
		initByMean(tree.getRoot());
		resetMeanDown(tree.getRoot());
		
		if (scaleByBranchLength) {
			for (int i = 0; i < 50; i++) {
				double [][] oldPosition = new double[tree.getNodeCount()][2];
				for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
					oldPosition[j][0] = position[j][0];
					oldPosition[j][1] = position[j][1];
				}
				resetMeanDown(tree.getRoot());
				resetMeanUp(tree.getRoot());
				
				double max = 0;
				for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
					double delta0 = oldPosition[j][0] - position[j][0];
					double delta1 = oldPosition[j][1] - position[j][1];
					max = Math.max(max, Math.max(Math.abs(delta0), Math.abs(delta1)));
				}
				if (max < EPSILON) {
					// System.err.println(i + " maxdelta = " + max);
					break;
				}
			}
		}
	}
	
	void initByMean(Node node) {
		if (!node.isLeaf()) {
			initByMean(node.getLeft());
			initByMean(node.getRight());
			
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			setHalfWayPosition(nodeNr, child1, child2);
		}
	}		
	
	/** bottom up recalculation **/
	void resetMeanUp(Node node) {
		if (!node.isLeaf()) {
			resetMeanUp(node.getLeft());
			resetMeanUp(node.getRight());
			
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (node.isRoot()) {
				setHalfWayPosition(nodeNr, child1, child2);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent);
			}
		}
	}

	/** top down recalculation **/
	void resetMeanDown(Node node) {
		if (!node.isLeaf()) {
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (node.isRoot()) {
				//setHalfWayPosition(nodeNr, child1, child2);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent);
			}

			resetMeanDown(node.getLeft());
			resetMeanDown(node.getRight());
		}
	}

	void setHalfWayPosition(int nodeNr, int child1, int child2) {
		// start in weighted middle of the children
		

		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double len = b1 + b2;
			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2) / len;
			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2) / len;
			//position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
			//position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
		} else{
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			if (tree.getNode(nodeNr).isRoot()) {
				double len = b1 + b2;
				b1 = b1 / len;
				b2 = b2 / len;
				double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
				b1 /= w;
				b2 /= w;
			} else {
				double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
				double len = b1 + b2 + p;
				b1 /= len;
				b2 /= len;
				p /= len;
				double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
				b1 /= w;
				b2 /= w;
				p /= w;
				parentweight[nodeNr] = p;
			}
			position[nodeNr][0] = position[child1][0] * b1 + position[child2][0] * b2;
			position[nodeNr][1] = position[child1][1] * b1 + position[child2][1] * b2;
		}
		if (MAX_ITER <= 0) {return;}

		// optimise for subst model
		double [] newPos = new double[2];
		newPos[0] = position[nodeNr][0]; 
		newPos[1] = position[nodeNr][1]; 

		double logP = 
				substModel.getLogLikelihood(position[nodeNr], position[child1], branchLengths[child1]) +
				substModel.getLogLikelihood(position[nodeNr], position[child2], branchLengths[child2]);
		
		// optimise by random walk
		int i = 0;
		double epsilon = 0.5;
		while (i < MAX_ITER && epsilon > MIN_EPSILON) {
//			boolean progress = false;
//			newPos[0] += epsilon;
//			double newLogP = 
//					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
//					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]);
//			double dX = newLogP - logP;
//			if (dX > 0) {
//				logP = newLogP;
//				progress = true;
//			} else {
//				newPos[0] -= epsilon;
//			}
//			newPos[1] += epsilon;
//			newLogP = 
//					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
//					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]);
//			double dY= newLogP - logP;
//			if (dY > 0) {
//				logP = newLogP;
//				progress = true;
//			} else {
//				newPos[1] -= epsilon;
//			}
//			if (progress) {
//				position[nodeNr][0] = newPos[0]; 
//				position[nodeNr][1] = newPos[1];
//				epsilon *= 1.5;
//			} else {
//				epsilon /= 1.5;
//			}
			
			newPos[0] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
			newPos[1] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
			double newLogP = 
					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]);
			if (newLogP > logP) {
				logP = newLogP;
				position[nodeNr][0] = newPos[0]; 
				position[nodeNr][1] = newPos[1]; 
			} else {
				newPos[0] = position[nodeNr][0]; 
				newPos[1] = position[nodeNr][1]; 
			}
			i++;
		}
	}

	void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
			double len = b1 + b2 + p;
			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2 + position[parent][0] * p) / len;
			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2 + position[parent][1] * p) / len;

			//position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2] + position[parent][0] * branchLengths[nodeNr]) / sumLengths[nodeNr];
			//position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2] + position[parent][1] * branchLengths[nodeNr]) / sumLengths[nodeNr];
		} else {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
			double len = b1 + b2 + p;
			b1 /= len;
			b2 /= len;
			p /= len;
			double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
			p /= w;
			position[nodeNr][0] += position[parent][0] * p;
			position[nodeNr][1] += position[parent][1] * p;

		}

		if (MAX_ITER <= 0) {return;}
		
		// optimise for subst model
		double [] newPos = new double[2];
		newPos[0] = position[nodeNr][0]; 
		newPos[1] = position[nodeNr][1]; 

		double logP = 
				substModel.getLogLikelihood(position[nodeNr], position[child1], branchLengths[child1]) +
				substModel.getLogLikelihood(position[nodeNr], position[child2], branchLengths[child2]) +
				substModel.getLogLikelihood(position[parent], position[nodeNr], branchLengths[nodeNr]);
		
		// optimise by random walk
		int i = 0;
		double epsilon = 0.5;
		while (i < MAX_ITER && epsilon > MIN_EPSILON) {
//			boolean progress = false;
//			newPos[0] += epsilon;
//			double newLogP = 
//					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
//					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]) +
//					substModel.getLogLikelihood(position[parent], newPos, branchLengths[nodeNr]);
//			double dX = newLogP - logP;
//			if (dX > 0) {
//				logP = newLogP;
//				progress = true;
//			} else {
//				newPos[0] -= epsilon;
//			}
//			newPos[1] += epsilon;
//			newLogP = 
//					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
//					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]) +
//					substModel.getLogLikelihood(position[parent], newPos, branchLengths[nodeNr]);
//			double dY= newLogP - logP;
//			if (dY > 0) {
//				logP = newLogP;
//				progress = true;
//			} else {
//				newPos[1] -= epsilon;
//			}
//			if (progress) {
//				position[nodeNr][0] = newPos[0]; 
//				position[nodeNr][1] = newPos[1];
//				epsilon *= 1.5;
//			} else {
//				epsilon /= 1.5;
//			}

			newPos[0] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
			newPos[1] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
			double newLogP = 
					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]) +
					substModel.getLogLikelihood(position[parent], newPos, branchLengths[nodeNr]);
			if (newLogP > logP) {
				logP = newLogP;
				position[nodeNr][0] = newPos[0]; 
				position[nodeNr][1] = newPos[1]; 
			} else {
				newPos[0] = position[nodeNr][0]; 
				newPos[1] = position[nodeNr][1]; 
			}
			i++;
		}
	}
	final static int MAX_ITER = 0;
	final static double MIN_EPSILON = 0.001;
	
	@Override
	public boolean isStochastic() {
		return true;
	}

	public double[] getPostion(int iDim) {
		if (needsUpdate) {
			try {
				calculateLogP();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return position[iDim];
	}
	
	@Override
	public void store() {
		needsUpdate = true;
		super.store();
	}
	@Override
	public void restore() {
		needsUpdate = true;
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		return super.requiresRecalculation();
	}
	
}
