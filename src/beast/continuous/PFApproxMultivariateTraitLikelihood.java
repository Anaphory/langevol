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

@Description("Approximate likelihood by particle filter approximation")
public class PFApproxMultivariateTraitLikelihood extends GenericTreeLikelihood {
	public Input<Integer> nrOfParticlesInput = new Input<Integer>("nrOfParticles", "number of particles to use", 10);
	public Input<Integer> nrOfIterationsInput = new Input<Integer>("nrOfIterations", "number of iterations to run the particle filter", 10);
	
	double epsilon = 0.5;

	class LeafParticleSet {
		int particleCount;
		int nodeNr;
		double [] initialPosition;
		
		LeafParticleSet(int nodeNr, int particleCount, double [] initialPosition) {
			this.nodeNr = nodeNr;
			this.particleCount = particleCount;
			this.initialPosition = initialPosition;
		}
		double [] getPosition(int iParticle) {
			return initialPosition;
		}
		void init(double [] position) {
			throw new RuntimeException("thou shalt not call init on LeafParticleSet");
		}
		
		void resample(ContinuousSubstitutionModel substModel) {}
		public void randomize() {}
	} // class LeafParticleSet

	
	class ParticleSet extends LeafParticleSet {
		double [][] pposition;
		double [] logP;
		
		ParticleSet(int nodeNr, int particleCount, double [] initialPosition) {
			super(nodeNr, particleCount, initialPosition);
			pposition = new double[particleCount][2];
			logP = new double[particleCount];
		}
		
		@Override
		double [] getPosition(int iParticle) {
			return pposition[iParticle];
		}

		@Override
		void init(double [] position) {
			for (double [] pos : pposition) {
				pos[0] = position[0];
				pos[1] = position[1];
			}
		}

		@Override
		void resample(ContinuousSubstitutionModel substModel) {
			Node node = tree.getNode(nodeNr);
			int iChild1 = node.getLeft().getNr();
			LeafParticleSet child1 = particleSets[iChild1];
			int iChild2 = node.getRight().getNr();
			LeafParticleSet child2 = particleSets[iChild2];
			
			if (node.isRoot()) {
				for (int i = 0; i < particleCount; i++) {
					double [] nodePosition = getPosition(i);
					logP[i] = substModel.getLogLikelihood(child1.getPosition(i), nodePosition, branchLengths[iChild1]) +
							substModel.getLogLikelihood(child2.getPosition(i), nodePosition, branchLengths[iChild2]);
				}
			} else {
				int iParent = node.getParent().getNr();
				LeafParticleSet parent = particleSets[iParent];
				for (int i = 0; i < particleCount; i++) {
					double [] nodePosition = getPosition(i);
					logP[i] = substModel.getLogLikelihood(child1.getPosition(i), nodePosition, branchLengths[iChild1]) +
							substModel.getLogLikelihood(child2.getPosition(i), nodePosition, branchLengths[iChild2]) +
							substModel.getLogLikelihood(nodePosition, parent.getPosition(i), branchLengths[node.getNr()]);
				}
			}
			double max = logP[0];
			for (double d : logP) {
				max = Math.max(d,  max);
			}
			for (int i = 0; i < particleCount; i++) {
				logP[i] = Math.exp(logP[i] - max);
			}
			double [][] newposition = new double[particleCount][2];
			for (int i = 0; i < particleCount; i++) {
				int k = Randomizer.randomChoicePDF(logP);
				newposition[i][0] = pposition[k][0];
				newposition[i][1] = pposition[k][1];
			}
			pposition = newposition;
		}
		
		@Override
		public void randomize() {
			for (int i = 0; i < particleCount; i++) {
				pposition[i][0] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
				pposition[i][1] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
			}			
		}

	} // class ParticleSet
	
	
	
	
	ContinuousSubstitutionModel substModel;
	TreeInterface tree;
	BranchRateModel clockModel;
	int particleCount;
	int iterationCount;
	
	double [][] position;
	double [] branchLengths;
	double [] sumLengths;
	boolean needsUpdate = true;
	
	LeafParticleSet [] particleSets;
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		clockModel = branchRateModelInput.get();
		SiteModel siteModel = (SiteModel) siteModelInput.get();
		substModel = (ContinuousSubstitutionModel) siteModel.substModelInput.get();
		tree = treeInput.get();
		
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
		
		
		particleCount = nrOfParticlesInput.get();
		iterationCount = nrOfIterationsInput.get();

		// set up particles
		particleSets = new LeafParticleSet[tree.getNodeCount()];
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			particleSets[i] = new LeafParticleSet(i, particleCount, position[i]);
		}
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			particleSets[i] = new ParticleSet(i, particleCount, position[i]);
		}
		
	}
	
	@Override
	public double calculateLogP() throws Exception {
		logP = 0.0;
		calcBranchLengths();
		setUpInitialPositions();
		for (int i = 0; i < iterationCount; i++) {
			resampleUp(tree.getRoot());
			resampleDown(tree.getRoot());
		}
		
		
		logP = calcLogP();
		needsUpdate = false;
		return logP;
	}

	void resampleUp(Node node) {
		if (!node.isLeaf()) {
			particleSets[node.getNr()].randomize();
			particleSets[node.getNr()].resample(substModel);

			resampleUp(node.getLeft());
			resampleUp(node.getRight());
		}
	}

	void resampleDown(Node node) {
		if (!node.isLeaf()) {
			resampleDown(node.getLeft());
			resampleDown(node.getRight());
			particleSets[node.getNr()].randomize();
			particleSets[node.getNr()].resample(substModel);

		}
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
		// logP is mean of logP's of the particles(?)
		for (int i = 0; i < particleCount;  i++) {
			for (Node node : tree.getNodesAsArray()) {
				if (!node.isRoot()) {
					logP += substModel.getLogLikelihood(
							particleSets[node.getParent().getNr()].getPosition(i), 
							particleSets[node.getNr()].getPosition(i), 
							branchLengths[node.getNr()]);
				}
			}
		}
		logP /= particleCount;
		return logP;
	}
	
	void setUpInitialPositions() {
		initByMean(tree.getRoot());
		resetMean2(tree.getRoot());
		
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			particleSets[i].init(position[i]);
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
	void resetMean(Node node) {
		if (!node.isLeaf()) {
			resetMean(node.getLeft());
			resetMean(node.getRight());
			
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
	void resetMean2(Node node) {
		if (!node.isLeaf()) {
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (node.isRoot()) {
				setHalfWayPosition(nodeNr, child1, child2);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent);
			}

			resetMean2(node.getLeft());
			resetMean2(node.getRight());
		}
	}

	private void setHalfWayPosition(int nodeNr, int child1, int child2) {
		// start in weighted middle of the children
		position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
		position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
		
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

	private void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
		position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2] + position[parent][0] * branchLengths[nodeNr]) / sumLengths[nodeNr];
		position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2] + position[parent][1] * branchLengths[nodeNr]) / sumLengths[nodeNr];

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
