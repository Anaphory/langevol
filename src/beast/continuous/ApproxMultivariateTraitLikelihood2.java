package beast.continuous;


import beast.core.Description;
import sphericalGeo.SphericalDiffusionModel;

@Description("Approximate likelihood by MAP approximation of internal states")
public class ApproxMultivariateTraitLikelihood2 extends ApproxMultivariateTraitLikelihood {
	
	double [][] sphereposition;

	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		sphereposition = new double[position.length][3];
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
		}
		
		parentweight = new double[position.length];
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			parentweight[i] = 0.0;
		}
	}
	
	
	void caclPositions() {
		final double EPSILON = 1e-8;
		initByMean(tree.getRoot());
		resetMeanDown(tree.getRoot());
		
			
		if (scaleByBranchLength) {
		double [][] oldPosition = new double[tree.getNodeCount()][3];
		for (int i = 0; i < 50; i++) {
			if (i % 5 == 0)
			for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
				oldPosition[j][0] = sphereposition[j][0];
				oldPosition[j][1] = sphereposition[j][1];
				oldPosition[j][2] = sphereposition[j][2];
			}

				resetMeanDown(tree.getRoot());
				resetMeanUp(tree.getRoot());

				if (i % 5 == 0) {
				double max = 0;
				for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
					double delta0 = oldPosition[j][0] - sphereposition[j][0];
					double delta1 = oldPosition[j][1] - sphereposition[j][1];
					double delta2 = oldPosition[j][2] - sphereposition[j][2];
					max = Math.max(max, Math.max(Math.abs(delta0), Math.max(Math.abs(delta1), Math.abs(delta2))));
				}
				if (max < EPSILON) {
					break;
				}
				}
					
			}
		}
		
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			position[i] = SphericalDiffusionModel.cartesian2Sperical(sphereposition[i], true);
		}

//		System.err.println("maxdelta2 = " + max);
		
	}
	

	
	@Override
	void setHalfWayPosition(int nodeNr, int child1, int child2) {
		// start in weighted middle of the children
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double len = b1 + b2;
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2) / len;
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2) / len;
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2) / len;
//			double len = 1.0/branchLengths[child1] + 1.0/branchLengths[child2];
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] / branchLengths[child1] + sphereposition[child2][0] / branchLengths[child2]) / len;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] / branchLengths[child1] + sphereposition[child2][1] / branchLengths[child2]) / len;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] / branchLengths[child1] + sphereposition[child2][2] / branchLengths[child2]) / len;
		} else{
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			if (tree.getNode(nodeNr).isRoot()) {
				double len = b1 + b2;
				b1 = b1 / len;
				b2 = b2 / len;
				double w = (1.0 + b1 * parentweight[child1] + b2 * parentweight[child2]);
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
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2);
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2);
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2);
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] + sphereposition[child2][0]) / 2.0;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] + sphereposition[child2][1]) / 2.0;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] + sphereposition[child2][2]) / 2.0;
		}
		normalise(sphereposition[nodeNr]);
	}

	@Override
	void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
			double len = b1 + b2 + p;
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2 + sphereposition[parent][0] * p) / len;
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2 + sphereposition[parent][1] * p) / len;
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2 + sphereposition[parent][2] * p) / len;

//			double len = 1.0/branchLengths[child1] + 1.0/branchLengths[child2] + 1.0/branchLengths[nodeNr];
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] / branchLengths[child1] + sphereposition[child2][0] / branchLengths[child2] + sphereposition[parent][0] / branchLengths[nodeNr]) / len;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] / branchLengths[child1] + sphereposition[child2][1] / branchLengths[child2] + sphereposition[parent][1] / branchLengths[nodeNr]) / len;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] / branchLengths[child1] + sphereposition[child2][2] / branchLengths[child2] + sphereposition[parent][2] / branchLengths[nodeNr]) / len;
		} else {
			double b1 = branchLengths[child1];
			double b2 = branchLengths[child2];
			double p = branchLengths[nodeNr];
			double len = b1 + b2 + p;
			b1 /= len;
			b2 /= len;
			p /= len;
			double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
			p /= w;
			sphereposition[nodeNr][0] += sphereposition[parent][0] * p;
			sphereposition[nodeNr][1] += sphereposition[parent][1] * p;
			sphereposition[nodeNr][2] += sphereposition[parent][2] * p;
		}
		normalise(sphereposition[nodeNr]);
	}


	private void normalise(double[] position) {
		double len = Math.sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
		position[0] /= len;
		position[1] /= len;
		position[2] /= len;
	}

}
