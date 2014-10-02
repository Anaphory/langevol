package beast.continuous;



import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

@Description("Helper class for logging locations from ApproxMultivariateTraitLikelihood")
public class TraitFunction extends RealParameter {
	public Input<ApproxMultivariateTraitLikelihood> likelihoodInput = new Input<ApproxMultivariateTraitLikelihood>("likelihood", "trait likelihood to be logged", Validate.REQUIRED);

	ApproxMultivariateTraitLikelihood likelihood;
	Tree tree;
	
	@Override
	public void initAndValidate() throws Exception {
		likelihood = likelihoodInput.get();
        tree =  ((Tree) (likelihood.treeInput.get()));
	}

	@Override
	public int getDimension() {
		return tree.getNodeCount() * 2;
	}

	@Override
	public int getMinorDimension1() {
		return 2;
	}

	@Override
	public double getArrayValue() {
		return 0;
	}

	@Override
	public Double getMatrixValue(int i, int j) {
		return likelihood.getPostion(i)[j];
	}

}
