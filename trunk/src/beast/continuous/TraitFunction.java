package beast.continuous;


import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Tree;

@Description("Helper class for logging locations from ApproxMultivariateTraitLikelihood")
public class TraitFunction extends BEASTObject implements Function {
	public Input<ApproxMultivariateTraitLikelihood> likelihoodInput = new Input<ApproxMultivariateTraitLikelihood>("likelihood", "trait likelihood to be logged", Validate.REQUIRED);
	public Input<Integer> posInput = new Input<Integer>("pos", "position of trait to be logged", 0);

	ApproxMultivariateTraitLikelihood likelihood;
	Tree tree;
	int pos;
	
	@Override
	public void initAndValidate() throws Exception {
		likelihood = likelihoodInput.get();
        tree =  ((Tree) (likelihood.treeInput.get()));
        pos = posInput.get();
	}

	@Override
	public int getDimension() {
		return tree.getNodeCount();
	}

	@Override
	public double getArrayValue() {
		return 0;
	}

	@Override
	public double getArrayValue(int iDim) {
		return likelihood.position[iDim][pos];
	}
}
