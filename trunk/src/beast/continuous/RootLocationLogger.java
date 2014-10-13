package beast.continuous;

import java.io.PrintStream;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;

public class RootLocationLogger extends BEASTObject implements Loggable {
	public Input<ApproxMultivariateTraitLikelihood> likelihoodInput = new Input<ApproxMultivariateTraitLikelihood>("likelihood","geographical likelihood that provides locations", Validate.REQUIRED);

	ApproxMultivariateTraitLikelihood likelihood;
	@Override
	public void initAndValidate() throws Exception {
		likelihood = likelihoodInput.get();
	}

	@Override
	public void init(PrintStream out) throws Exception {
		String name = getID();
		if (name == null) {
			name = "root";
		}
		out.append(name+".latitude\t" + name + ".longitude\t");
	}

	@Override
	public void log(int nSample, PrintStream out) {
		double [] position = likelihood.getPostion(likelihood.tree.getRoot().getNr());

		out.append(position[0] +"\t" + position[1] +"\t");
	}

	@Override
	public void close(PrintStream out) {}

}
