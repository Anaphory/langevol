package beast.evolution.likelihood;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.TreeInterface;

@Description("Mixture model of two likelihoods, either goverend by an indicator variable, or a mixture weight")
public class MixtureModelTreeLikelihood extends Distribution {
	public Input<List<TreeLikelihood>> likelihoodInput = new Input<List<TreeLikelihood>>("likelihood", 
			"likelihood, which one is used depends on the indicator variable", new ArrayList<>());
	public Input<IntegerParameter> indicatorInput = new Input<IntegerParameter>("indicator",
			"indicator, one for each site, which determines which likelihood is used at that site. "
			+ "At site i, likelihood number indicator[i] is used from the list");
	public Input<RealParameter> mixtureInput = new Input<RealParameter>("mixture", "weight assigned to "
			+ "first model. Either use this single mixture weight, or use indicators", 
			Validate.XOR, indicatorInput);

	
	IntegerParameter indicator;
	RealParameter mixture;
	List<TreeLikelihood> likelihoods;
	
	@Override
	public void initAndValidate() throws Exception {
		likelihoods = likelihoodInput.get();
		int nrOfSites = likelihoods.get(0).dataInput.get().getSiteCount();
		
		// make sure all likelihoods are on the same tree
		TreeInterface tree = likelihoods.get(0).treeInput.get();
		for (TreeLikelihood likelihood : likelihoods) {
			if (likelihood.treeInput.get() != tree) {
				throw new RuntimeException("All likelihoods should use the same tree");
			}
		}
		
		indicator = indicatorInput.get();
		if (indicator != null) {
			indicator.setDimension(nrOfSites);
		}

		mixture = mixtureInput.get();
	}
	
	
	@Override
	public double calculateLogP() throws Exception {
		
		
		double [][] patternLogLikelihoods = new double[indicator.getUpper()+1][];
		Alignment [] data = new Alignment[indicator.getUpper()+1];
		double [] ascertainmentCorrection = new double[indicator.getUpper()+1];
		
		for (int i = 0; i < likelihoods.size(); i++) {
			TreeLikelihood likelihood = likelihoods.get(i);
			likelihood.calculateLogP();
			if (likelihood.beagle != null) {
				likelihood.beagle.beagle.getSiteLogLikelihoods(likelihood.beagle.patternLogLikelihoods);
			}
			patternLogLikelihoods[i] =  (likelihood.beagle == null ? likelihood.patternLogLikelihoods : 
				likelihood.beagle.patternLogLikelihoods);
			data[i] = likelihood.dataInput.get();
			if (data[i].isAscertained) {
				ascertainmentCorrection[i]= data[i].getAscertainmentCorrection(patternLogLikelihoods[i]);
			}
		}

		
		logP = 0;
		if (mixture != null) {
			for (int i = 0; i < indicator.getDimension(); i++) {
				int index = indicator.getValue(i);
				int pattern = data[index].getPatternIndex(i);
				double logSiteP = patternLogLikelihoods[index][pattern] - ascertainmentCorrection[index];
				logP += logSiteP;
			}
		} else {
			double beta = mixture.getValue();
			for (int i = 0; i < data[0].getSiteCount(); i++) {
				int pattern = data[0].getPatternIndex(i);
				double logSiteP = beta * patternLogLikelihoods[0][pattern] + (1.0-beta) * patternLogLikelihoods[0][pattern]/
						(1.0 - beta * ascertainmentCorrection[0] - (1.0-beta) * ascertainmentCorrection[1]);
				logP += logSiteP;
			}
		}
		return logP;
	}
	
	
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

	@Override
	public void store() {
		super.store();
	}

	@Override
	public void restore() {
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		return super.requiresRecalculation();
	}
	
}
