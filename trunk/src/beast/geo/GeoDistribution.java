package beast.geo;

import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.substitutionmodel.Frequencies;

public class GeoDistribution extends Distribution {
	public Input<Taxon> taxonInput = new Input<Taxon>("taxon", "prior applies to this taxon");
	public Input<IntegerParameter> locationInput = new Input<IntegerParameter>("location", "locations of taxa on the graph", Validate.REQUIRED);
	public Input<TaxonSet> taxonSetInput = new Input<TaxonSet>("taxonset", "taxon comes from this set of taxa", Validate.REQUIRED);
	public Input<Boolean> onrootInput = new Input<Boolean>("onroot","whether the distribution is on the root of the tree or a leaf node", false);

	public Input<Frequencies> freqsInput = new Input<Frequencies>("frequencies", "frequencies specifying a distribution over graph nodes", Validate.REQUIRED);
	
	IntegerParameter location;
	boolean onroot;
	int taxonNr;
	Frequencies frequencies;
	
	
	@Override
	public void initAndValidate() throws Exception {
		location = locationInput.get();

		onroot = onrootInput.get();
		TaxonSet taxonSet = taxonSetInput.get();
		if (onroot) {
			taxonNr = taxonSet.getTaxonCount() * 2 - 1 - 1; 
		} else {
			Taxon taxon = taxonInput.get();
			if (taxon == null) {
				throw new RuntimeException("if onroot='false' (which it is by default) a leaf taxon must be specified");
			}
			taxonNr = taxonSet.asStringList().indexOf(taxon.getID());
			if (taxonNr < 0) {
				throw new RuntimeException("Could not find taxon \""
						+ taxon.getID() + "\" in taxonset (id=" + taxonSet.getID()
						+ "). Typo perhaps?");
			}
		}
		frequencies = freqsInput.get();
	}
	
	
	@Override
	public double calculateLogP() throws Exception {
		int i = location.getValue(taxonNr);
		logP = Math.log(frequencies.getFreqs()[i]);
		return logP;
	}

	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}

}
