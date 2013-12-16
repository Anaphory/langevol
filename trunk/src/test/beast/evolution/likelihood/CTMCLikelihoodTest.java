package test.beast.evolution.likelihood;

import java.io.File;
import java.util.List;

import org.junit.Test;

import test.beast.BEASTTestCase;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.tree.Tree;
import beast.util.LogAnalyser;
import beast.util.NexusParser;

public class CTMCLikelihoodTest extends LanguageTreeLikelihoodTest {
	@Test
	public void testCTMCLikelihood() throws Exception {
		Alignment data = getBinaryAlignment();
		testCTMCLikelihood(data, "ringe-ctmc-strict-no-asc");
	}

	@Test
	public void testCTMCLikelihoodWithAscertainementCorrection() throws Exception {
		Alignment data = getBinaryAlignment();
		data = getAscertainedAlignment(data, 0, 1);
		testCTMCLikelihood(data, "ringe-ctmc-strict");
	}

	public void testCTMCLikelihood(Alignment data, String logFile) throws Exception {
		LogAnalyser analyser = new LogAnalyser(new String[] { "examples/testdata/" + logFile + ".log" }, 2000, 0);
		Double[] likelihoods = analyser.getTrace("likelihood");
		Double[] clocks = analyser.getTrace("clock.rate");
		Double[] freqs1s = analyser.getTrace("binary.frequencies1");
		Double[] freqs2s = analyser.getTrace("binary.frequencies2");

		NexusParser parser = new NexusParser();
		parser.parseFile(new File("examples/testdata/" + logFile + ".trees"));
		List<Tree> trees = parser.trees;

		assertTrue(trees.size() == likelihoods.length);

		for (int i = 0; i < likelihoods.length; i++) {
			Tree tree = trees.get(i);

			RealParameter frequencies = new RealParameter(freqs1s[i] + " " + freqs2s[i]);
			Frequencies freqs = new Frequencies();
			freqs.initByName("frequencies", frequencies);

			GeneralSubstitutionModel CTMC = new GeneralSubstitutionModel();
			CTMC.initByName("frequencies", freqs, "rates", new RealParameter("1.0 1.0"));

			SiteModel siteModel = new SiteModel();
			siteModel.initByName("mutationRate", new RealParameter(clocks[i] + ""), "gammaCategoryCount", 1, "substModel", CTMC);

			TreeLikelihood likelihood = newTreeLikelihood();
			likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);

			double fLogP = 0;
			likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "useAmbiguities", true);
			fLogP = likelihood.calculateLogP();
			assertEquals(likelihoods[i], fLogP, BEASTTestCase.PRECISION);
		}
	}

}
