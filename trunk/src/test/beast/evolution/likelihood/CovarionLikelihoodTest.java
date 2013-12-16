package test.beast.evolution.likelihood;

import java.io.File;
import java.util.List;

import org.junit.Test;

import test.beast.BEASTTestCase;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.BinaryCovarion;
import beast.evolution.tree.Tree;
import beast.util.LogAnalyser;
import beast.util.NexusParser;

public class CovarionLikelihoodTest extends LanguageTreeLikelihoodTest {

	 @Test
	 public void testCovarionLikelihood() throws Exception {
		 Alignment data = new Alignment();
		 data.initByName("sequence", oldnorse, "sequence", avestan, "sequence",
		 gothic, "sequence", luvian, "sequence", oldpersian, "sequence", vedic,
		 "sequence", umbrian, "sequence", oldhighgerman, "sequence", oldprussian,
		 "sequence", latin, "sequence", welsh, "sequence",
		 lithuanian, "sequence", oldenglish, "sequence", classicalarmenian,
		 "sequence", hittite, "sequence", oldirish, "sequence", albanian,
		 "sequence", oldchurchslavonic, "sequence", ancientgreek, "sequence",
		 lycian, "sequence", latvian, "sequence", tocharian_b,
		 "sequence", tocharian_a, "sequence", oscan, "dataType",
		 "twoStateCovarion", "strip", true);
		
		 testCovarionLikelihood(data, "ringe-covarion-strict-no-asc");
	 }
	
	 @Test
	 public void testCovarionLikelihoodWithAscertainementCorrection() throws
	 Exception {
		 Alignment data = new Alignment();
		 data.initByName("sequence", oldnorse, "sequence", avestan, "sequence",
		 gothic, "sequence", luvian, "sequence", oldpersian, "sequence", vedic,
		 "sequence", umbrian, "sequence", oldhighgerman, "sequence", oldprussian,
		 "sequence", latin, "sequence", welsh, "sequence",
		 lithuanian, "sequence", oldenglish, "sequence", classicalarmenian,
		 "sequence", hittite, "sequence", oldirish, "sequence", albanian,
		 "sequence", oldchurchslavonic, "sequence", ancientgreek, "sequence",
		 lycian, "sequence", latvian, "sequence", tocharian_b,
		 "sequence", tocharian_a, "sequence", oscan, "dataType",
		 "twoStateCovarion", "strip", true);
		 data = getAscertainedAlignment(data, 0, 1);
		 testCovarionLikelihood(data, "ringe-covarion-strict");
	 }

	public void testCovarionLikelihood(Alignment data, String logFile) throws Exception {
		LogAnalyser analyser = new LogAnalyser(new String[] { "examples/testdata/" + logFile + ".log" }, 2000, 0);
		Double[] likelihoods = analyser.getTrace("likelihood");
		Double[] clocks = analyser.getTrace("clock.rate");
		Double[] freqs1s = analyser.getTrace("frequencies1");
		Double[] freqs2s = analyser.getTrace("frequencies2");
		Double[] bcov_alpha = analyser.getTrace("bcov.alpha");
		Double[] bcov_s = analyser.getTrace("bcov.s");

		NexusParser parser = new NexusParser();
		parser.parseFile(new File("examples/testdata/" + logFile + ".trees"));
		List<Tree> trees = parser.trees;

		assertTrue(trees.size() == likelihoods.length);

		for (int i = 0; i < likelihoods.length; i++) {
			Tree tree = trees.get(i);

			RealParameter alpha = new RealParameter(bcov_alpha[i] + "");
			RealParameter switchRate = new RealParameter(bcov_s[i] + "");
			RealParameter frequencies = new RealParameter(freqs1s[i] + " " + freqs2s[i]);
			RealParameter hfrequencies = new RealParameter("0.5 0.5");
			BinaryCovarion covarion = new BinaryCovarion();
			covarion.initByName("alpha", alpha, "switchRate", switchRate, "vfrequencies", frequencies, "hfrequencies", hfrequencies);

			SiteModel siteModel = new SiteModel();
			siteModel.initByName("mutationRate", new RealParameter(clocks[i] + ""), "gammaCategoryCount", 1, "substModel", covarion);

			TreeLikelihood likelihood = newTreeLikelihood();
			likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);

			double fLogP = 0;
			likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel);
			fLogP = likelihood.calculateLogP();
			String s = tree.getRoot().toNewick();
			assertEquals(likelihoods[i], fLogP, BEASTTestCase.PRECISION);
		}
	}

}
