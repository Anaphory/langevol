package test.beast.evolution.likelihood;

import java.io.File;
import java.util.List;

import org.junit.Test;

import test.beast.BEASTTestCase;
import beast.app.BeastMCMC;
import beast.core.Distribution;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.ALSTreeLikelihood;
import beast.evolution.likelihood.AnyTipObservationProcess;
import beast.evolution.likelihood.ThreadedALSTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SimpleDeathModel;
import beast.evolution.tree.Tree;
import beast.util.LogAnalyser;
import beast.util.NexusParser;

public class SDolloLikelihoodTest extends LanguageTreeLikelihoodTest {

	static boolean threaded  = false;
	
	@Test
	public void testSDolloLikelihood() throws Exception {
		threaded  = false;
		doTest();
	}
	
//	@Test
//	public void testThreadedSDolloLikelihood() throws Exception {
//		BeastMCMC.m_nThreads = 4;		
//		threaded  = true;
//		doTest();
//	}
	
	
	public void doTest() throws Exception {
		Alignment data = getSDolloAlignment();
		String logFile = "ringe-dollo-strict";

		LogAnalyser analyser = new LogAnalyser(new String[] { "examples/testdata/" + logFile + ".log" }, /*2000,*/ 0);
		Double[] likelihoods = analyser.getTrace("likelihood");
		Double[] cognate_loss = analyser.getTrace("cognate.loss");

		NexusParser parser = new NexusParser();
		parser.parseFile(new File("examples/testdata/" + logFile + ".trees"));
		List<Tree> trees = parser.trees;

		assertTrue(trees.size() == likelihoods.length);

		for (int i = 0; i < likelihoods.length; i++) {
			Tree tree = trees.get(i);

			RealParameter frequencies = new RealParameter("1 0");
			Frequencies freqs = new Frequencies();
			freqs.initByName("frequencies", frequencies);

			RealParameter deathrate = new RealParameter(cognate_loss[i] + "");
			SimpleDeathModel SDollo = new SimpleDeathModel();
			SDollo.initByName("deathprob", deathrate, "frequencies", freqs);

			SiteModel siteModel = new SiteModel();
			siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", SDollo);

			AnyTipObservationProcess observationProcess = new AnyTipObservationProcess();
			observationProcess.initByName("tree", tree, "siteModel", siteModel, "data", data,
					"integrateGainRate", "true", "mu", deathrate);
			
			System.setProperty("java.only", "true");
			
			Distribution likelihood = getTreeLikelihood();
			likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "observationprocess", observationProcess);

			double fLogP = 0;
			likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "useAmbiguities", true);
			fLogP = likelihood.calculateLogP();
			assertEquals(likelihoods[i], fLogP, BEASTTestCase.PRECISION);
		}
	}

	private Distribution getTreeLikelihood() {
		if (threaded) {
			return new ThreadedALSTreeLikelihood();
		} else {
			// System.setProperty("java.only", "false"); subst-model does not support BEAGLE
			System.setProperty("java.only", "false");
			return new ALSTreeLikelihood();
		}
	}

}
