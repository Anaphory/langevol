package test.beast.evolution.operators;

import org.junit.Test;

import beast.core.MCMC;
import beast.util.LogAnalyser;
import beast.util.XMLParser;

import com.sun.org.apache.xerces.internal.xni.parser.XMLParseException;

import junit.framework.TestCase;

public class DeltaExchangeScalerTest extends TestCase {
	
	String testXML = "<beast version='2.0'\n" +
"       namespace='beast.core:beast.core.util:beast.evolution.operators'>\n" +
"\n" +
"\n" +
"    <run spec='MCMC' id='mcmc' chainLength='10000000'>\n" +
"        <state>\n" +
"            <stateNode spec='parameter.RealParameter' id='freqs1' lower='0.0' upper='1.0' value='0.5 0.5'/>\n" +
"            <stateNode spec='parameter.RealParameter' id='clock1' value='1.0'/>\n" +
"\n" +
"            <stateNode spec='parameter.RealParameter' id='freqs2' lower='0.0' upper='1.0' value='0.5 0.5'/>\n" +
"            <stateNode spec='parameter.RealParameter' id='clock2' value='1.0'/>\n" +
"        </state>\n" +
"\n" +
"        <distribution spec='CompoundDistribution' id='posterior'>\n" +
"            <distribution spec='beast.math.distributions.Prior' x='@freqs1'>\n" +
"                <distr spec='beast.math.distributions.Uniform' lower='0.0' upper='1.0'/>\n" +
"            </distribution>\n" +
"            <distribution spec='beast.math.distributions.Prior' x='@clock1'>\n" +
"                <distr spec='beast.math.distributions.Gamma' alpha='2.0' beta='2.0'/>\n" +
"            </distribution>\n" +
"\n" +
"            <distribution spec='beast.math.distributions.Prior' x='@freqs2'>\n" +
"                <distr spec='beast.math.distributions.Uniform' lower='0.0' upper='1.0'/>\n" +
"            </distribution>\n" +
"            <distribution spec='beast.math.distributions.Prior' x='@clock2'>\n" +
"                <distr spec='beast.math.distributions.Gamma' alpha='2.0' beta='2.0'/>\n" +
"            </distribution>\n" +
"        </distribution>\n" +
"\n" +
"        <operator delta='0.01' id='FrequenciesExchanger' spec='DeltaExchangeOperator' weight='1' parameter='@freqs1'/>\n" +
"        <operator id='clockScaler' spec='ScaleOperator' scaleFactor='0.5' weight='1' parameter='@clock1'/>\n" +
"\n" +
"\n" +
"        <operator delta='0.01' id='FrequenciesExchanger2' spec='DeltaExchangeOperator' weight='1' parameter='@freqs2'/>\n" +
"        <operator id='clockScaler2' spec='ScaleOperator' scaleFactor='0.5' weight='1' parameter='@clock2'/>\n" +
"        <operator id='DeltaExchangeScaler2' spec='DeltaExchangeScaler' delta='0.01' weight='1' deltaexchange='@freqs2' scale='@clock2'/>\n" +
"\n" +
"\n" +
"        <logger logEvery='500' fileName='testDeltaExchangeScaler.log'>\n" +
"            <log idref='clock1'/>\n" +
"            <log idref='freqs1'/>\n" +
"            <log idref='clock2'/>\n" +
"            <log idref='freqs2'/>\n" +
"        </logger>\n" +
"\n" +
"        <logger logEvery='10000'>\n" +
"            <log idref='clock1'/>\n" +
"            <log idref='clock2'/>\n" +
"        </logger>\n" +
"    </run>\n" +
"\n" +
"</beast>";
	
	
	@Test
	public void testDeltaExchangeScaler() throws Exception {
		XMLParser parser = new XMLParser();
		MCMC mcmc = (MCMC) parser.parseFragment(testXML, true);
		//mcmc.run();
		LogAnalyser analyser = new LogAnalyser(new String[] {"testDeltaExchangeScaler.log"});
		double clock2 = analyser.getMean("clock2");
		double freqs21 = analyser.getMean("freqs21");
		assertEquals(4.0, clock2, 0.01);
		assertEquals(0.5, freqs21, 0.01);
	}

}
