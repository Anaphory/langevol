/*
 * ALSTreeLikelihood.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond,
 * Andrew Rambaut, Marc Suchard and Alexander V. Alekseyenko
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.likelihood.ThreadedTreeLikelihood;
import beast.evolution.sitemodel.SiteModelInterface;


@Description("Threaded Treelikelihood for running the Multi-State Stochastic Dollo process")
public class ThreadedALSTreeLikelihood extends ThreadedTreeLikelihood implements PartialsProvider {
    public Input<AbstractObservationProcess> op = new Input<AbstractObservationProcess>("observationprocess", "description here");

    protected AbstractObservationProcess observationProcess;

    @Override
    public void initAndValidate() {
        observationProcess = op.get();
        // ensure TreeLikelihood initialises the partials for tips
        useAmbiguitiesInput.setValue(true, this);
        super.initAndValidate();
    }

    @Override
    public double calculateLogP() {
        // Calculate the partial likelihoods
        super.calculateLogP();
        // get the frequency model
        double[] freqs = ((SiteModelInterface.Base) siteModelInput.get()).substModelInput.get().getFrequencies();
        // let the observationProcess handle the rest
        logP = observationProcess.nodePatternLikelihood(freqs, this);
        return logP;
    }

    
	@Override
	public void getNodePartials(int iNode, double[] fPartials) {
		// TODO: Fix this method
		// Was:
		// m_likelihoodCore.getNodePartials(iNode, fPartials);
		// but the property m_likelihoodCore does not exist
	}
}
