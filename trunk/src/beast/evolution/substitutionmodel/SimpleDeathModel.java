package beast.evolution.substitutionmodel;


import java.lang.reflect.Constructor;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;

@Description("Simple Death substitution model, can be used as Stochastic Dollo model.")
public class SimpleDeathModel extends SubstitutionModel.Base {

    public Input<RealParameter> delParameter = new Input<RealParameter>("deathprob", "rate of death, used to calculate death probability", Validate.REQUIRED);

    /**
     * transition matrix for live states *
     */
    protected double[] trMatrix;
    /**
     * number of states *
     */
    int nrOfStates;
    boolean updateMatrix = true;

    protected EigenSystem eigenSystem;

    protected EigenDecomposition eigenDecomposition;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        double[] freqs = getFrequencies();
        nrOfStates = freqs.length;
        trMatrix = new double[(nrOfStates - 1) * (nrOfStates - 1)];
        
        eigenSystem = new DefaultEigenSystem(2);
    }
    
    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        synchronized (this) {
            if (updateMatrix) {
                eigenDecomposition = eigenSystem.decomposeMatrix(getRateMatrix());
                updateMatrix = false;
            }
        }
        return eigenDecomposition;
    }
    
	protected double[][] getRateMatrix() {
		double [][] rateMatrix = new double[2][2];
		rateMatrix[0][0] = - delParameter.get().getValue();
		rateMatrix[0][1] = delParameter.get().getValue();
		rateMatrix[1][0] = 0;
		rateMatrix[1][1] = 0;
		return rateMatrix;
	}
    
    @Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
        double distance = (fStartTime - fEndTime) * fRate;
        int i, j;
        // assuming that expected number of changes in CTMCModel is 1 per unit time
        // we are contributing s*deathRate number of changes per unit of time
        double deathProb = Math.exp(-distance * delParameter.get().getValue());
        trMatrix[0] = 1.0;

        for (i = 0; i < nrOfStates - 1; ++i) {
            for (j = 0; j < nrOfStates - 1; j++) {
                matrix[i * (nrOfStates) + j] = trMatrix[i * (nrOfStates - 1) + j] * deathProb;
            }
            matrix[i * (nrOfStates) + j] = (1.0 - deathProb);
        }

        for (j = 0; j < nrOfStates - 1; ++j) {
            matrix[nrOfStates * (nrOfStates - 1) + j] = 0.0;
        }

        matrix[nrOfStates * nrOfStates - 1] = 1.0;
        
        
//        eigenDecomposition = getEigenDecomposition(node);
//        double temp;
//        double[] iexp = new double[nrOfStates * nrOfStates];
//        // Eigen vectors
//        double[] Evec = eigenDecomposition.getEigenVectors();
//        // inverse Eigen vectors
//        double[] Ievc = eigenDecomposition.getInverseEigenVectors();
//        // Eigen values
//        double[] Eval = eigenDecomposition.getEigenValues();
//        for (i = 0; i < nrOfStates; i++) {
//            temp = Math.exp(distance * Eval[i]);
//            for (j = 0; j < nrOfStates; j++) {
//                iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
//            }
//        }
//
//        int u = 0;
//        for (i = 0; i < nrOfStates; i++) {
//            for (j = 0; j < nrOfStates; j++) {
//                temp = 0.0;
//                for (int k = 0; k < nrOfStates; k++) {
//                    temp += Evec[i * nrOfStates + k] * iexp[k * nrOfStates + j];
//                }
//
//                matrix[u] = Math.abs(temp);
//                u++;
//            }
//        }

    } // getTransitionProbabilities

    /**
     * CalculationNode implementation *
     */
    @Override
    protected boolean requiresRecalculation() {
        // we only get here if delParameter or mutationRate is dirty
    	updateMatrix = true;
        return true;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) throws Exception {
   		return dataType.getStateCount() == 2;
    }

} // class SimpleMutationDeathModel
