package beast.geo;

import java.util.Arrays;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.EigenSystem;


/**
 * Matrix represented by flat arrays based on Apache commons implementation of org.apache.commons.math3.linear.EigenDecomposition
**/
public class FlatMatrixEigenSystem2 implements EigenSystem {
    /** Internally used epsilon criteria. */
    private static final double EPSILON = 1e-12;
    /** Maximum number of iterations accepted in the implicit QL transformation */
    private byte maxIter = 30;
    /** Main diagonal of the tridiagonal matrix. */
    private double[] main;
    /** Secondary diagonal of the tridiagonal matrix. */
    private double[] secondary;
    /**
     * Transformer to tridiagonal (may be null if matrix is already
     * tridiagonal).
     */
    private TriDiagonalTransformer transformer;
    /** Real part of the realEigenvalues. */
    private double[] realEigenvalues;
    /** Imaginary part of the realEigenvalues. */
    private double[] imagEigenvalues;
    /** Eigenvectors. */
    //private ArrayRealVector[] eigenvectors;
    private double[] flateigenvectors;

    @Override
	public EigenDecomposition decomposeMatrix(double[][] matrix) {
		return null;
	}
	
	// dimension of matrix
	int n;
	
    public EigenDecomposition decomposeMatrix(double[] matrix, boolean isSymmetric) throws MathArithmeticException  {
    	System.out.println("Start EigneDecomposition");
        long start = System.currentTimeMillis();
        long start0 = start;
    	n = (int) Math.sqrt(matrix.length);
    	if (n*n != matrix.length) {
    		throw new RuntimeException("expected a square matrix");
    	}
    	
        //final double symTol = 10 * matrix.length * Precision.EPSILON;

        double [] Eval = new double[n];
        double [] Evec = new double[(n*n)];
        double [] Ievc = new double[(n*n)];

        if (isSymmetric) {
            transformToTridiagonal(matrix);
            long end = System.currentTimeMillis();System.out.println("SymStep 1 " + ((end-start)/100)/10.0+" seconds");start=end;
            
            findEigenVectors(transformer.getQ().getData());
            //for (int i  = 0; i < n; i++) {
            //	System.arraycopy(eigenvectors[i].toArray(), 0, Evec, i * n, n);
            //}
            System.arraycopy(flateigenvectors, 0, Evec, 0, n*n);
            end = System.currentTimeMillis();System.out.println("SymStep 2 " + ((end-start)/100)/10.0+" seconds");start=end;
            //findEigenVectors(transformer.getQ().transpose().getData());
            //for (int i  = 0; i < n; i++) {
            //	System.arraycopy(eigenvectors[i].toArray(), 0, Ievc, i * n, n);
            //}
        } else {
            final SchurTransformer t = transformToSchur(matrix);
            long end = System.currentTimeMillis();System.out.println("AStep 1 " + ((end-start)/100)/10.0+" seconds");start=end;
            findEigenVectorsFromSchur(t);
            //for (int i  = 0; i < n; i++) {
            //	System.arraycopy(eigenvectors[i].toArray(), 0, Evec, i * n, n);
            //}
            //System.arraycopy(flateigenvectors, 0, Evec, 0, n*n);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                	Evec[i*n+j] =  flateigenvectors[i + j*n];
                }
            }
            end = System.currentTimeMillis();System.out.println("AStep 2 " + ((end-start)/100)/10.0+" seconds");start=end;
        }
        
        RealMatrix m = new Array2DRowRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
            	m.setEntry(i,  j, flateigenvectors[i + j*n]);
            }
        }
        RealMatrix inverse = new LUDecomposition(m).getSolver().getInverse();
        for (int i  = 0; i < n; i++) {
        	System.arraycopy(inverse.getColumn(i), 0, Ievc, i * n, n);
        }
    	

        long end = System.currentTimeMillis();System.out.println("Step 3 " + ((end-start)/100)/10.0+" seconds");start=end;
        System.out.println("Total EigenDecomposition time " + ((end-start0)/100)/10.0+" seconds");
        return new EigenDecomposition(Evec, Ievc, realEigenvalues, imagEigenvalues);
    }

    /**
     * Transforms the matrix to tridiagonal form.
     *
     * @param matrix Matrix to transform.
     */
    private void transformToTridiagonal(final double [] matrix) {
        // transform the matrix to tridiagonal
        transformer = new TriDiagonalTransformer(matrix);
        main = transformer.getMainDiagonalRef();
        secondary = transformer.getSecondaryDiagonalRef();
    }

    /**
     * Find eigenvalues and eigenvectors (Dubrulle et al., 1971)
     *
     * @param householderMatrix Householder matrix of the transformation
     * to tridiagonal form.
     */
    private void findEigenVectors(final double[][] householderMatrix) {
        final double[][]z = householderMatrix.clone();
        //final int n = main.length;
        realEigenvalues = new double[n];
        imagEigenvalues = new double[n];
        final double[] e = new double[n];
        for (int i = 0; i < n - 1; i++) {
            realEigenvalues[i] = main[i];
            e[i] = secondary[i];
        }
        realEigenvalues[(n - 1)] = main[(n - 1)];
        e[(n - 1)] = 0;

        // Determine the largest main and secondary value in absolute term.
        double maxAbsoluteValue = 0;
        for (int i = 0; i < n; i++) {
            if (FastMath.abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = FastMath.abs(realEigenvalues[i]);
            }
            if (FastMath.abs(e[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = FastMath.abs(e[i]);
            }
        }
        // Make null any main and secondary value too small to be significant
        if (maxAbsoluteValue != 0) {
            for (int i=0; i < n; i++) {
                if (FastMath.abs(realEigenvalues[i]) <= Precision.EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0;
                }
                if (FastMath.abs(e[i]) <= Precision.EPSILON * maxAbsoluteValue) {
                    e[i]=0;
                }
            }
        }

        for (int j = 0; j < n; j++) {
            int its = 0;
            int m;
            do {
                for (m = j; m < n - 1; m++) {
                    double delta = FastMath.abs(realEigenvalues[(m)]) +
                        FastMath.abs(realEigenvalues[(m + 1)]);
                    if (FastMath.abs(e[(m)]) + delta == delta) {
                        break;
                    }
                }
                if (m != j) {
                    if (its == maxIter) {
                        throw new MaxCountExceededException(LocalizedFormats.CONVERGENCE_FAILED,
                                                            maxIter);
                    }
                    its++;
                    double q = (realEigenvalues[(j + 1)] - realEigenvalues[j]) / (2 * e[j]);
                    double t = FastMath.sqrt(1 + q * q);
                    if (q < 0.0) {
                        q = realEigenvalues[(m)] - realEigenvalues[j] + e[j] / (q - t);
                    } else {
                        q = realEigenvalues[(m)] - realEigenvalues[j] + e[j] / (q + t);
                    }
                    double u = 0.0;
                    double s = 1.0;
                    double c = 1.0;
                    int i;
                    for (i = m - 1; i >= j; i--) {
                        double p = s * e[i];
                        double h = c * e[i];
                        if (FastMath.abs(p) >= FastMath.abs(q)) {
                            c = q / p;
                            t = FastMath.sqrt(c * c + 1.0);
                            e[(i + 1)] = p * t;
                            s = 1.0 / t;
                            c *= s;
                        } else {
                            s = p / q;
                            t = FastMath.sqrt(s * s + 1.0);
                            e[(i + 1)] = q * t;
                            c = 1.0 / t;
                            s *= c;
                        }
                        if (e[(i + 1)] == 0.0) {
                            realEigenvalues[(i + 1)] -= u;
                            e[(m)] = 0.0;
                            break;
                        }
                        q = realEigenvalues[(i + 1)] - u;
                        t = (realEigenvalues[i] - q) * s + 2.0 * c * h;
                        u = s * t;
                        realEigenvalues[(i + 1)] = q + u;
                        q = c * t - h;
                        for (int ia = 0; ia < n; ia++) {
                            p = z[ia][(i + 1)];
                            z[ia][(i + 1)] = s * z[ia][i] + c * p;
                            z[ia][i] = c * z[ia][i] - s * p;
                        }
                    }
                    if (t == 0.0 && i >= j) {
                        continue;
                    }
                    realEigenvalues[j] -= u;
                    e[j] = q;
                    e[(m)] = 0.0;
                }
            } while (m != j);
        }

        //Sort the eigen values (and vectors) in increase order
        for (int i = 0; i < n; i++) {
            int k = i;
            double p = realEigenvalues[i];
            for (int j = i + 1; j < n; j++) {
                if (realEigenvalues[j] > p) {
                    k = j;
                    p = realEigenvalues[j];
                }
            }
            if (k != i) {
                realEigenvalues[k] = realEigenvalues[i];
                realEigenvalues[i] = p;
                for (int j = 0; j < n; j++) {
                    p = z[j][i];
                    z[j][i] = z[j][k];
                    z[j][k] = p;
                }
            }
        }

        // Determine the largest eigen value in absolute term.
        maxAbsoluteValue = 0;
        for (int i = 0; i < n; i++) {
            if (FastMath.abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue=FastMath.abs(realEigenvalues[i]);
            }
        }
        // Make null any eigen value too small to be significant
        if (maxAbsoluteValue != 0.0) {
            for (int i=0; i < n; i++) {
                if (FastMath.abs(realEigenvalues[i]) < Precision.EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0;
                }
            }
        }
//        eigenvectors = new ArrayRealVector[n];
//        final double[] tmp = new double[n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                tmp[j] = z[j][i];
//            }
//            eigenvectors[i] = new ArrayRealVector(tmp);
//        }
    }
    
    /**
     * Find eigenvectors from a matrix transformed to Schur form.
     *
     * @param schur the schur transformation of the matrix
     * @throws MathArithmeticException if the Schur form has a norm of zero
     */
    private void findEigenVectorsFromSchur(final SchurTransformer schur)
        throws MathArithmeticException {
        final double[] matrixT = schur.getTraw();
        final double[] matrixP = schur.getPraw();


        // compute matrix norm
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
           for (int j = FastMath.max(i - 1, 0); j < n; j++) {
               norm += FastMath.abs(matrixT[i*n+j]);
           }
        }

        // we can not handle a matrix with zero norm
        if (Precision.equals(norm, 0.0, EPSILON)) {
           throw new MathArithmeticException(LocalizedFormats.ZERO_NORM);
        }

        // Backsubstitute to find vectors of upper triangular form

        double r = 0.0;
        double s = 0.0;
        double z = 0.0;

        for (int idx = n - 1; idx >= 0; idx--) {
            double p = realEigenvalues[(idx)];
            double q = imagEigenvalues[(idx)];

            if (Precision.equals(q, 0.0)) {
                // Real vector
                int l = idx;
                matrixT[(idx)*n+(idx)] = 1.0;
                for (int i = idx - 1; i >= 0; i--) {
                    double w = matrixT[i*n+i] - p;
                    r = 0.0;
                    for (int j = l; j <= idx; j++) {
                        r += matrixT[i*n+j] * matrixT[j*n+(idx)];
                    }
                    if (Precision.compareTo(imagEigenvalues[i], 0.0, EPSILON) < 0) {
                        z = w;
                        s = r;
                    } else {
                        l = i;
                        if (Precision.equals(imagEigenvalues[i], 0.0)) {
                            if (w != 0.0) {
                                matrixT[i*n+(idx)] = -r / w;
                            } else {
                                matrixT[i*n+(idx)] = -r / (Precision.EPSILON * norm);
                            }
                        } else {
                            // Solve real equations
                            double x = matrixT[i*n+(i + 1)];
                            double y = matrixT[(i + 1)*n+i];
                            q = (realEigenvalues[i] - p) * (realEigenvalues[i] - p) +
                                imagEigenvalues[i] * imagEigenvalues[i];
                            double t = (x * s - z * r) / q;
                            matrixT[i*n+(idx)] = t;
                            if (FastMath.abs(x) > FastMath.abs(z)) {
                                matrixT[(i + 1)*n+(idx)] = (-r - w * t) / x;
                            } else {
                                matrixT[(i + 1)*n+(idx)] = (-s - y * t) / z;
                            }
                        }

                        // Overflow control
                        double t = FastMath.abs(matrixT[i*n+(idx)]);
                        if ((Precision.EPSILON * t) * t > 1) {
                            for (int j = i; j <= idx; j++) {
                                matrixT[j*n+(idx)] /= t;
                            }
                        }
                    }
                }
            } else if (q < 0.0) {
                // Complex vector
                int l = idx - 1;

                // Last vector component imaginary so matrix is triangular
                if (FastMath.abs(matrixT[(idx)*n+(idx - 1)]) > FastMath.abs(matrixT[(idx - 1)*n+(idx)])) {
                    matrixT[(idx - 1)*n+(idx - 1)] = q / matrixT[(idx)*n+(idx - 1)];
                    matrixT[(idx - 1)*n+(idx)]     = -(matrixT[(idx)*n+(idx)] - p) / matrixT[(idx)*n+(idx - 1)];
                } else {
                    final Complex result = cdiv(0.0, -matrixT[(idx - 1)*n+(idx)],
                                                matrixT[(idx - 1)*n+(idx - 1)] - p, q);
                    matrixT[(idx - 1)*n+(idx - 1)] = result.getReal();
                    matrixT[(idx - 1)*n+(idx)]     = result.getImaginary();
                }

                matrixT[(idx)*n+(idx - 1)] = 0.0;
                matrixT[(idx)*n+(idx)]     = 1.0;

                for (int i = idx - 2; i >= 0; i--) {
                    double ra = 0.0;
                    double sa = 0.0;
                    for (int j = l; j <= idx; j++) {
                        ra += matrixT[i*n+j] * matrixT[j*n+(idx - 1)];
                        sa += matrixT[i*n+j] * matrixT[j*n+(idx)];
                    }
                    double w = matrixT[i*n+i] - p;

                    if (Precision.compareTo(imagEigenvalues[i], 0.0, EPSILON) < 0) {
                        z = w;
                        r = ra;
                        s = sa;
                    } else {
                        l = i;
                        if (Precision.equals(imagEigenvalues[i], 0.0)) {
                            final Complex c = cdiv(-ra, -sa, w, q);
                            matrixT[i*n+(idx - 1)] = c.getReal();
                            matrixT[i*n+(idx)] = c.getImaginary();
                        } else {
                            // Solve complex equations
                            double x = matrixT[i*n+(i + 1)];
                            double y = matrixT[(i + 1)*n+i];
                            double vr = (realEigenvalues[i] - p) * (realEigenvalues[i] - p) +
                                        imagEigenvalues[i] * imagEigenvalues[i] - q * q;
                            final double vi = (realEigenvalues[i] - p) * 2.0 * q;
                            if (Precision.equals(vr, 0.0) && Precision.equals(vi, 0.0)) {
                                vr = Precision.EPSILON * norm *
                                     (FastMath.abs(w) + FastMath.abs(q) + FastMath.abs(x) +
                                      FastMath.abs(y) + FastMath.abs(z));
                            }
                            final Complex c     = cdiv(x * r - z * ra + q * sa,
                                                       x * s - z * sa - q * ra, vr, vi);
                            matrixT[i*n+(idx - 1)] = c.getReal();
                            matrixT[i*n+(idx)]     = c.getImaginary();

                            if (FastMath.abs(x) > (FastMath.abs(z) + FastMath.abs(q))) {
                                matrixT[(i + 1)*n+(idx - 1)] = (-ra - w * matrixT[i*n+(idx - 1)] +
                                                           q * matrixT[i*n+(idx)]) / x;
                                matrixT[(i + 1)*n+(idx)]     = (-sa - w * matrixT[i*n+(idx)] -
                                                           q * matrixT[i*n+(idx - 1)]) / x;
                            } else {
                                final Complex c2        = cdiv(-r - y * matrixT[i*n+(idx - 1)],
                                                               -s - y * matrixT[i*n+(idx)], z, q);
                                matrixT[(i + 1)*n+(idx - 1)] = c2.getReal();
                                matrixT[(i + 1)*n+(idx)]     = c2.getImaginary();
                            }
                        }

                        // Overflow control
                        double t = FastMath.max(FastMath.abs(matrixT[i*n+(idx - 1)]),
                                                FastMath.abs(matrixT[i*n+(idx)]));
                        if ((Precision.EPSILON * t) * t > 1) {
                            for (int j = i; j <= idx; j++) {
                                matrixT[j*n+(idx - 1)] /= t;
                                matrixT[j*n+(idx)] /= t;
                            }
                        }
                    }
                }
            }
        }

        // Back transformation to get eigenvectors of original matrix
        for (int j = n - 1; j >= 0; j--) {
            for (int i = 0; i <= n - 1; i++) {
                z = 0.0;
                for (int k = 0; k <= FastMath.min(j, n - 1); k++) {
                    z += matrixP[i+n*k] * matrixT[k*n+j];
                }
                matrixP[i+n*j] = z;
            }
        }

        
        flateigenvectors = matrixP;
        
//        eigenvectors = new ArrayRealVector[n];
//        final double[] tmp = new double[n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                tmp[j] = matrixP[j+n*i];
//            }
//            eigenvectors[i] = new ArrayRealVector(tmp);
//        }
    }

    /**
     * Performs a division of two complex numbers.
     *
     * @param xr real part of the first number
     * @param xi imaginary part of the first number
     * @param yr real part of the second number
     * @param yi imaginary part of the second number
     * @return result of the complex division
     */
    private Complex cdiv(final double xr, final double xi,
                         final double yr, final double yi) {
        return new Complex(xr, xi).divide(new Complex(yr, yi));
    }
    
    class TriDiagonalTransformer {
        /** Householder vectors. */
        private final double householderVectors[][];
        /** Main diagonal. */
        private final double[] main;
        /** Secondary diagonal. */
        private final double[] secondary;
        /** Cached value of Q. */
        private RealMatrix cachedQ;
        /** Cached value of Qt. */
        private RealMatrix cachedQt;
        /** Cached value of T. */
        private RealMatrix cachedT;

        /**
         * Build the transformation to tridiagonal shape of a symmetrical matrix.
         * <p>The specified matrix is assumed to be symmetrical without any check.
         * Only the upper triangular part of the matrix is used.</p>
         *
         * @param matrix Symmetrical matrix to transform.
         * @throws NonSquareMatrixException if the matrix is not square.
         */
        public TriDiagonalTransformer(double[] matrix) {
//            if (!matrix.isSquare()) {
//                throw new NonSquareMatrixException(matrix.getRowDimension(),
//                                                   matrix.getColumnDimension());
//            }

            final int m = n;
            householderVectors= new double[n][n];
            for (int i = 0; i < n; i++) {
            	System.arraycopy(matrix, i*n, householderVectors[i], 0, n);
            }
            main      = new double[(m)];
            secondary = new double[(m - 1)];
            cachedQ   = null;
            cachedQt  = null;
            cachedT   = null;

            // transform matrix
            transform();
        }

        /**
         * Returns the matrix Q of the transform.
         * <p>Q is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
         * @return the Q matrix
         */
        public RealMatrix getQ() {
            if (cachedQ == null) {
                cachedQ = getQT().transpose();
            }
            return cachedQ;
        }

        /**
         * Returns the transpose of the matrix Q of the transform.
         * <p>Q is an orthogonal matrix, i.e. its transpose is also its inverse.</p>
         * @return the Q matrix
         */
        public RealMatrix getQT() {
            if (cachedQt == null) {
                final int m = n;//householderVectors.length;
                double[][] qta = new double[(m)][(m)];

                // build up first part of the matrix by applying Householder transforms
                for (int k = m - 1; k >= 1; --k) {
                    final double[] hK = householderVectors[(k - 1)];
                    qta[k][k] = 1;
                    if (hK[k] != 0.0) {
                        final double inv = 1.0 / (secondary[(k - 1)] * hK[k]);
                        double beta = 1.0 / secondary[(k - 1)];
                        qta[k][k] = 1 + beta * hK[k];
                        for (int i = k + 1; i < m; ++i) {
                            qta[k][i] = beta * hK[i];
                        }
                        for (int j = k + 1; j < m; ++j) {
                            beta = 0;
                            for (int i = k + 1; i < m; ++i) {
                                beta += qta[j][i] * hK[i];
                            }
                            beta *= inv;
                            qta[j][k] = beta * hK[k];
                            for (int i = k + 1; i < m; ++i) {
                                qta[j][i] += beta * hK[i];
                            }
                        }
                    }
                }
                qta[(0)][(0)] = 1;
                cachedQt = MatrixUtils.createRealMatrix(qta);
            }

            // return the cached matrix
            return cachedQt;
        }

        /**
         * Returns the tridiagonal matrix T of the transform.
         * @return the T matrix
         */
        public RealMatrix getT() {
            if (cachedT == null) {
                final int m = n;//main.length;
                double[][] ta = new double[(m)][(m)];
                for (int i = 0; i < m; ++i) {
                    ta[i][i] = main[i];
                    if (i > 0) {
                        ta[i][(i - 1)] = secondary[(i - 1)];
                    }
                    if (i < n - 1) {
                        ta[i][(i + 1)] = secondary[i];
                    }
                }
                cachedT = MatrixUtils.createRealMatrix(ta);
            }

            // return the cached matrix
            return cachedT;
        }

        /**
         * Get the Householder vectors of the transform.
         * <p>Note that since this class is only intended for internal use,
         * it returns directly a reference to its internal arrays, not a copy.</p>
         * @return the main diagonal elements of the B matrix
         */
        double[][] getHouseholderVectorsRef() {
            return householderVectors;
        }

        /**
         * Get the main diagonal elements of the matrix T of the transform.
         * <p>Note that since this class is only intended for internal use,
         * it returns directly a reference to its internal arrays, not a copy.</p>
         * @return the main diagonal elements of the T matrix
         */
        double[] getMainDiagonalRef() {
            return main;
        }

        /**
         * Get the secondary diagonal elements of the matrix T of the transform.
         * <p>Note that since this class is only intended for internal use,
         * it returns directly a reference to its internal arrays, not a copy.</p>
         * @return the secondary diagonal elements of the T matrix
         */
        double[] getSecondaryDiagonalRef() {
            return secondary;
        }

        /**
         * Transform original matrix to tridiagonal form.
         * <p>Transformation is done using Householder transforms.</p>
         */
        private void transform() {
            final int m = n;//householderVectors.length;
            final double[] z = new double[(m)];
            for (int k = 0; k < m - 1; k++) {

                //zero-out a row and a column simultaneously
                final double[] hK = householderVectors[k];
                main[k] = hK[k];
                double xNormSqr = 0;
                for (int j = k + 1; j < m; ++j) {
                    final double c = hK[j];
                    xNormSqr += c * c;
                }
                final double a = (hK[(k + 1)] > 0) ? -FastMath.sqrt(xNormSqr) : FastMath.sqrt(xNormSqr);
                secondary[k] = a;
                if (a != 0.0) {
                    // apply Householder transform from left and right simultaneously

                    hK[(k + 1)] -= a;
                    final double beta = -1 / (a * hK[(k + 1)]);

                    // compute a = beta A v, where v is the Householder vector
                    // this loop is written in such a way
                    //   1) only the upper triangular part of the matrix is accessed
                    //   2) access is cache-friendly for a matrix stored in rows
                    Arrays.fill(z, k + 1, m, 0);
                    for (int i = k + 1; i < m; ++i) {
                        final double[] hI = householderVectors[i];
                        final double hKI = hK[i];
                        double zI = hI[i] * hKI;
                        for (int j = i + 1; j < m; ++j) {
                            final double hIJ = hI[j];
                            zI   += hIJ * hK[j];
                            z[j] += hIJ * hKI;
                        }
                        z[i] = beta * (z[i] + zI);
                    }

                    // compute gamma = beta vT z / 2
                    double gamma = 0;
                    for (int i = k + 1; i < m; ++i) {
                        gamma += z[i] * hK[i];
                    }
                    gamma *= beta / 2;

                    // compute z = z - gamma v
                    for (int i = k + 1; i < m; ++i) {
                        z[i] -= gamma * hK[i];
                    }

                    // update matrix: A = A - v zT - z vT
                    // only the upper triangular part of the matrix is updated
                    for (int i = k + 1; i < m; ++i) {
                        final double[] hI = householderVectors[i];
                        for (int j = i; j < m; ++j) {
                            hI[j] -= hK[i] * z[j] + z[i] * hK[j];
                        }
                    }
                }
            }
            main[(m - 1)] = householderVectors[(m - 1)][(m - 1)];
        }
    }

    /**
     * Transforms the matrix to Schur form and calculates the eigenvalues.
     *
     * @param matrix Matrix to transform.
     * @return the {@link SchurTransformer Shur transform} for this matrix
     */
    private SchurTransformer transformToSchur(final double [] matrix) {
        final SchurTransformer schurTransform = new SchurTransformer(matrix);
        final double[][] matT = schurTransform.getT().getData();

        realEigenvalues = new double[n];
        imagEigenvalues = new double[n];

        for (int i = 0; i < n; i++) {
            if (i == (n - 1) ||
                Precision.equals(matT[(i + 1)][i], 0.0, EPSILON)) {
                realEigenvalues[i] = matT[i][i];
            } else {
                final double x = matT[(i + 1)][(i + 1)];
                final double p = 0.5 * (matT[i][i] - x);
                final double z = FastMath.sqrt(FastMath.abs(p * p + matT[(i + 1)][i] * matT[i][(i + 1)]));
                realEigenvalues[i] = x + p;
                imagEigenvalues[i] = z;
                realEigenvalues[(i + 1)] = x + p;
                imagEigenvalues[(i + 1)] = -z;
                i++;
            }
        }
        return schurTransform;
    }


    /**
     * Class transforming a general real matrix to Schur form.
     * <p>A m &times; m matrix A can be written as the product of three matrices: A = P
     * &times; T &times; P<sup>T</sup> with P an orthogonal matrix and T an quasi-triangular
     * matrix. Both P and T are m &times; m matrices.</p>
     * <p>Transformation to Schur form is often not a goal by itself, but it is an
     * intermediate step in more general decomposition algorithms like
     * {@link EigenDecomposition eigen decomposition}. This class is therefore
     * intended for internal use by the library and is not public. As a consequence
     * of this explicitly limited scope, many methods directly returns references to
     * internal arrays, not copies.</p>
     * <p>This class is based on the method hqr2 in class EigenvalueDecomposition
     * from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</p>
     *
     * @see <a href="http://mathworld.wolfram.com/SchurDecomposition.html">Schur Decomposition - MathWorld</a>
     * @see <a href="http://en.wikipedia.org/wiki/Schur_decomposition">Schur Decomposition - Wikipedia</a>
     * @see <a href="http://en.wikipedia.org/wiki/Householder_transformation">Householder Transformations</a>
     * @version $Id: SchurTransformer.java 1538368 2013-11-03 13:57:37Z erans $
     * @since 3.1
     */
    class SchurTransformer {
        /** Maximum allowed iterations for convergence of the transformation. */
        private static final int MAX_ITERATIONS = 100;

        /** P matrix. */
        private final double matrixP[];
        /** T matrix. */
        private final double matrixT[];
        /** Cached value of P. */
        private RealMatrix cachedP;
        /** Cached value of T. */
        private RealMatrix cachedT;
        /** Cached value of PT. */
        private RealMatrix cachedPt;

        /** Epsilon criteria taken from JAMA code (originally was 2^-52). */
        private final double epsilon = Precision.EPSILON;

        /**
         * Build the transformation to Schur form of a general real matrix.
         *
         * @param matrix matrix to transform
         * @throws NonSquareMatrixException if the matrix is not square
         */
        public SchurTransformer(final double [] matrix) {
//            if (!matrix.isSquare()) {
//                throw new NonSquareMatrixException(matrix.getRowDimension(),
//                                                   matrix.getColumnDimension());
//            }
			long start = System.currentTimeMillis();

            HessenbergTransformer transformer = new HessenbergTransformer(matrix);
			long end = System.currentTimeMillis();
			System.out.println("ST 1 Done in " + ((end-start)/100)/10.0+" seconds");
			start = end;

			matrixT = transformer.getH().getDataR();
            matrixP = transformer.getP().getDataR();
            cachedT = null;
            cachedP = null;
            cachedPt = null;

			end = System.currentTimeMillis();
			System.out.println("ST 1 Done in " + ((end-start)/100)/10.0+" seconds");
			start = end;
			
            // transform matrix
            transform();
            end = System.currentTimeMillis();
			System.out.println("ST 2 Done in " + ((end-start)/100)/10.0+" seconds");
        }

        /**
         * Returns the matrix P of the transform.
         * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.</p>
         *
         * @return the P matrix
         */
        public RealMatrix getP() {
            if (cachedP == null) {
                //cachedP = MatrixUtils.createRealMatrix(matrixP);
                cachedP = new Array2DRowRealMatrix(n,  n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                    	cachedP.setEntry(i, j, matrixP[(i+n*j)]);
                    }
                }

            }
            return cachedP;
        }
        
        public double [] getPraw() {
        	return matrixP;
        }
        /**
         * Returns the transpose of the matrix P of the transform.
         * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.</p>
         *
         * @return the transpose of the P matrix
         */
        public RealMatrix getPT() {
            if (cachedPt == null) {
                cachedPt = getP().transpose();
            }

            // return the cached matrix
            return cachedPt;
        }

        /**
         * Returns the quasi-triangular Schur matrix T of the transform.
         *
         * @return the T matrix
         */
        public RealMatrix getT() {
            if (cachedT == null) {
                //cachedT = MatrixUtils.createRealMatrix(matrixT);
                cachedT = new Array2DRowRealMatrix(n,n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                    	cachedT.setEntry(i, j, matrixT[(i*n+j)]);
                    }
                }
            }

            // return the cached matrix
            return cachedT;
        }
        public double [] getTraw() {
        	return matrixT;
        }

        /**
         * Transform original matrix to Schur form.
         * @throws MaxCountExceededException if the transformation does not converge
         */
        private void transform() {

            // compute matrix norm
            final double norm = getNorm();

            // shift information
            final ShiftInfo shift = new ShiftInfo();

            // Outer loop over eigenvalue index
            int iteration = 0;
            int iu = n - 1;
            while (iu >= 0) {

                // Look for single small sub-diagonal element
                final int il = findSmallSubDiagonalElement(iu, norm);

                // Check for convergence
                if (il == iu) {
                    // One root found
                    matrixT[iu*n+iu] += shift.exShift;
                    iu--;
                    iteration = 0;
                } else if (il == iu - 1) {
                    // Two roots found
                    double p = (matrixT[(iu - 1)*n+(iu - 1)] - matrixT[iu*n+iu]) / 2.0;
                    double q = p * p + matrixT[iu*n+(iu - 1)] * matrixT[(iu - 1)*n+iu];
                    matrixT[iu*n+iu] += shift.exShift;
                    matrixT[(iu - 1)*n+(iu - 1)] += shift.exShift;

                    if (q >= 0) {
                        double z = FastMath.sqrt(FastMath.abs(q));
                        if (p >= 0) {
                            z = p + z;
                        } else {
                            z = p - z;
                        }
                        final double x = matrixT[iu*n+(iu - 1)];
                        final double s = FastMath.abs(x) + FastMath.abs(z);
                        p = x / s;
                        q = z / s;
                        final double r = FastMath.sqrt(p * p + q * q);
                        p /= r;
                        q /= r;

                        // Row modification
                        for (int j = iu - 1; j < n; j++) {
                            z = matrixT[(iu - 1)*n+j];
                            matrixT[(iu - 1)*n+j] = q * z + p * matrixT[iu*n+j];
                            matrixT[iu*n+j] = q * matrixT[iu*n+j] - p * z;
                        }

                        // Column modification
                        for (int i = 0; i <= iu; i++) {
                            z = matrixT[i*n+(iu - 1)];
                            matrixT[i*n+(iu - 1)] = q * z + p * matrixT[i*n+iu];
                            matrixT[i*n+iu] = q * matrixT[i*n+iu] - p * z;
                        }

                        // Accumulate transformations
                        for (int i = 0; i <= n - 1; i++) {
                            z = matrixP[(i+n*(iu - 1))];
                            matrixP[(i+n*(iu - 1))] = q * z + p * matrixP[(i+n*iu)];
                            matrixP[(i+n*iu)] = q * matrixP[(i+n*iu)] - p * z;
                        }
                    }
                    iu -= 2;
                    iteration = 0;
                } else {
                    // No convergence yet
                    computeShift(il, iu, iteration, shift);

                    // stop transformation after too many iterations
                    if (++iteration > MAX_ITERATIONS) {
                        throw new MaxCountExceededException(LocalizedFormats.CONVERGENCE_FAILED,
                                                            MAX_ITERATIONS);
                    }

                    // the initial houseHolder vector for the QR step
                    final double[] hVec = new double[(3)];

                    final int im = initQRStep(il, iu, shift, hVec);
                    performDoubleQRStep(il, im, iu, shift, hVec);
                }
            }
        }

        /**
         * Computes the L1 norm of the (quasi-)triangular matrix T.
         *
         * @return the L1 norm of matrix T
         */
        private double getNorm() {
            double norm = 0.0;
            for (int i = 0; i < n; i++) {
                // as matrix T is (quasi-)triangular, also take the sub-diagonal element into account
                for (int j = FastMath.max(i - 1, 0); j < n; j++) {
                    norm += FastMath.abs(matrixT[i*n+j]);
                }
            }
            return norm;
        }

        /**
         * Find the first small sub-diagonal element and returns its index.
         *
         * @param startIdx the starting index for the search
         * @param norm the L1 norm of the matrix
         * @return the index of the first small sub-diagonal element
         */
        private int findSmallSubDiagonalElement(final int startIdx, final double norm) {
            int l = startIdx;
            while (l > 0) {
                double s = FastMath.abs(matrixT[(l - 1)*n+(l - 1)]) + FastMath.abs(matrixT[(l)*n+(l)]);
                if (s == 0.0) {
                    s = norm;
                }
                if (FastMath.abs(matrixT[(l)*n+(l - 1)]) < epsilon * s) {
                    break;
                }
                l--;
            }
            return l;
        }

        /**
         * Compute the shift for the current iteration.
         *
         * @param l the index of the small sub-diagonal element
         * @param idx the current eigenvalue index
         * @param iteration the current iteration
         * @param shift holder for shift information
         */
        private void computeShift(final int l, final int idx, final int iteration, final ShiftInfo shift) {
            // Form shift
            shift.x = matrixT[(idx)*n+(idx)];
            shift.y = shift.w = 0.0;
            if (l < idx) {
                shift.y = matrixT[(idx - 1)*n+(idx - 1)];
                shift.w = matrixT[(idx)*n+(idx - 1)] * matrixT[(idx - 1)*n+(idx)];
            }

            // Wilkinson's original ad hoc shift
            if (iteration == 10) {
                shift.exShift += shift.x;
                for (int i = 0; i <= idx; i++) {
                    matrixT[i*n+(i)] -= shift.x;
                }
                final double s = FastMath.abs(matrixT[(idx)*n+(idx - 1)]) + FastMath.abs(matrixT[(idx - 1)*n+(idx - 2)]);
                shift.x = 0.75 * s;
                shift.y = 0.75 * s;
                shift.w = -0.4375 * s * s;
            }

            // MATLAB's new ad hoc shift
            if (iteration == 30) {
                double s = (shift.y - shift.x) / 2.0;
                s = s * s + shift.w;
                if (s > 0.0) {
                    s = FastMath.sqrt(s);
                    if (shift.y < shift.x) {
                        s = -s;
                    }
                    s = shift.x - shift.w / ((shift.y - shift.x) / 2.0 + s);
                    for (int i = 0; i <= idx; i++) {
                        matrixT[i*n+(i)] -= s;
                    }
                    shift.exShift += s;
                    shift.x = shift.y = shift.w = 0.964;
                }
            }
        }

        /**
         * Initialize the householder vectors for the QR step.
         *
         * @param il the index of the small sub-diagonal element
         * @param iu the current eigenvalue index
         * @param shift shift information holder
         * @param hVec the initial houseHolder vector
         * @return the start index for the QR step
         */
        private int initQRStep(int il, final int iu, final ShiftInfo shift, double[] hVec) {
            // Look for two consecutive small sub-diagonal elements
            int im = iu - 2;
            while (im >= il) {
                final double z = matrixT[(im)*n+(im)];
                final double r = shift.x - z;
                double s = shift.y - z;
                hVec[(0)] = (r * s - shift.w) / matrixT[(im + 1)*n+(im)] + matrixT[(im)*n+(im + 1)];
                hVec[(1)] = matrixT[(im + 1)*n+(im + 1)] - z - r - s;
                hVec[(2)] = matrixT[(im + 2)*n+(im + 1)];

                if (im == il) {
                    break;
                }

                final double lhs = FastMath.abs(matrixT[(im)*n+(im - 1)]) * (FastMath.abs(hVec[(1)]) + FastMath.abs(hVec[(2)]));
                final double rhs = FastMath.abs(hVec[(0)]) * (FastMath.abs(matrixT[(im - 1)*n+(im - 1)]) +
                                                            FastMath.abs(z) +
                                                            FastMath.abs(matrixT[(im + 1)*n+(im + 1)]));

                if (lhs < epsilon * rhs) {
                    break;
                }
                im--;
            }

            return im;
        }

        /**
         * Perform a double QR step involving rows l:idx and columns m:n
         *
         * @param il the index of the small sub-diagonal element
         * @param im the start index for the QR step
         * @param iu the current eigenvalue index
         * @param shift shift information holder
         * @param hVec the initial houseHolder vector
         */
        private void performDoubleQRStep(final int il, final int im, final int iu,
                                         final ShiftInfo shift, final double[] hVec) {

            double p = hVec[(0)];
            double q = hVec[(1)];
            double r = hVec[(2)];

            for (int k = im; k <= iu - 1; k++) {
                boolean notlast = k != (iu - 1);
                if (k != im) {
                    p = matrixT[k*n+(k - 1)];
                    q = matrixT[(k + 1)*n+(k - 1)];
                    r = notlast ? matrixT[(k + 2)*n+(k - 1)] : 0.0;
                    shift.x = FastMath.abs(p) + FastMath.abs(q) + FastMath.abs(r);
                    if (Precision.equals(shift.x, 0.0, epsilon)) {
                        continue;
                    }
                    p /= shift.x;
                    q /= shift.x;
                    r /= shift.x;
                }
                double s = FastMath.sqrt(p * p + q * q + r * r);
                if (p < 0.0) {
                    s = -s;
                }
                if (s != 0.0) {
                    if (k != im) {
                        matrixT[k*n+(k - 1)] = -s * shift.x;
                    } else if (il != im) {
                        matrixT[k*n+(k - 1)] = -matrixT[k*n+(k - 1)];
                    }
                    p += s;
                    shift.x = p / s;
                    shift.y = q / s;
                    double z = r / s;
                    q /= p;
                    r /= p;

                    // Row modification
                    for (int j = k; j < n; j++) {
                        p = matrixT[k*n+j] + q * matrixT[(k + 1)*n+j];
                        if (notlast) {
                            p += r * matrixT[(k + 2)*n+j];
                            matrixT[(k + 2)*n+j] -= p * z;
                        }
                        matrixT[k*n+j] -= p * shift.x;
                        matrixT[(k + 1)*n+j] -= p * shift.y;
                    }

                    // Column modification
                    for (int i = 0; i <= FastMath.min(iu, k + 3); i++) {
                        p = shift.x * matrixT[i*n+k] + shift.y * matrixT[i*n+(k + 1)];
                        if (notlast) {
                            p += z * matrixT[i*n+(k + 2)];
                            matrixT[i*n+(k + 2)] -= p * r;
                        }
                        matrixT[i*n+k] -= p;
                        matrixT[i*n+(k + 1)] -= p * q;
                    }

                    // Accumulate transformations
                    final int high = n - 1;
                    for (int i = 0; i <= high; i++) {
                        p = shift.x * matrixP[(i+n*k)] + shift.y * matrixP[(i+n*(k + 1))];
                        if (notlast) {
                            p += z * matrixP[(i+n*(k + 2))];
                            matrixP[(i+n*(k + 2))] -= p * r;
                        }
                        matrixP[(i+n*k)] -= p;
                        matrixP[(i+n*(k + 1))] -= p * q;
                    }
                }  // (s != 0)
            }  // k loop

            // clean up pollution due to round-off errors
            for (int i = im + 2; i <= iu; i++) {
                matrixT[i*n+(i-2)] = 0.0;
                if (i > im + 2) {
                    matrixT[i*n+(i-3)] = 0.0;
                }
            }
        }

    }

    /**
     * Internal data structure holding the current shift information.
     * Contains variable names as present in the original JAMA code.
     */
    private static class ShiftInfo {
        // CHECKSTYLE: stop all

        /** x shift info */
        double x;
        /** y shift info */
        double y;
        /** w shift info */
        double w;
        /** Indicates an exceptional shift. */
        double exShift;

        // CHECKSTYLE: resume all
    }
    
    
	class MyRealMatrix {
		double [] data;
		MyRealMatrix(double [] data) {
			this.data = data;
		}
		public double[] getDataR() {
			return data;
		}
		public double [][] getData() {
			double [][] matrix = new double[n][n];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					matrix[i][j] = data[(i * n + j)];
				}
			}
			return matrix;
		}
	}
    //householderVectors = new double[(n*n)];
   //	System.arraycopy(matrix, 0, householderVectors, 0, n * n);

	
	/**
	 * Class transforming a general real matrix to Hessenberg form.
	 * <p>A m &times; m matrix A can be written as the product of three matrices: A = P
	 * &times; H &times; P<sup>T</sup> with P an orthogonal matrix and H a Hessenberg
	 * matrix. Both P and H are m &times; m matrices.</p>
	 * <p>Transformation to Hessenberg form is often not a goal by itself, but it is an
	 * intermediate step in more general decomposition algorithms like
	 * {@link EigenDecomposition eigen decomposition}. This class is therefore
	 * intended for internal use by the library and is not public. As a consequence
	 * of this explicitly limited scope, many methods directly returns references to
	 * internal arrays, not copies.</p>
	 * <p>This class is based on the method orthes in class EigenvalueDecomposition
	 * from the <a href="http://math.nist.gov/javanumerics/jama/">JAMA</a> library.</p>
	 *
	 * @see <a href="http://mathworld.wolfram.com/HessenbergDecomposition.html">MathWorld</a>
	 * @see <a href="http://en.wikipedia.org/wiki/Householder_transformation">Householder Transformations</a>
	 * @version $Id: HessenbergTransformer.java 1538368 2013-11-03 13:57:37Z erans $
	 * @since 3.1
	 */
	class HessenbergTransformer {
	    /** Householder vectors. */
	    private final double householderVectors[];
	    /** Temporary storage vector. */
	    private final double ort[];
	    /** Cached value of P. */
	    private MyRealMatrix cachedP;
	    /** Cached value of Pt. */
	    private MyRealMatrix cachedPt;
	    /** Cached value of H. */
	    private MyRealMatrix cachedH;

	    /**
	     * Build the transformation to Hessenberg form of a general matrix.
	     *
	     * @param matrix matrix to transform
	     * @throws NonSquareMatrixException if the matrix is not square
	     */
	    public HessenbergTransformer(final double[] matrix) {
//	        if (!matrix.isSquare()) {
//	            throw new NonSquareMatrixException(matrix.getRowDimension(),
//	                    matrix.getColumnDimension());
//	        }

	        householderVectors = new double[(n*n)];
        	System.arraycopy(matrix, 0, householderVectors, 0, n * n);
	        
	        ort = new double[n];
	        cachedP = null;
	        cachedPt = null;
	        cachedH = null;

	        // transform matrix
	        transform();
	    }

	    /**
	     * Returns the matrix P of the transform.
	     * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.</p>
	     *
	     * @return the P matrix
	     */
	    public MyRealMatrix getP() {
	        if (cachedP == null) {
	            final int high = n - 1;
	            final double[] pa = new double[(n*n)];

	            for (int i = 0; i < n; i++) {
	                //for (int j = 0; j < n; j++) {
	                //    pa[(i+n*j)] = (i == j) ? 1 : 0;
	                //}
                    pa[(i+n*i)] = 1;
	            }

	            for (int m = high - 1; m >= 1; m--) {
	                if (householderVectors[(m*n+m - 1)] != 0.0) {
	                    for (int i = m + 1; i <= high; i++) {
	                        ort[i] = householderVectors[i*n + m - 1];
	                    }

	                    for (int j = m; j <= high; j++) {
	                        double g = 0.0;

	                        for (int i = m; i <= high; i++) {
	                            g += ort[i] * pa[(i+n* j)];
	                        }

	                        // Double division avoids possible underflow
	                        g = (g / ort[(m)]) / householderVectors[(m*n + m - 1)];

	                        for (int i = m; i <= high; i++) {
	                            pa[(i+n * j)] += g * ort[i];
	                        }
	                    }
	                }
	            }

	            cachedP = new MyRealMatrix(pa);//.createRealMatrix(pa);
	        }
	        return cachedP;
	    }

	    /**
	     * Returns the transpose of the matrix P of the transform.
	     * <p>P is an orthogonal matrix, i.e. its inverse is also its transpose.</p>
	     *
	     * @return the transpose of the P matrix
	     */
//	    public RealMatrix getPT() {
//	        if (cachedPt == null) {
//	            cachedPt = getP().transpose();
//	        }
//
//	        // return the cached matrix
//	        return cachedPt;
//	    }

	    /**
	     * Returns the Hessenberg matrix H of the transform.
	     *
	     * @return the H matrix
	     */
	    public MyRealMatrix getH() {
	        if (cachedH == null) {
	            final int m = n;
	            final double[] h = new double[n*n];
	            for (int i = 0; i < m; ++i) {
	                if (i > 0) {
	                    // copy the entry of the lower sub-diagonal
	                    h[i*n+i - 1] = householderVectors[i*n+i - 1];
	                }

	                // copy upper triangular part of the matrix
	                //for (int j = i; j < m; ++j) {
	                //    h[i*n+j] = householderVectors[i*n+j];
	                //}
	                System.arraycopy(householderVectors, i*n+i, h, i*n+i, m - i);
	            }
	            cachedH = new MyRealMatrix(h);//.createRealMatrix(h);
	        }

	        // return the cached matrix
	        return cachedH;
	    }

	    /**
	     * Get the Householder vectors of the transform.
	     * <p>Note that since this class is only intended for internal use, it returns
	     * directly a reference to its internal arrays, not a copy.</p>
	     *
	     * @return the main diagonal elements of the B matrix
	     */
	    double[] getHouseholderVectorsRef() {
	        return householderVectors;
	    }

	    /**
	     * Transform original matrix to Hessenberg form.
	     * <p>Transformation is done using Householder transforms.</p>
	     */
	    private void transform() {
	        final int high = n - 1;

	        for (int m = 1; m <= high - 1; m++) {
	            // Scale column.
	            double scale = 0;
	            for (int i = m; i <= high; i++) {
	                scale += FastMath.abs(householderVectors[i*n+m - 1]);
	            }

	            if (!Precision.equals(scale, 0)) {
	                // Compute Householder transformation.
	                double h = 0;
	                for (int i = high; i >= m; i--) {
	                    ort[i] = householderVectors[i*n+m - 1] / scale;
	                    h += ort[i] * ort[i];
	                }
	                final double g = (ort[(m)] > 0) ? -FastMath.sqrt(h) : FastMath.sqrt(h);

	                h -= ort[(m)] * g;
	                ort[(m)] -= g;

	                // Apply Householder similarity transformation
	                // H = (I - u*u' / h) * H * (I - u*u' / h)

	                for (int j = m; j < n; j++) {
	                    double f = 0;
	                    for (int i = high; i >= m; i--) {
	                        f += ort[i] * householderVectors[i*n+j];
	                    }
	                    f /= h;
	                    for (int i = m; i <= high; i++) {
	                        householderVectors[i*n+j] -= f * ort[i];
	                    }
	                }

	                for (int i = 0; i <= high; i++) {
	                    double f = 0;
	                    for (int j = high; j >= m; j--) {
	                        f += ort[j] * householderVectors[i*n+j];
	                    }
	                    f /= h;
	                    for (int j = m; j <= high; j++) {
	                        householderVectors[i*n+j] -= f * ort[j];
	                    }
	                }

	                ort[(m)] = scale * ort[(m)];
	                householderVectors[m*n+m - 1] = scale * g;
	            }
	        }
	    }
	}

    
}
