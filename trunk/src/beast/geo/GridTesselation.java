package beast.geo;


import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

import beast.continuous.SphericalDiffusionModel;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.substitutionmodel.DefaultEigenSystem;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.EigenSystem;
import beast.util.Randomizer;

@Description("Tesselates (part of) a sphere with equal sized quadrangle on a grid")
public class GridTesselation extends SphereTesselation {

		public Input<File> mapFileInput = new Input<File>("map","geographic world map in Mercator projection -- red coloured items will have distance defined by 'reddistance'. "
				+ "If not specified, all distances are 1.0",
				new File("/home/remco/data/geo/aboriginal25.bmp"));
		public Input<Double> reddistanceInput = new Input<Double>("reddistance","distance fpr red coloured items", 1.0);
		
		GraphNode [][] grid;
		
		/** flag indicating we use native code **/
		boolean m_bNative = false;
		
		public GridTesselation() {
			bboxInput.setRule(Validate.REQUIRED);
		}
	
		@Override
		public void initAndValidate() throws Exception {
			
			parseBBox();
			Quadrangle q = new Quadrangle(minLat, minLong, maxLat, maxLong);
			double [] center = q.getCenter();
			
			int depth = depthInput.get();
			if (depth <= 0) {
				throw new RuntimeException("depth must be a positive number");
			}
			
			// create vertices
			double long0 = minLong - center[1];
			double longStep = (maxLong - minLong)/depth;
			double lat0 = minLat - center[0];
			double latStep = (maxLat - minLat)/depth;
			Vertex [][] vertex = new Vertex[depth][depth];
			for (int i = 0; i < depth; i++) {
				double lat1 = lat0 + latStep * i;
				for (int j = 0; j < depth; j++) {
					double long1 = long0 + longStep * j;
					double [] point = SphericalDiffusionModel.reverseMap(lat1, long1, center[0], center[1]);
					vertex[i][j] = new Vertex(point[0], point[1]);
				}
			}
			
			grid = new GraphNode[depth-1][depth-1];
			nodes = new ArrayList<GraphNode>();
			for (int i = 0; i < depth - 1; i++) {
				for (int j = 0; j < depth - 1; j++) {
					Quadrangle q1 = new Quadrangle(vertex[i][j],vertex[i+1][j], vertex[i+1][j+1],vertex[i][j+1]);
					nodes.add(q1);
					grid[i][j] = q1;
				}
			}
			
			
			boolean all = allNeighborsInput.get();
			for (int i = 0; i < depth - 1; i++) {
				for (int j = 0; j < depth - 1; j++) {
					List<GraphNode> neighbors = new ArrayList<GraphNode>();
					if (i > 0) {
						neighbors.add(grid[i-1][j]);
					}
					if (all && i > 0 && j < depth-2) {
						neighbors.add(grid[i-1][j + 1]);
					}
					if (j < depth-2) {
						neighbors.add(grid[i][j+1]);
					}
					if (all && i < depth - 2 && j < depth-2) {
						neighbors.add(grid[i+1][j + 1]);
					}
					if (i < depth-2) {
						neighbors.add(grid[i+1][j]);
					}
					if (all && i < depth - 2 && j > 0) {
						neighbors.add(grid[i+1][j - 1]);
					}
					if (j > 0) {
						neighbors.add(grid[i][j-1]);
					}
					if (all && i > 0  && j > 0) {
						neighbors.add(grid[i-1][j - 1]);
					}
					grid[i][j].neighbours = neighbors.toArray(new GraphNode[]{});
					grid[i][j].setUpDistances(useGreatCircleInput.get());
				}
			}
			
			System.err.println("#nodes = " + nodes.size());
	
			// renumber remaining quadrangles
			renumber();
	
			// adjust distances to map
			final BufferedImage image = ImageIO.read(mapFileInput.get());
			int w = image.getWidth();
			int h = image.getHeight();

			double redDistance = reddistanceInput.get();
			for (GraphNode t : nodes) {
				double [] c = t.getCenter();
				int x =(int)( w * (c[1]+180) / 360.0);
				int y =(int)( h * (c[0]+90) / 180.0);
				int color = image.getRGB(x, y) & 0xFFFFFF;
				if (color == 0x00FF00) {
					t.scaleDistance(redDistance);
				}
			}

			setUpLatLongMap();
			
			// log some stats
			System.err.println("#nodes= " + nodes.size());
			
			
			
			
			int n = nodes.size();
			RealMatrix m = new Array2DRowRealMatrix(n, n);
			//RealMatrix m = new OpenMapRealMatrix(n,n);
			int [][] matrix = new int[n][n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					matrix[t.id][node.id] = 1; 
					m.setEntry(t.id, node.id, 1);
				}
				matrix[t.id][t.id] = -t.neighbours.length; 
				m.setEntry(t.id, t.id, -t.neighbours.length);
			}
//			int [][] matrix2 = new int[n][n];
//			for (int i = 0; i < n; i++) {
//				for (int j = 0; j < n; j++) {
//					int sum = 0;
//					for (int k = 0; k < n; k++) {
//						sum += matrix[i][k] * matrix[k][j];
//					}
//					matrix2[i][j] = sum;
//				}
//			}
//			System.out.println();
//			for (int i = 0; i < n; i++) {
//				System.out.println(Arrays.toString(matrix[i]));
//			}
//			System.out.println();
//			for (int i = 0; i < n; i++) {
//				System.out.println(Arrays.toString(matrix2[i]));
//			}
			long start = System.currentTimeMillis();
			//singularValueDecomposition5();
			//SingularValueDecomposition svd = new SingularValueDecomposition(m);
//			RealMatrix U = svd.getU();
//			RealMatrix S = svd.getS();
//			RealMatrix V = svd.getV();
			long end = System.currentTimeMillis();
			
//			System.out.println(Arrays.toString(U));
//			System.out.println(Arrays.toString(singularValues));
//			System.out.println(Arrays.toString(V));
			//System.out.println("Done in " + (end-start)/1000+" seconds");

			start = System.currentTimeMillis();
			singularValueDecomposition6();
			end = System.currentTimeMillis();
//
//			System.out.println(Arrays.toString(U));
			System.out.println(Arrays.toString(S));
//			System.out.println(Arrays.toString(V));
			System.out.println("Done in " + ((end-start)/100)/10.0+" seconds");
			
			
			// sanity check
			RealMatrix Um = new Array2DRowRealMatrix(n, n); 
			RealMatrix Sm = new Array2DRowRealMatrix(n, n); 
			RealMatrix Vm = new Array2DRowRealMatrix(n, n);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					Um.setEntry(i,  j, U[j*n+i]);
					Vm.setEntry(i,  j, V[i*n+j]);
				}
				Sm.setEntry(i, i, S[i]);
			}
			
			RealMatrix P = Vm.multiply(Sm).multiply(Um.transpose());
			double max = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					max = Math.max(max, Math.abs(matrix[i][j] - P.getEntry(i, j)));
				}
			}
			System.out.println("Max abs differnce: " + max);
			
			//P = Vm.multiply(Um.transpose());
			//System.out.println(P);
			
			
			System.exit(0);
			
			
			
		}
		
	
		
	    /** Relative threshold for small singular values. */
	    //private static final double EPS = 0x1.0p-52;
	    private static final double EPS = 0x1.0p-26;
	    /** Absolute threshold for small singular values. */
	    //private static final double TINY = 0x1.0p-966;
	    private static final double TINY = 0x1.0p-52;

	    /** Relative threshold for small singular values. */
	    private static final float EPSF = 0x1.0p-26f;
	    /** Absolute threshold for small singular values. */
	    private static final float TINYF = 0x1.0p-52f;

	    
	    /** Computed singular values. */
	    private double[] S;
	    double[] U, V;
	    private float[] singularValuesf;
	    float[] Uf, Vf;
	    /** row dimension = column dimension */
	    private int n;
	    /** Indicator for transposed matrix. */
	    //private boolean transposed;
	    /**
	     * Tolerance value for small singular values, calculated once we have
	     * populated "singularValues".
	     **/
	    private double tol;


		/* adapted from org.apache.commons.math3.linear.SingularValueDecomposition */
	    public void singularValueDecomposition() {
	    	System.err.println("java apache matrix=array of array SVD");
			long start = System.currentTimeMillis();
			int n = nodes.size();
	        final double[][] A = new double[n][n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id][node.id] = 1; 
				}
				A[t.id][t.id] = -t.neighbours.length; 
			}


	         // "m" is always the largest dimension.
            //transposed = false;

	        S = new double[n];
	        final double[][] U = new double[n][n];
	        final double[][] V = new double[n][n];
	        final double[] e = new double[n];
	        final double[] work = new double[n];
	        // Reduce A to bidiagonal form, storing the diagonal elements
	        // in s and the super-diagonal elements in e.
	        final int nct = FastMath.min(n - 1, n);
	        final int nrt = FastMath.max(0, n - 2);
	        for (int k = 0; k < FastMath.max(nct, nrt); k++) {
	            if (k < nct) {
	                // Compute the transformation for the k-th column and
	                // place the k-th diagonal in s[k].
	                // Compute 2-norm of k-th column without under/overflow.
	                S[k] = 0;
	                for (int i = k; i < n; i++) {
	                    S[k] = FastMath.hypot(S[k], A[i][k]);
	                }
	                if (S[k] != 0) {
	                    if (A[k][k] < 0) {
	                        S[k] = -S[k];
	                    }
	                    for (int i = k; i < n; i++) {
	                        A[i][k] /= S[k];
	                    }
	                    A[k][k] += 1;
	                }
	                S[k] = -S[k];
	            }
	            for (int j = k + 1; j < n; j++) {
	                if (k < nct &&
	                    S[k] != 0) {
	                    // Apply the transformation.
	                    double t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += A[i][k] * A[i][j];
	                    }
	                    t = -t / A[k][k];
	                    for (int i = k; i < n; i++) {
	                        A[i][j] += t * A[i][k];
	                    }
	                }
	                // Place the k-th row of A into e for the
	                // subsequent calculation of the row transformation.
	                e[j] = A[k][j];
	            }
	            if (k < nct) {
	                // Place the transformation in U for subsequent back
	                // multiplication.
	                for (int i = k; i < n; i++) {
	                    U[i][k] = A[i][k];
	                }
	            }
	            if (k < nrt) {
	                // Compute the k-th row transformation and place the
	                // k-th super-diagonal in e[k].
	                // Compute 2-norm without under/overflow.
	                e[k] = 0;
	                for (int i = k + 1; i < n; i++) {
	                    e[k] = FastMath.hypot(e[k], e[i]);
	                }
	                if (e[k] != 0) {
	                    if (e[k + 1] < 0) {
	                        e[k] = -e[k];
	                    }
	                    for (int i = k + 1; i < n; i++) {
	                        e[i] /= e[k];
	                    }
	                    e[k + 1] += 1;
	                }
	                e[k] = -e[k];
	                if (k + 1 < n &&
	                    e[k] != 0) {
	                    // Apply the transformation.
	                    for (int i = k + 1; i < n; i++) {
	                        work[i] = 0;
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        for (int i = k + 1; i < n; i++) {
	                            work[i] += e[j] * A[i][j];
	                        }
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        final double t = -e[j] / e[k + 1];
	                        for (int i = k + 1; i < n; i++) {
	                            A[i][j] += t * work[i];
	                        }
	                    }
	                }

	                // Place the transformation in V for subsequent
	                // back multiplication.
	                for (int i = k + 1; i < n; i++) {
	                    V[i][k] = e[i];
	                }
	            }
	        }
	        
	        
			long end = System.currentTimeMillis();			
			System.out.println("Step 1 in " + (end-start)/1000+" seconds");
			start = end;
			
	        
	        // Set up the final bidiagonal matrix or order p.
	        int p = n;
	        if (nct < n) {
	            S[nct] = A[nct][nct];
	        }
	        if (n < p) {
	            S[p - 1] = 0;
	        }
	        if (nrt + 1 < p) {
	            e[nrt] = A[nrt][p - 1];
	        }
	        e[p - 1] = 0;

	        // Generate U.
	        for (int j = nct; j < n; j++) {
	            for (int i = 0; i < n; i++) {
	                U[i][j] = 0;
	            }
	            U[j][j] = 1;
	        }
	        for (int k = nct - 1; k >= 0; k--) {
	            if (S[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    double t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += U[i][k] * U[i][j];
	                    }
	                    t = -t / U[k][k];
	                    for (int i = k; i < n; i++) {
	                        U[i][j] += t * U[i][k];
	                    }
	                }
	                for (int i = k; i < n; i++) {
	                    U[i][k] = -U[i][k];
	                }
	                U[k][k] = 1 + U[k][k];
	                for (int i = 0; i < k - 1; i++) {
	                    U[i][k] = 0;
	                }
	            } else {
	                for (int i = 0; i < n; i++) {
	                    U[i][k] = 0;
	                }
	                U[k][k] = 1;
	            }
	        }

	        // Generate V.
	        for (int k = n - 1; k >= 0; k--) {
	            if (k < nrt &&
	                e[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    double t = 0;
	                    for (int i = k + 1; i < n; i++) {
	                        t += V[i][k] * V[i][j];
	                    }
	                    t = -t / V[k + 1][k];
	                    for (int i = k + 1; i < n; i++) {
	                        V[i][j] += t * V[i][k];
	                    }
	                }
	            }
	            for (int i = 0; i < n; i++) {
	                V[i][k] = 0;
	            }
	            V[k][k] = 1;
	        }

	        // Main iteration loop for the singular values.
	        final int pp = p - 1;
	        int iter = 0;
	        while (p > 0) {
	            int k;
	            int kase;
	            // Here is where a test for too many iterations would go.
	            // This section of the program inspects for
	            // negligible elements in the s and e arrays.  On
	            // completion the variables kase and k are set as follows.
	            // kase = 1     if s(p) and e[k-1] are negligible and k<p
	            // kase = 2     if s(k) is negligible and k<p
	            // kase = 3     if e[k-1] is negligible, k<p, and
	            //              s(k), ..., s(p) are not negligible (qr step).
	            // kase = 4     if e(p-1) is negligible (convergence).
	            for (k = p - 2; k >= 0; k--) {
	                final double threshold
	                    = TINY + EPS * (FastMath.abs(S[k]) +
	                                    FastMath.abs(S[k + 1]));

	                // the following condition is written this way in order
	                // to break out of the loop when NaN occurs, writing it
	                // as "if (FastMath.abs(e[k]) <= threshold)" would loop
	                // indefinitely in case of NaNs because comparison on NaNs
	                // always return false, regardless of what is checked
	                // see issue MATH-947
	                if (!(FastMath.abs(e[k]) > threshold)) {
	                    e[k] = 0;
	                    break;
	                }

	            }

	            if (k == p - 2) {
	                kase = 4;
	            } else {
	                int ks;
	                for (ks = p - 1; ks >= k; ks--) {
	                    if (ks == k) {
	                        break;
	                    }
	                    final double t = (ks != p ? FastMath.abs(e[ks]) : 0) +
	                        (ks != k + 1 ? FastMath.abs(e[ks - 1]) : 0);
	                    if (FastMath.abs(S[ks]) <= TINY + EPS * t) {
	                        S[ks] = 0;
	                        break;
	                    }
	                }
	                if (ks == k) {
	                    kase = 3;
	                } else if (ks == p - 1) {
	                    kase = 1;
	                } else {
	                    kase = 2;
	                    k = ks;
	                }
	            }
	            k++;
	            // Perform the task indicated by kase.
	            switch (kase) {
	                // Deflate negligible s(p).
	                case 1: {
	                    double f = e[p - 2];
	                    e[p - 2] = 0;
	                    for (int j = p - 2; j >= k; j--) {
	                        double t = FastMath.hypot(S[j], f);
	                        final double cs = S[j] / t;
	                        final double sn = f / t;
	                        S[j] = t;
	                        if (j != k) {
	                            f = -sn * e[j - 1];
	                            e[j - 1] = cs * e[j - 1];
	                        }

	                        for (int i = 0; i < n; i++) {
	                            t = cs * V[i][j] + sn * V[i][p - 1];
	                            V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
	                            V[i][j] = t;
	                        }
	                    }
	                }
	                break;
	                // Split at negligible s(k).
	                case 2: {
	                    double f = e[k - 1];
	                    e[k - 1] = 0;
	                    for (int j = k; j < p; j++) {
	                        double t = FastMath.hypot(S[j], f);
	                        final double cs = S[j] / t;
	                        final double sn = f / t;
	                        S[j] = t;
	                        f = -sn * e[j];
	                        e[j] = cs * e[j];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * U[i][j] + sn * U[i][k - 1];
	                            U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
	                            U[i][j] = t;
	                        }
	                    }
	                }
	                break;
	                // Perform one qr step.
	                case 3: {
	                    // Calculate the shift.
	                    final double maxPm1Pm2 = FastMath.max(FastMath.abs(S[p - 1]),
	                                                          FastMath.abs(S[p - 2]));
	                    final double scale = FastMath.max(FastMath.max(FastMath.max(maxPm1Pm2,
	                                                                                FastMath.abs(e[p - 2])),
	                                                                   FastMath.abs(S[k])),
	                                                      FastMath.abs(e[k]));
	                    final double sp = S[p - 1] / scale;
	                    final double spm1 = S[p - 2] / scale;
	                    final double epm1 = e[p - 2] / scale;
	                    final double sk = S[k] / scale;
	                    final double ek = e[k] / scale;
	                    final double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
	                    final double c = (sp * epm1) * (sp * epm1);
	                    double shift = 0;
	                    if (b != 0 ||
	                        c != 0) {
	                        shift = FastMath.sqrt(b * b + c);
	                        if (b < 0) {
	                            shift = -shift;
	                        }
	                        shift = c / (b + shift);
	                    }
	                    double f = (sk + sp) * (sk - sp) + shift;
	                    double g = sk * ek;
	                    // Chase zeros.
	                    for (int j = k; j < p - 1; j++) {
	                        double t = FastMath.hypot(f, g);
	                        double cs = f / t;
	                        double sn = g / t;
	                        if (j != k) {
	                            e[j - 1] = t;
	                        }
	                        f = cs * S[j] + sn * e[j];
	                        e[j] = cs * e[j] - sn * S[j];
	                        g = sn * S[j + 1];
	                        S[j + 1] = cs * S[j + 1];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * V[i][j] + sn * V[i][j + 1];
	                            V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
	                            V[i][j] = t;
	                        }
	                        t = FastMath.hypot(f, g);
	                        cs = f / t;
	                        sn = g / t;
	                        S[j] = t;
	                        f = cs * e[j] + sn * S[j + 1];
	                        S[j + 1] = -sn * e[j] + cs * S[j + 1];
	                        g = sn * e[j + 1];
	                        e[j + 1] = cs * e[j + 1];
	                        if (j < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = cs * U[i][j] + sn * U[i][j + 1];
	                                U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
	                                U[i][j] = t;
	                            }
	                        }
	                    }
	                    e[p - 2] = f;
	                    iter++;
	                }
	                break;
	                // Convergence.
	                default: {
	                    // Make the singular values positive.
	                    if (S[k] <= 0) {
	                        S[k] = S[k] < 0 ? -S[k] : 0;

	                        for (int i = 0; i <= pp; i++) {
	                            V[i][k] = -V[i][k];
	                        }
	                    }
	                    // Order the singular values.
	                    while (k < pp) {
	                        if (S[k] >= S[k + 1]) {
	                            break;
	                        }
	                        double t = S[k];
	                        S[k] = S[k + 1];
	                        S[k + 1] = t;
	                        if (k < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = V[i][k + 1];
	                                V[i][k + 1] = V[i][k];
	                                V[i][k] = t;
	                            }
	                        }
	                        if (k < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = U[i][k + 1];
	                                U[i][k + 1] = U[i][k];
	                                U[i][k] = t;
	                            }
	                        }
	                        k++;
	                    }
	                    iter = 0;
	                    p--;
	                }
	                break;
	            }
	        }

	        // Set the small value tolerance used to calculate rank and pseudo-inverse
	        tol = FastMath.max(n * S[0] * EPS,
	                           FastMath.sqrt(Precision.SAFE_MIN));
	        
	        
	        
			end = System.currentTimeMillis();			
			System.out.println("Step 2 in " + (end-start)/1000+" seconds");

	        
	        double max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(Math.abs(U[i][j])- Math.abs(V[i][j])));
		       }
	        }
	        System.out.println("max diff = " + max);
	        max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(U[i][j] + V[i][j]));
		       }
	        }
	        System.out.println("max diff = " + max);
//        	System.out.println(Arrays.toString(singularValues));
//	        System.out.println();
//	        for (int i = 0; i < n; i++) {
//	        	System.out.println(Arrays.toString(V[i]));
//	        }
	        System.out.println("tolerance = " + tol);


	    }
		
	    public void singularValueDecomposition2() {
	    	System.err.println("java apache matrix = array SVD");
			long start = System.currentTimeMillis();
			int n = nodes.size();
	        final double[] A = new double[n*n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id*n+node.id] = 1; 
				}
				A[t.id*n+t.id] = -t.neighbours.length; 
			}


	         // "m" is always the largest dimension.
            //transposed = false;

	        S = new double[n];
	        final double[] U = new double[n*n];
	        final double[] V = new double[n*n];
	        final double[] e = new double[n];
	        final double[] work = new double[n];
	        // Reduce A to bidiagonal form, storing the diagonal elements
	        // in s and the super-diagonal elements in e.
	        final int nct = FastMath.min(n - 1, n);
	        final int nrt = FastMath.max(0, n - 2);
	        for (int k = 0; k < FastMath.max(nct, nrt); k++) {
	            if (k < nct) {
	                // Compute the transformation for the k-th column and
	                // place the k-th diagonal in s[k].
	                // Compute 2-norm of k-th column without under/overflow.
	                S[k] = 0;
	                for (int i = k; i < n; i++) {
	                    S[k] = FastMath.hypot(S[k], A[i*n+k]);
	                }
	                if (S[k] != 0) {
	                    if (A[k*n+k] < 0) {
	                        S[k] = -S[k];
	                    }
	                    for (int i = k; i < n; i++) {
	                        A[i*n+k] /= S[k];
	                    }
	                    A[k*n+k] += 1;
	                }
	                S[k] = -S[k];
	            }
	            for (int j = k + 1; j < n; j++) {
	                if (k < nct &&
	                    S[k] != 0) {
	                    // Apply the transformation.
	                    double t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += A[i*n+k] * A[i*n+j];
	                    }
	                    t = -t / A[k*n+k];
	                    for (int i = k; i < n; i++) {
	                        A[i*n+j] += t * A[i*n+k];
	                    }
	                }
	                // Place the k-th row of A into e for the
	                // subsequent calculation of the row transformation.
	                e[j] = A[k*n+j];
	            }
	            if (k < nct) {
	                // Place the transformation in U for subsequent back
	                // multiplication.
	                for (int i = k; i < n; i++) {
	                    U[i*n+k] = A[i*n+k];
	                }
	            }
	            if (k < nrt) {
	                // Compute the k-th row transformation and place the
	                // k-th super-diagonal in e[k].
	                // Compute 2-norm without under/overflow.
	                e[k] = 0;
	                for (int i = k + 1; i < n; i++) {
	                    e[k] = FastMath.hypot(e[k], e[i]);
	                }
	                if (e[k] != 0) {
	                    if (e[k + 1] < 0) {
	                        e[k] = -e[k];
	                    }
	                    for (int i = k + 1; i < n; i++) {
	                        e[i] /= e[k];
	                    }
	                    e[k + 1] += 1;
	                }
	                e[k] = -e[k];
	                if (k + 1 < n &&
	                    e[k] != 0) {
	                    // Apply the transformation.
	                    for (int i = k + 1; i < n; i++) {
	                        work[i] = 0;
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        for (int i = k + 1; i < n; i++) {
	                            work[i] += e[j] * A[i*n+j];
	                        }
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        final double t = -e[j] / e[k + 1];
	                        for (int i = k + 1; i < n; i++) {
	                            A[i*n+j] += t * work[i];
	                        }
	                    }
	                }

	                // Place the transformation in V for subsequent
	                // back multiplication.
	                for (int i = k + 1; i < n; i++) {
	                    V[i*n+k] = e[i];
	                }
	            }
	        }
	        
	        
			long end = System.currentTimeMillis();			
			System.out.println("Step 1 in " + (end-start)/1000+" seconds");
			start = end;
			
	        
	        // Set up the final bidiagonal matrix or order p.
	        int p = n;
	        if (nct < n) {
	            S[nct] = A[nct*n+nct];
	        }
	        if (n < p) {
	            S[p - 1] = 0;
	        }
	        if (nrt + 1 < p) {
	            e[nrt] = A[nrt*n+p - 1];
	        }
	        e[p - 1] = 0;

	        // Generate U.
	        for (int j = nct; j < n; j++) {
	            for (int i = 0; i < n; i++) {
	                U[i*n+j] = 0;
	            }
	            U[j*n+j] = 1;
	        }
	        for (int k = nct - 1; k >= 0; k--) {
	            if (S[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    double t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += U[i*n+k] * U[i*n+j];
	                    }
	                    t = -t / U[k*n+k];
	                    for (int i = k; i < n; i++) {
	                        U[i*n+j] += t * U[i*n+k];
	                    }
	                }
	                for (int i = k; i < n; i++) {
	                    U[i*n+k] = -U[i*n+k];
	                }
	                U[k*n+k] = 1 + U[k*n+k];
	                for (int i = 0; i < k - 1; i++) {
	                    U[i*n+k] = 0;
	                }
	            } else {
	                for (int i = 0; i < n; i++) {
	                    U[i*n+k] = 0;
	                }
	                U[k*n+k] = 1;
	            }
	        }

	        // Generate V.
	        for (int k = n - 1; k >= 0; k--) {
	            if (k < nrt &&
	                e[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    double t = 0;
	                    for (int i = k + 1; i < n; i++) {
	                        t += V[i*n+k] * V[i*n+j];
	                    }
	                    t = -t / V[(k + 1)*n+k];
	                    for (int i = k + 1; i < n; i++) {
	                        V[i*n+j] += t * V[i*n+k];
	                    }
	                }
	            }
	            for (int i = 0; i < n; i++) {
	                V[i*n+k] = 0;
	            }
	            V[k*n+k] = 1;
	        }

	        // Main iteration loop for the singular values.
	        final int pp = p - 1;
	        int iter = 0;
	        while (p > 0) {
	            int k;
	            int kase;
	            // Here is where a test for too many iterations would go.
	            // This section of the program inspects for
	            // negligible elements in the s and e arrays.  On
	            // completion the variables kase and k are set as follows.
	            // kase = 1     if s(p) and e[k-1] are negligible and k<p
	            // kase = 2     if s(k) is negligible and k<p
	            // kase = 3     if e[k-1] is negligible, k<p, and
	            //              s(k), ..., s(p) are not negligible (qr step).
	            // kase = 4     if e(p-1) is negligible (convergence).
	            for (k = p - 2; k >= 0; k--) {
	                final double threshold
	                    = TINY + EPS * (FastMath.abs(S[k]) +
	                                    FastMath.abs(S[k + 1]));

	                // the following condition is written this way in order
	                // to break out of the loop when NaN occurs, writing it
	                // as "if (FastMath.abs(e[k]) <= threshold)" would loop
	                // indefinitely in case of NaNs because comparison on NaNs
	                // always return false, regardless of what is checked
	                // see issue MATH-947
	                if (!(FastMath.abs(e[k]) > threshold)) {
	                    e[k] = 0;
	                    break;
	                }

	            }

	            if (k == p - 2) {
	                kase = 4;
	            } else {
	                int ks;
	                for (ks = p - 1; ks >= k; ks--) {
	                    if (ks == k) {
	                        break;
	                    }
	                    final double t = (ks != p ? FastMath.abs(e[ks]) : 0) +
	                        (ks != k + 1 ? FastMath.abs(e[ks - 1]) : 0);
	                    if (FastMath.abs(S[ks]) <= TINY + EPS * t) {
	                        S[ks] = 0;
	                        break;
	                    }
	                }
	                if (ks == k) {
	                    kase = 3;
	                } else if (ks == p - 1) {
	                    kase = 1;
	                } else {
	                    kase = 2;
	                    k = ks;
	                }
	            }
	            k++;
	            // Perform the task indicated by kase.
	            switch (kase) {
	                // Deflate negligible s(p).
	                case 1: {
	                    double f = e[p - 2];
	                    e[p - 2] = 0;
	                    for (int j = p - 2; j >= k; j--) {
	                        double t = FastMath.hypot(S[j], f);
	                        final double cs = S[j] / t;
	                        final double sn = f / t;
	                        S[j] = t;
	                        if (j != k) {
	                            f = -sn * e[j - 1];
	                            e[j - 1] = cs * e[j - 1];
	                        }

	                        for (int i = 0; i < n; i++) {
	                            t = cs * V[i*n+j] + sn * V[i*n+p - 1];
	                            V[i*n+p - 1] = -sn * V[i*n+j] + cs * V[i*n+p - 1];
	                            V[i*n+j] = t;
	                        }
	                    }
	                }
	                break;
	                // Split at negligible s(k).
	                case 2: {
	                    double f = e[k - 1];
	                    e[k - 1] = 0;
	                    for (int j = k; j < p; j++) {
	                        double t = FastMath.hypot(S[j], f);
	                        final double cs = S[j] / t;
	                        final double sn = f / t;
	                        S[j] = t;
	                        f = -sn * e[j];
	                        e[j] = cs * e[j];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * U[i*n+j] + sn * U[i*n+k - 1];
	                            U[i*n+k - 1] = -sn * U[i*n+j] + cs * U[i*n+k - 1];
	                            U[i*n+j] = t;
	                        }
	                    }
	                }
	                break;
	                // Perform one qr step.
	                case 3: {
	                    // Calculate the shift.
	                    final double maxPm1Pm2 = FastMath.max(FastMath.abs(S[p - 1]),
	                                                          FastMath.abs(S[p - 2]));
	                    final double scale = FastMath.max(FastMath.max(FastMath.max(maxPm1Pm2,
	                                                                                FastMath.abs(e[p - 2])),
	                                                                   FastMath.abs(S[k])),
	                                                      FastMath.abs(e[k]));
	                    final double sp = S[p - 1] / scale;
	                    final double spm1 = S[p - 2] / scale;
	                    final double epm1 = e[p - 2] / scale;
	                    final double sk = S[k] / scale;
	                    final double ek = e[k] / scale;
	                    final double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
	                    final double c = (sp * epm1) * (sp * epm1);
	                    double shift = 0;
	                    if (b != 0 ||
	                        c != 0) {
	                        shift = FastMath.sqrt(b * b + c);
	                        if (b < 0) {
	                            shift = -shift;
	                        }
	                        shift = c / (b + shift);
	                    }
	                    double f = (sk + sp) * (sk - sp) + shift;
	                    double g = sk * ek;
	                    // Chase zeros.
	                    for (int j = k; j < p - 1; j++) {
	                        double t = FastMath.hypot(f, g);
	                        double cs = f / t;
	                        double sn = g / t;
	                        if (j != k) {
	                            e[j - 1] = t;
	                        }
	                        f = cs * S[j] + sn * e[j];
	                        e[j] = cs * e[j] - sn * S[j];
	                        g = sn * S[j + 1];
	                        S[j + 1] = cs * S[j + 1];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * V[i*n+j] + sn * V[i*n+j + 1];
	                            V[i*n+j + 1] = -sn * V[i*n+j] + cs * V[i*n+j + 1];
	                            V[i*n+j] = t;
	                        }
	                        t = FastMath.hypot(f, g);
	                        cs = f / t;
	                        sn = g / t;
	                        S[j] = t;
	                        f = cs * e[j] + sn * S[j + 1];
	                        S[j + 1] = -sn * e[j] + cs * S[j + 1];
	                        g = sn * e[j + 1];
	                        e[j + 1] = cs * e[j + 1];
	                        if (j < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = cs * U[i*n+j] + sn * U[i*n+j + 1];
	                                U[i*n+j + 1] = -sn * U[i*n+j] + cs * U[i*n+j + 1];
	                                U[i*n+j] = t;
	                            }
	                        }
	                    }
	                    e[p - 2] = f;
	                    iter++;
	                }
	                break;
	                // Convergence.
	                default: {
	                    // Make the singular values positive.
	                    if (S[k] <= 0) {
	                        S[k] = S[k] < 0 ? -S[k] : 0;

	                        for (int i = 0; i <= pp; i++) {
	                            V[i*n+k] = -V[i*n+k];
	                        }
	                    }
	                    // Order the singular values.
	                    while (k < pp) {
	                        if (S[k] >= S[k + 1]) {
	                            break;
	                        }
	                        double t = S[k];
	                        S[k] = S[k + 1];
	                        S[k + 1] = t;
	                        if (k < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = V[i*n+k + 1];
	                                V[i*n+k + 1] = V[i*n+k];
	                                V[i*n+k] = t;
	                            }
	                        }
	                        if (k < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = U[i*n+k + 1];
	                                U[i*n+k + 1] = U[i*n+k];
	                                U[i*n+k] = t;
	                            }
	                        }
	                        k++;
	                    }
	                    iter = 0;
	                    p--;
	                }
	                break;
	            }
	        }

	        // Set the small value tolerance used to calculate rank and pseudo-inverse
	        tol = FastMath.max(n * S[0] * EPS,
	                           FastMath.sqrt(Precision.SAFE_MIN));
	        
	        
	        
			end = System.currentTimeMillis();			
			System.out.println("Step 2 in " + (end-start)/1000+" seconds");

	        
	        double max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(Math.abs(U[i*n+j])- Math.abs(V[i*n+j])));
		       }
	        }
	        System.out.println("max diff = " + max);
	        max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(U[i*n+j] + V[i*n+j]));
		       }
	        }
	        System.out.println("max diff = " + max);
//        	System.out.println(Arrays.toString(singularValues));
//	        System.out.println();
//	        for (int i = 0; i < n; i++) {
//	        	System.out.println(Arrays.toString(V[i]));
//	        }
	        System.out.println("tolerance = " + tol);


	    }

	    public void singularValueDecomposition4() {
	    	System.err.println("native SVD");
			long start = System.currentTimeMillis();
			int n = nodes.size();
	        final double[] A = new double[n*n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id+n*node.id] = 1; 
				}
				A[t.id+n*t.id] = -t.neighbours.length; 
			}
			double [] R = singularValueDecompositionNativeDP(A, n);
	        U = new double[n*n];
	        V = new double[n*n];
	        S = new double[n];
	        System.arraycopy(R, 0, U, 0, n * n);
	        System.arraycopy(R, n*n, S, 0, n);
	        System.arraycopy(R, n*n + n, V, 0, n * n);

	    }
	    public void singularValueDecomposition4SP() {
			long start = System.currentTimeMillis();
			int n = nodes.size();
	        final float[] A = new float[n*n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id+n*node.id] = 1; 
				}
				A[t.id+n*t.id] = -t.neighbours.length; 
			}
			float [] R = singularValueDecompositionNativeSP(A, n);
	        U = new double[n*n];
	        V = new double[n*n];
	        S = new double[n];
	        for (int i = 0;i < n*n; i++) {
	        	U[i] = R[i];
	        }
	        for (int i = 0;i < n; i++) {
	        	S[i] = R[i + n*n];
	        }
	        for (int i = 0;i < n*n; i++) {
	        	V[i] = R[i+n*(n+1)];
	        }
	    }
	    
		/** double precision **/
		native double[] singularValueDecompositionNativeDP(double[] a, int n);
		/** single precision **/
		native float[] singularValueDecompositionNativeSP(float[] a, int n);

		/* adapted from org.apache.commons.math3.linear.SingularValueDecomposition */
	    public void singularValueDecomposition3() {
	    	System.err.println("optimised java SVD");
			long start = System.currentTimeMillis();
			int n = nodes.size();
	        final double[] A = new double[n*n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id+n*node.id] = 1; 
				}
				A[t.id+n*t.id] = -t.neighbours.length; 
			}


	         // "m" is always the largest dimension.
            //transposed = false;

	        S = new double[n];
	        U = new double[n*n];
	        V = new double[n*n];
	        final double[] e = new double[n];
	        final double[] work = new double[n];
	        // Reduce A to bidiagonal form, storing the diagonal elements
	        // in s and the super-diagonal elements in e.
	        final int nct = FastMath.min(n - 1, n);
	        final int nrt = FastMath.max(0, n - 2);
	        final int nMaxNctNrt = FastMath.max(nct, nrt);
	        for (int k = 0; k < nMaxNctNrt; k++) {
	            if (k < nct) {
	                // Compute the transformation for the k-th column and
	                // place the k-th diagonal in s[k].
	                // Compute 2-norm of k-th column without under/overflow.
	            	double svk = 0; // contains singularValues[k] 
	            	svk = 0;
	                for (int i = k; i < n; i++) {
	                	svk = FastMath.hypot(svk, A[n*k + i]);
	                }
	                if (svk != 0) {
	                    if (A[n*k + k] < 0) {
	                    	svk = -svk;
	                    }
	                    for (int i = k; i < n; i++) {
	                        A[n*k + i] /= svk;
	                    }
	                    A[n*k + k] += 1;
	                }
	                S[k] = -svk;
	            }
	            for (int j = k + 1; j < n; j++) {
	                if (k < nct && S[k] != 0) {
	                    // Apply the transformation.
	                    double t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += A[n*k + i] * A[n*j+i];
	                    }
	                    t = -t / A[k+n*k];
	                    for (int i = k; i < n; i++) {
	                        A[n*j + i] += t * A[n*k + i];
	                    }
	                }
	                // Place the k-th row of A into e for the
	                // subsequent calculation of the row transformation.
	                e[j] = A[k+n*j];
	            }
	            if (k < nct) {
	                // Place the transformation in U for subsequent back
	                // multiplication.
	                //for (int i = k; i < n; i++) {
	                //    U[n*k + i] = A[n*k + i];
	                //}
	                System.arraycopy(A, k+n*k, U, k+n*k, n - k);
	            }
	            if (k < nrt) {
	                // Compute the k-th row transformation and place the
	                // k-th super-diagonal in e[k].
	                // Compute 2-norm without under/overflow.
	                e[k] = 0;
	                for (int i = k + 1; i < n; i++) {
	                    e[k] = FastMath.hypot(e[k], e[i]);
	                }
	                if (e[k] != 0) {
	                    if (e[k + 1] < 0) {
	                        e[k] = -e[k];
	                    }
	                    for (int i = k + 1; i < n; i++) {
	                        e[i] /= e[k];
	                    }
	                    e[k + 1] += 1;
	                }
	                e[k] = -e[k];
	                if (k + 1 < n &&
	                    e[k] != 0) {
	                    // Apply the transformation.
	                    for (int i = k + 1; i < n; i++) {
	                        work[i] = 0;
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        for (int i = k + 1; i < n; i++) {
	                            work[i] += e[j] * A[n*j+i];
	                        }
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        final double t = -e[j] / e[k + 1];
	                        for (int i = k + 1; i < n; i++) {
	                            A[n*j+i] += t * work[i];
	                        }
	                    }
	                }

	                // Place the transformation in V for subsequent
	                // back multiplication.
	                //for (int i = k + 1; i < n; i++) {
	                //    V[n*k + i] = e[i];
	                //}
	                System.arraycopy(e, k+1, V, n*k + k + 1, n - k - 1);
	            }
	        }
	        
	        
			long end = System.currentTimeMillis();			
			System.out.println("Step 1 in " + (end-start)/1000+" seconds");
			start = end;
			
	        
	        // Set up the final bidiagonal matrix or order p.
	        int p = n;
	        if (nct < n) {
	            S[nct] = A[nct+n*nct];
	        }
	        if (n < p) {
	            S[p - 1] = 0;
	        }
	        if (nrt + 1 < p) {
	            e[nrt] = A[nrt+n*(p - 1)];
	        }
	        e[p - 1] = 0;

	        // Generate U.
	        for (int j = nct; j < n; j++) {
	            //for (int i = 0; i < n; i++) {
	            //    U[n*j+i] = 0;
	            //}
                Arrays.fill(U, n*j, n*j+n, 0);
	            U[j+n*j] = 1;
	        }
	        for (int k = nct - 1; k >= 0; k--) {
	            if (S[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    double t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += U[n*k + i] * U[n*j+i];
	                    }
	                    t = -t / U[k+n*k];
	                    for (int i = k; i < n; i++) {
	                        U[n*j+i] += t * U[n*k + i];
	                    }
	                }
	                for (int i = k; i < n; i++) {
	                    U[n*k + i] = -U[n*k + i];
	                }
	                U[k+n*k] = 1 + U[k+n*k];
	                //for (int i = 0; i < k - 1; i++) {
	                //    U[n*k + i] = 0;
	                //}
	                if (k > 0)
	                	Arrays.fill(U, n*k, n*k+k - 1, 0);
	            } else {
	                //for (int i = 0; i < n; i++) {
	                //    U[n*k + i] = 0;
	                //}
		            Arrays.fill(U, n*k, n*k+n, 0);
	                U[k+n*k] = 1;
	            }
	        }

	        // Generate V.
	        for (int k = n - 1; k >= 0; k--) {
	            if (k < nrt &&
	                e[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    double t = 0;
	                    for (int i = k + 1; i < n; i++) {
	                        t += V[n*k + i] * V[n*j+i];
	                    }
	                    t = -t / V[k + 1+n*k];
	                    for (int i = k + 1; i < n; i++) {
	                        V[n*j+i] += t * V[n*k + i];
	                    }
	                }
	            }
	            //for (int i = 0; i < n; i++) {
	            //    V[n*k + i] = 0;
	            //}
	            Arrays.fill(V, n*k, n*k+n, 0);
	            V[k+n*k] = 1;
	        }

	        // Main iteration loop for the singular values.
	        final int pp = p - 1;
	        int iter = 0;
	        while (p > 0) {
	            int k;
	            int kase;
	            // Here is where a test for too many iterations would go.
	            // This section of the program inspects for
	            // negligible elements in the s and e arrays.  On
	            // completion the variables kase and k are set as follows.
	            // kase = 1     if s(p) and e[k-1] are negligible and k<p
	            // kase = 2     if s(k) is negligible and k<p
	            // kase = 3     if e[k-1] is negligible, k<p, and
	            //              s(k), ..., s(p) are not negligible (qr step).
	            // kase = 4     if e(p-1) is negligible (convergence).
	            for (k = p - 2; k >= 0; k--) {
	                final double threshold
	                    = TINY + EPS * (FastMath.abs(S[k]) +
	                                    FastMath.abs(S[k + 1]));

	                // the following condition is written this way in order
	                // to break out of the loop when NaN occurs, writing it
	                // as "if (FastMath.abs(e[k]) <= threshold)" would loop
	                // indefinitely in case of NaNs because comparison on NaNs
	                // always return false, regardless of what is checked
	                // see issue MATH-947
	                if (!(FastMath.abs(e[k]) > threshold)) {
	                    e[k] = 0;
	                    break;
	                }

	            }

	            if (k == p - 2) {
	                kase = 4;
	            } else {
	                int ks;
	                for (ks = p - 1; ks >= k; ks--) {
	                    if (ks == k) {
	                        break;
	                    }
	                    final double t = (ks != p ? FastMath.abs(e[ks]) : 0) +
	                        (ks != k + 1 ? FastMath.abs(e[ks - 1]) : 0);
	                    if (FastMath.abs(S[ks]) <= TINY + EPS * t) {
	                        S[ks] = 0;
	                        break;
	                    }
	                }
	                if (ks == k) {
	                    kase = 3;
	                } else if (ks == p - 1) {
	                    kase = 1;
	                } else {
	                    kase = 2;
	                    k = ks;
	                }
	            }
	            k++;
	            // Perform the task indicated by kase.
	            switch (kase) {
	                // Deflate negligible s(p).
	                case 1: {
	                    double f = e[p - 2];
	                    e[p - 2] = 0;
	                    for (int j = p - 2; j >= k; j--) {
	                        double t = FastMath.hypot(S[j], f);
	                        final double cs = S[j] / t;
	                        final double sn = f / t;
	                        S[j] = t;
	                        if (j != k) {
	                            f = -sn * e[j - 1];
	                            e[j - 1] = cs * e[j - 1];
	                        }

	                        for (int i = 0; i < n; i++) {
	                            t = cs * V[n*j+i] + sn * V[n*(p - 1) + i];
	                            V[n*(p - 1) + i] = -sn * V[n*j+i] + cs * V[n*(p - 1) + i];
	                            V[n*j+i] = t;
	                        }
	                    }
	                }
	                break;
	                // Split at negligible s(k).
	                case 2: {
	                    double f = e[k - 1];
	                    e[k - 1] = 0;
	                    for (int j = k; j < p; j++) {
	                        double t = FastMath.hypot(S[j], f);
	                        final double cs = S[j] / t;
	                        final double sn = f / t;
	                        S[j] = t;
	                        f = -sn * e[j];
	                        e[j] = cs * e[j];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * U[n*j+i] + sn * U[n*(k - 1) + i];
	                            U[n*(k - 1) + i] = -sn * U[n*j+i] + cs * U[n*(k - 1) + i];
	                            U[n*j+i] = t;
	                        }
	                    }
	                }
	                break;
	                // Perform one qr step.
	                case 3: {
	                    // Calculate the shift.
	                    final double maxPm1Pm2 = FastMath.max(FastMath.abs(S[p - 1]),
	                                                          FastMath.abs(S[p - 2]));
	                    final double scale = FastMath.max(FastMath.max(FastMath.max(maxPm1Pm2,
	                                                                                FastMath.abs(e[p - 2])),
	                                                                   FastMath.abs(S[k])),
	                                                      FastMath.abs(e[k]));
	                    final double sp = S[p - 1] / scale;
	                    final double spm1 = S[p - 2] / scale;
	                    final double epm1 = e[p - 2] / scale;
	                    final double sk = S[k] / scale;
	                    final double ek = e[k] / scale;
	                    final double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
	                    final double c = (sp * epm1) * (sp * epm1);
	                    double shift = 0;
	                    if (b != 0 ||
	                        c != 0) {
	                        shift = FastMath.sqrt(b * b + c);
	                        if (b < 0) {
	                            shift = -shift;
	                        }
	                        shift = c / (b + shift);
	                    }
	                    double f = (sk + sp) * (sk - sp) + shift;
	                    double g = sk * ek;
	                    // Chase zeros.
	                    for (int j = k; j < p - 1; j++) {
	                        double t = FastMath.hypot(f, g);
	                        double cs = f / t;
	                        double sn = g / t;
	                        if (j != k) {
	                            e[j - 1] = t;
	                        }
	                        f = cs * S[j] + sn * e[j];
	                        e[j] = cs * e[j] - sn * S[j];
	                        g = sn * S[j + 1];
	                        S[j + 1] = cs * S[j + 1];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * V[n*j+i] + sn * V[n*(j + 1) + i];
	                            V[n*(j + 1) + i] = -sn * V[n*j+i] + cs * V[n*(j + 1) + i];
	                            V[n*j+i] = t;
	                        }
	                        t = FastMath.hypot(f, g);
	                        cs = f / t;
	                        sn = g / t;
	                        S[j] = t;
	                        f = cs * e[j] + sn * S[j + 1];
	                        S[j + 1] = -sn * e[j] + cs * S[j + 1];
	                        g = sn * e[j + 1];
	                        e[j + 1] = cs * e[j + 1];
	                        if (j < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = cs * U[n*j+i] + sn * U[n*(j + 1) + i];
	                                U[n*(j + 1) + i] = -sn * U[n*j+i] + cs * U[n*(j + 1) + i];
	                                U[n*j+i] = t;
	                            }
	                        }
	                    }
	                    e[p - 2] = f;
	                    iter++;
	                }
	                break;
	                // Convergence.
	                default: {
//	                    // Make the singular values positive.
	                	
// RRB: no need for ordering when doing exponentiation	                	
//	                    if (S[k] <= 0) {
//	                        S[k] = S[k] < 0 ? -S[k] : 0;
//
//	                        for (int i = 0; i <= pp; i++) {
//	                            V[n*k + i] = -V[n*k + i];
//	                        }
//	                    }
//	                    // Order the singular values.
//	                    while (k < pp) {
//	                        if (S[k] >= S[k + 1]) {
//	                            break;
//	                        }
//	                        double t = S[k];
//	                        S[k] = S[k + 1];
//	                        S[k + 1] = t;
//	                        if (k < n - 1) {
//	                            for (int i = 0; i < n; i++) {
//	                                t = V[i+n*(k + 1)];
//	                                V[i+n*(k + 1)] = V[n*k + i];
//	                                V[n*k + i] = t;
//	                            }
//	                        }
//	                        if (k < n - 1) {
//	                            for (int i = 0; i < n; i++) {
//	                                t = U[i+n*(k + 1)];
//	                                U[i+n*(k + 1)] = U[n*k + i];
//	                                U[n*k + i] = t;
//	                            }
//	                        }
//	                        k++;
//	                    }
	                    iter = 0;
	                    p--;
	                }
	                break;
	            }
	        }

	        // Set the small value tolerance used to calculate rank and pseudo-inverse
	        tol = FastMath.max(n * S[0] * EPS,
	                           FastMath.sqrt(Precision.SAFE_MIN));
	        
	        
	        
			end = System.currentTimeMillis();			
			System.out.println("Step 2 in " + (end-start)/1000+" seconds");

	        
	        double max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(Math.abs(U[n*j+i])- Math.abs(V[n*j+i])));
		       }
	        }
	        System.out.println("max diff = " + max);
	        max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(U[n*j+i] + V[n*j+i]));
		       }
	        }
	        System.out.println("max diff = " + max);
//        	System.out.println(Arrays.toString(singularValues));
//	        System.out.println();
//	        for (int i = 0; i < n; i++) {
//	        	System.out.println(Arrays.toString(V[i]));
//	        }
	        System.out.println("tolerance = " + tol);


	    }
		

		/* adapted from org.apache.commons.math3.linear.SingularValueDecomposition */
	    public void singularValueDecomposition3f() {
			long start = System.currentTimeMillis();
			int n = nodes.size();
	        final float[] A = new float[n*n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id+n*node.id] = 1; 
				}
				A[t.id+n*t.id] = -t.neighbours.length; 
			}


	         // "m" is always the largest dimension.
            //transposed = false;

	        singularValuesf = new float[n];
	        Uf = new float[n*n];
	        Vf = new float[n*n];
	        final float[] e = new float[n];
	        final float[] work = new float[n];
	        // Reduce A to bidiagonal form, storing the diagonal elements
	        // in s and the super-diagonal elements in e.
	        final int nct = FastMath.min(n - 1, n);
	        final int nrt = FastMath.max(0, n - 2);
	        for (int k = 0; k < FastMath.max(nct, nrt); k++) {
	            if (k < nct) {
	                // Compute the transformation for the k-th column and
	                // place the k-th diagonal in s[k].
	                // Compute 2-norm of k-th column without under/overflow.
	                singularValuesf[k] = 0;
	                for (int i = k; i < n; i++) {
	                    singularValuesf[k] = (float) FastMath.hypot(singularValuesf[k], A[n*k + i]);
	                }
	                if (singularValuesf[k] != 0) {
	                    if (A[k+n*k] < 0) {
	                        singularValuesf[k] = -singularValuesf[k];
	                    }
	                    for (int i = k; i < n; i++) {
	                        A[n*k + i] /= singularValuesf[k];
	                    }
	                    A[k+n*k] += 1;
	                }
	                singularValuesf[k] = -singularValuesf[k];
	            }
	            for (int j = k + 1; j < n; j++) {
	                if (k < nct &&
	                    singularValuesf[k] != 0) {
	                    // Apply the transformation.
	                    float t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += A[n*k + i] * A[n*j+i];
	                    }
	                    t = -t / A[k+n*k];
	                    for (int i = k; i < n; i++) {
	                        A[n*j+i] += t * A[n*k + i];
	                    }
	                }
	                // Place the k-th row of A into e for the
	                // subsequent calculation of the row transformation.
	                e[j] = A[k+n*j];
	            }
	            if (k < nct) {
	                // Place the transformation in U for subsequent back
	                // multiplication.
	                for (int i = k; i < n; i++) {
	                    Uf[n*k + i] = A[n*k + i];
	                }
	            }
	            if (k < nrt) {
	                // Compute the k-th row transformation and place the
	                // k-th super-diagonal in e[k].
	                // Compute 2-norm without under/overflow.
	                e[k] = 0;
	                for (int i = k + 1; i < n; i++) {
	                    e[k] = (float) FastMath.hypot(e[k], e[i]);
	                }
	                if (e[k] != 0) {
	                    if (e[k + 1] < 0) {
	                        e[k] = -e[k];
	                    }
	                    for (int i = k + 1; i < n; i++) {
	                        e[i] /= e[k];
	                    }
	                    e[k + 1] += 1;
	                }
	                e[k] = -e[k];
	                if (k + 1 < n &&
	                    e[k] != 0) {
	                    // Apply the transformation.
	                    for (int i = k + 1; i < n; i++) {
	                        work[i] = 0;
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        for (int i = k + 1; i < n; i++) {
	                            work[i] += e[j] * A[n*j+i];
	                        }
	                    }
	                    for (int j = k + 1; j < n; j++) {
	                        final float t = -e[j] / e[k + 1];
	                        for (int i = k + 1; i < n; i++) {
	                            A[n*j+i] += t * work[i];
	                        }
	                    }
	                }

	                // Place the transformation in V for subsequent
	                // back multiplication.
	                //for (int i = k + 1; i < n; i++) {
	                //    V[n*k + i] = e[i];
	                //}
	                System.arraycopy(e, k+1, Vf, n*k + k + 1, n - k - 1);
	            }
	        }
	        
	        
			long end = System.currentTimeMillis();			
			System.out.println("Step 1 in " + (end-start)/1000+" seconds");
			start = end;
			
	        
	        // Set up the final bidiagonal matrix or order p.
	        int p = n;
	        if (nct < n) {
	            singularValuesf[nct] = A[nct+n*nct];
	        }
	        if (n < p) {
	            singularValuesf[p - 1] = 0;
	        }
	        if (nrt + 1 < p) {
	            e[nrt] = A[nrt+n*(p - 1)];
	        }
	        e[p - 1] = 0;

	        // Generate U.
	        for (int j = nct; j < n; j++) {
	            //for (int i = 0; i < n; i++) {
	            //    U[n*j+i] = 0;
	            //}
                Arrays.fill(Uf, n*j, n*j+n, 0);
	            Uf[j+n*j] = 1;
	        }
	        for (int k = nct - 1; k >= 0; k--) {
	            if (singularValuesf[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    float t = 0;
	                    for (int i = k; i < n; i++) {
	                        t += Uf[n*k + i] * Uf[n*j+i];
	                    }
	                    t = -t / Uf[k+n*k];
	                    for (int i = k; i < n; i++) {
	                        Uf[n*j+i] += t * Uf[n*k + i];
	                    }
	                }
	                for (int i = k; i < n; i++) {
	                    Uf[n*k + i] = -Uf[n*k + i];
	                }
	                Uf[k+n*k] = 1 + Uf[k+n*k];
	                //for (int i = 0; i < k - 1; i++) {
	                //    U[n*k + i] = 0;
	                //}
	                if (k > 0)
	                	Arrays.fill(Uf, n*k, n*k+k - 1, 0);
	            } else {
	                //for (int i = 0; i < n; i++) {
	                //    U[n*k + i] = 0;
	                //}
		            Arrays.fill(Uf, n*k, n*k+n, 0);
	                Uf[k+n*k] = 1;
	            }
	        }

	        // Generate V.
	        for (int k = n - 1; k >= 0; k--) {
	            if (k < nrt &&
	                e[k] != 0) {
	                for (int j = k + 1; j < n; j++) {
	                    float t = 0;
	                    for (int i = k + 1; i < n; i++) {
	                        t += Vf[n*k + i] * Vf[n*j+i];
	                    }
	                    t = -t / Vf[k + 1+n*k];
	                    for (int i = k + 1; i < n; i++) {
	                        Vf[n*j+i] += t * Vf[n*k + i];
	                    }
	                }
	            }
	            //for (int i = 0; i < n; i++) {
	            //    V[n*k + i] = 0;
	            //}
	            Arrays.fill(Vf, n*k, n*k+n, 0);
	            Vf[k+n*k] = 1;
	        }

	        // Main iteration loop for the singular values.
	        final int pp = p - 1;
	        int iter = 0;
	        while (p > 0) {
	            int k;
	            int kase;
	            // Here is where a test for too many iterations would go.
	            // This section of the program inspects for
	            // negligible elements in the s and e arrays.  On
	            // completion the variables kase and k are set as follows.
	            // kase = 1     if s(p) and e[k-1] are negligible and k<p
	            // kase = 2     if s(k) is negligible and k<p
	            // kase = 3     if e[k-1] is negligible, k<p, and
	            //              s(k), ..., s(p) are not negligible (qr step).
	            // kase = 4     if e(p-1) is negligible (convergence).
	            for (k = p - 2; k >= 0; k--) {
	                final float threshold
	                    = TINYF + EPSF * (FastMath.abs(singularValuesf[k]) +
	                                    FastMath.abs(singularValuesf[k + 1]));

	                // the following condition is written this way in order
	                // to break out of the loop when NaN occurs, writing it
	                // as "if (FastMath.abs(e[k]) <= threshold)" would loop
	                // indefinitely in case of NaNs because comparison on NaNs
	                // always return false, regardless of what is checked
	                // see issue MATH-947
	                if (!(FastMath.abs(e[k]) > threshold)) {
	                    e[k] = 0;
	                    break;
	                }

	            }

	            if (k == p - 2) {
	                kase = 4;
	            } else {
	                int ks;
	                for (ks = p - 1; ks >= k; ks--) {
	                    if (ks == k) {
	                        break;
	                    }
	                    final float t = (ks != p ? FastMath.abs(e[ks]) : 0) +
	                        (ks != k + 1 ? FastMath.abs(e[ks - 1]) : 0);
	                    if (FastMath.abs(singularValuesf[ks]) <= TINY + EPS * t) {
	                        singularValuesf[ks] = 0;
	                        break;
	                    }
	                }
	                if (ks == k) {
	                    kase = 3;
	                } else if (ks == p - 1) {
	                    kase = 1;
	                } else {
	                    kase = 2;
	                    k = ks;
	                }
	            }
	            k++;
	            // Perform the task indicated by kase.
	            switch (kase) {
	                // Deflate negligible s(p).
	                case 1: {
	                    float f = e[p - 2];
	                    e[p - 2] = 0;
	                    for (int j = p - 2; j >= k; j--) {
	                        float t = (float) FastMath.hypot(singularValuesf[j], f);
	                        final float cs = singularValuesf[j] / t;
	                        final float sn = f / t;
	                        singularValuesf[j] = t;
	                        if (j != k) {
	                            f = -sn * e[j - 1];
	                            e[j - 1] = cs * e[j - 1];
	                        }

	                        for (int i = 0; i < n; i++) {
	                            t = cs * Vf[n*j+i] + sn * Vf[i+n*(p - 1)];
	                            Vf[i+n*(p - 1)] = -sn * Vf[n*j+i] + cs * Vf[i+n*(p - 1)];
	                            Vf[n*j+i] = t;
	                        }
	                    }
	                }
	                break;
	                // Split at negligible s(k).
	                case 2: {
	                    float f = e[k - 1];
	                    e[k - 1] = 0;
	                    for (int j = k; j < p; j++) {
	                        float t = (float) FastMath.hypot(singularValuesf[j], f);
	                        final float cs = singularValuesf[j] / t;
	                        final float sn = f / t;
	                        singularValuesf[j] = t;
	                        f = -sn * e[j];
	                        e[j] = cs * e[j];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * Uf[n*j+i] + sn * Uf[i+n*(k - 1)];
	                            Uf[i+n*(k - 1)] = -sn * Uf[n*j+i] + cs * Uf[i+n*(k - 1)];
	                            Uf[n*j+i] = t;
	                        }
	                    }
	                }
	                break;
	                // Perform one qr step.
	                case 3: {
	                    // Calculate the shift.
	                    final float maxPm1Pm2 = FastMath.max(FastMath.abs(singularValuesf[p - 1]),
	                                                          FastMath.abs(singularValuesf[p - 2]));
	                    final float scale = FastMath.max(FastMath.max(FastMath.max(maxPm1Pm2,
	                                                                                FastMath.abs(e[p - 2])),
	                                                                   FastMath.abs(singularValuesf[k])),
	                                                      FastMath.abs(e[k]));
	                    final float sp = singularValuesf[p - 1] / scale;
	                    final float spm1 = singularValuesf[p - 2] / scale;
	                    final float epm1 = e[p - 2] / scale;
	                    final float sk = singularValuesf[k] / scale;
	                    final float ek = e[k] / scale;
	                    final float b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0f;
	                    final float c = (sp * epm1) * (sp * epm1);
	                    float shift = 0;
	                    if (b != 0 ||
	                        c != 0) {
	                        shift = (float) FastMath.sqrt(b * b + c);
	                        if (b < 0) {
	                            shift = -shift;
	                        }
	                        shift = c / (b + shift);
	                    }
	                    float f = (sk + sp) * (sk - sp) + shift;
	                    float g = sk * ek;
	                    // Chase zeros.
	                    for (int j = k; j < p - 1; j++) {
	                        float t = (float) FastMath.hypot(f, g);
	                        float cs = f / t;
	                        float sn = g / t;
	                        if (j != k) {
	                            e[j - 1] = t;
	                        }
	                        f = cs * singularValuesf[j] + sn * e[j];
	                        e[j] = cs * e[j] - sn * singularValuesf[j];
	                        g = sn * singularValuesf[j + 1];
	                        singularValuesf[j + 1] = cs * singularValuesf[j + 1];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * Vf[n*j+i] + sn * Vf[n*(j + 1) + i];
	                            Vf[n*(j + 1) + i] = -sn * Vf[n*j+i] + cs * Vf[n*(j + 1) + i];
	                            Vf[n*j+i] = t;
	                        }
	                        t = (float) FastMath.hypot(f, g);
	                        cs = f / t;
	                        sn = g / t;
	                        singularValuesf[j] = t;
	                        f = cs * e[j] + sn * singularValuesf[j + 1];
	                        singularValuesf[j + 1] = -sn * e[j] + cs * singularValuesf[j + 1];
	                        g = sn * e[j + 1];
	                        e[j + 1] = cs * e[j + 1];
	                        if (j < n - 1) {
	                            for (int i = 0; i < n; i++) {
	                                t = cs * Uf[n*j+i] + sn * Uf[n*(j + 1) + i];
	                                Uf[n*(j + 1) + i] = -sn * Uf[n*j+i] + cs * Uf[n*(j + 1) + i];
	                                Uf[n*j+i] = t;
	                            }
	                        }
	                    }
	                    e[p - 2] = f;
	                    iter++;
	                }
	                break;
	                // Convergence.
	                default: {
//	                    // Make the singular values positive.
//	                    if (singularValues[k] <= 0) {
//	                        singularValues[k] = singularValues[k] < 0 ? -singularValues[k] : 0;
//
//	                        for (int i = 0; i <= pp; i++) {
//	                            V[n*k + i] = -V[n*k + i];
//	                        }
//	                    }
//	                    // Order the singular values.
//	                    while (k < pp) {
//	                        if (singularValues[k] >= singularValues[k + 1]) {
//	                            break;
//	                        }
//	                        float t = singularValues[k];
//	                        singularValues[k] = singularValues[k + 1];
//	                        singularValues[k + 1] = t;
//	                        if (k < n - 1) {
//	                            for (int i = 0; i < n; i++) {
//	                                t = V[i+n*(k + 1)];
//	                                V[i+n*(k + 1)] = V[n*k + i];
//	                                V[n*k + i] = t;
//	                            }
//	                        }
//	                        if (k < n - 1) {
//	                            for (int i = 0; i < n; i++) {
//	                                t = U[i+n*(k + 1)];
//	                                U[i+n*(k + 1)] = U[n*k + i];
//	                                U[n*k + i] = t;
//	                            }
//	                        }
//	                        k++;
//	                    }
	                    iter = 0;
	                    p--;
	                }
	                break;
	            }
	        }

	        // Set the small value tolerance used to calculate rank and pseudo-inverse
	        tol = FastMath.max(n * singularValuesf[0] * EPS,
	                           FastMath.sqrt(Precision.SAFE_MIN));
	        
	        
	        
			end = System.currentTimeMillis();			
			System.out.println("Step 2 in " + (end-start)/1000+" seconds");

	        
	        float max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(Math.abs(Uf[n*j+i])- Math.abs(Vf[n*j+i])));
		       }
	        }
	        System.out.println("max diff = " + max);
	        max = 0;
	       for (int i = 0; i < n; i++) {
		       for (int j = 0; j < n; j++) {
		    	   max = Math.max(max, Math.abs(Uf[n*j+i] + Vf[n*j+i]));
		       }
	        }
	        System.out.println("max diff = " + max);
//        	System.out.println(Arrays.toString(singularValues));
//	        System.out.println();
//	        for (int i = 0; i < n; i++) {
//	        	System.out.println(Arrays.toString(V[i]));
//	        }
	        System.out.println("tolerance = " + tol);


	    }

		/** return distance between node (x1, y1) and neighbouring node nr (x2, y2) **/
		double getDistance(int x1, int y1, int x2, int y2) {
			GraphNode v1 = grid[x1][y1];
			int depth = grid.length;
			if (x2 < 0 || x2 >= depth || y2 < 0 || y2 >= depth) {
				return Double.POSITIVE_INFINITY;
			}
			if (Math.abs(x1 - x2) > 1 || Math.abs(y1 - y2) > 1 ) {
				throw new RuntimeException("node (x1, y1) must be neighbouring node (x2, y2)");
			}
			GraphNode v2 = grid[x2][y2];
			int k = 0;
			for (GraphNode v : v1.neighbours) {
				if (v2.id == v.id) {
					return v1.distance[k];
				}
				k++;
			}
			throw new RuntimeException("node (x2, y2) is not in neighbor set of (x1, y1)");
		}
		
		public enum DISTANCE_METHOD {by_graph, parallel, by_native}; 
		
		public double[] distances(GraphNode t1, DISTANCE_METHOD method) {
			switch(method) { 
			case by_graph : 
				return super.distances(t1);
			case parallel : 
				return parallelDistances(t1);
			case by_native: 
				return nativeDistances(t1);
			}
			return null;
		}

		private double[] parallelDistances(GraphNode t1) {
			if (!allNeighborsInput.get()) {
				throw new RuntimeException("parallelDistances() assumes allNeighbours = true");
			}
			int depth = grid.length;
			double [] distances = new double[(depth+2)*(depth+2)];
			Arrays.fill(distances, 1e100);
			int id = t1.id;
			distances[(id/depth + 1) * (depth+2) + (id%depth) + 1] = 0;
			double [] pairwise = getPairWise();

			double [] distances2 = new double[(depth+2)*(depth+2)];
			System.arraycopy(distances, 0, distances2, 0, distances.length);
			
			for (int step = 0; step < 50; step++) {
				for (int i = 0; i < depth; i++) {
					for (int j = 0; j < depth; j++) {
						int k = ((i+1) * (depth + 2) + j +1);
						double d = distances[k];
						d = Math.min(d, distances[k - depth - 2] + pairwise[k*8]);
						d = Math.min(d, distances[k - depth - 1] + pairwise[k*8+1]);
						d = Math.min(d, distances[k + 1] + pairwise[k*8+2]);
						d = Math.min(d, distances[k + depth + 3] + pairwise[k*8+3]);
						d = Math.min(d, distances[k + depth + 2] + pairwise[k*8+4]);
						d = Math.min(d, distances[k + depth + 1] + pairwise[k*8+5]);
						d = Math.min(d, distances[k - 1] + pairwise[k*8+6]);
						d = Math.min(d, distances[k - depth - 3] + pairwise[k*8+7]);
						distances2[k] = d;
					}
				}
				double [] tmp = distances; 
				distances = distances2;
				distances2 = tmp;
			}
			
			
			return distances;
		}
		
		double [] getPairWise() {
			int depth = grid.length;
			double [] pairwise = new double[(depth+2)*(depth+2) * 8];
			Arrays.fill(pairwise, 1e100);
			for (int i = 0; i < depth; i++) {
				for (int j = 0; j < depth; j++) {
					int k = ((i+1) * (depth + 2) + j +1)* 8;
					pairwise[k++] = getDistance(i, j, i-1, j);
					pairwise[k++] = getDistance(i, j, i-1, j + 1);
					pairwise[k++] = getDistance(i, j, i, j+1);
					pairwise[k++] = getDistance(i, j, i+1, j + 1);
					pairwise[k++] = getDistance(i, j, i+1, j);
					pairwise[k++] = getDistance(i, j, i+1, j - 1);
					pairwise[k++] = getDistance(i, j, i, j-1);
					pairwise[k++] = getDistance(i, j, i-1, j - 1);
				}
			}
			return pairwise;
		}

		native void setDevice(int device);
		native double [] doNativeDistances(double [] distances, double [] pairwise, int depth, int iterations);
		native float [] doNativeDistancesSinglePrecission(float [] distances, float[] pairwise, int depth, int iterations);

		private double[] nativeDistances(GraphNode t1) {
			if (!allNeighborsInput.get()) {
				throw new RuntimeException("parallelDistances() assumes allNeighbours = true");
			}
			int depth = grid.length;
			double [] distances = new double[(depth+2)*(depth+2)];
			Arrays.fill(distances, 1e100);
			int id = t1.id;
			distances[(id/depth + 1) * (depth+2) + (id%depth) + 1] = 0;
			double [] pairwise = getPairWise();

			distances = doNativeDistances(distances, pairwise, depth, 50);
			
			return distances;
		}

	    public void singularValueDecomposition5() {
			long start = System.currentTimeMillis();
			int n = nodes.size();
			FlatMatrixEigenSystem es = new FlatMatrixEigenSystem(n);
	        final double[] A = new double[n*n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id+n*node.id] = 1; 
				}
				A[t.id+n*t.id] = -t.neighbours.length; 
			}
			
//			double [][] matrix2 = new double[n][n];
//			for (int i = 0; i < n; i++) {
//				for (int j = 0; j < n; j++) {
//					double sum = 0;
//					for (int k = 0; k < n; k++) {
//						sum += A[i*n+k] * A[k*n+j];
//					}
//					matrix2[i][j] = sum;
//				}
//			}
			EigenDecomposition e = es.decomposeMatrix(A);
			S = e.getEigenValues();
			System.out.println(Arrays.toString(S));
			U  = e.getEigenVectors();
	        V = e.getInverseEigenVectors();
	    }

	    public void singularValueDecomposition6() {
			long start = System.currentTimeMillis();
			int n = nodes.size();
			FlatMatrixEigenSystem2 es = new FlatMatrixEigenSystem2();
	        final double[] A = new double[n*n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id+n*node.id] = 1; 
				}
				A[t.id+n*t.id] = -t.neighbours.length; 
			}
			
//			double [][] matrix2 = new double[n][n];
//			for (int i = 0; i < n; i++) {
//				for (int j = 0; j < n; j++) {
//					double sum = 0;
//					for (int k = 0; k < n; k++) {
//						sum += A[i*n+k] * A[k*n+j];
//					}
//					matrix2[i][j] = sum;
//				}
//			}
			EigenDecomposition e = es.decomposeMatrix(A, false);
			S = e.getEigenValues();
			U  = e.getEigenVectors();
	        V = e.getInverseEigenVectors();
	    }
	
		void loadLibrary() {
			try {
				System.loadLibrary("Tesselate");
				m_bNative = true;
				//setDevice(1);
				System.err.println("Going native!!!");
			} catch (Throwable e) {
				System.err.println(e.getMessage());
				System.err.println("Going slow!!!");
			}
		}
		
		static public void getTransitionProbabilities(double distance, double[] probabilities, int nrOfStates, EigenDecomposition eigenDecomposition) {
	
		    int i, j, k;
		    double temp;
		
		    // is the following really necessary?
		    // implemented a pool of iexp matrices to support multiple threads
		    // without creating a new matrix each call. - AJD
		    // a quick timing experiment shows no difference - RRB
		    double[] iexp = new double[nrOfStates * nrOfStates];
		    // Eigen vectors
		    double[] Evec = eigenDecomposition.getEigenVectors();
		    // inverse Eigen vectors
		    double[] Ievc = eigenDecomposition.getInverseEigenVectors();
		    // Eigen values
		    double[] Eval = eigenDecomposition.getEigenValues();
		    for (i = 0; i < nrOfStates; i++) {
		        temp = Math.exp(distance * Eval[i]);
		        for (j = 0; j < nrOfStates; j++) {
		            iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
		        }
		    }
		
		    int u = 0;
		    for (i = 0; i < nrOfStates; i++) {
		        for (j = 0; j < nrOfStates; j++) {
		            temp = 0.0;
		            for (k = 0; k < nrOfStates; k++) {
		                temp += Evec[i * nrOfStates + k] * iexp[k * nrOfStates + j];
		            }
		
		            probabilities[u] = Math.abs(temp);
		            u++;
		        }
		    }
		} // getTransitionProbabilities
		
		public static void main(String[] args) throws Exception {
			final GridTesselation tessel = new GridTesselation();
			tessel.loadLibrary();
			tessel.initByName("depth", 51, "bbox", "5 120 47 154", "reddistance", "0.5", 
					"allNeighbors", true, "useGreatCircle", true);
			
			int n = tessel.nodes.size();
			Randomizer.setSeed(122L);
			final GraphNode t1 = tessel.getLowerLeftCorner();//tessel.nodes.get(Randomizer.nextInt(n));
			//double [] dist = new double [] {0.0};

			
			final double [] dist = tessel.distances(t1, DISTANCE_METHOD.by_graph);
			final double [] dist2 = tessel.distances(t1, DISTANCE_METHOD.by_native);
			double max = 0;
			int depth = tessel.grid.length;
			for (int i = 0; i < depth; i++) {
				for (int j = 0; j < depth; j++) {
					max = Math.max(max, Math.abs(dist[i * depth + j] - dist2[(i+1)*(depth+2) + j + 1]));
				}
			}
			System.err.println("max diff = " + max);

		}
}


/*
#nodes = 2500
Done in 68 seconds (native)

Step 1 in 24 seconds
Step 2 in 52 seconds
Done in 77 seconds (java)
Done in 78 seconds (java)
Max abs differnce: max|A-USV| = 1.1102230246251565E-13

Single precision
Step 1 in 23 seconds
Step 2 in 45 seconds
Done in 68 seconds
Max abs differnce: 6.056500684259447E-5

Double precision, same TINY & EPS as single precision
Step 1 in 26 seconds
Step 2 in 45 secondsDone in 72 seconds
Max abs differnce: 5.949777826601515E-9

native - large TINY
#nodes= 2500
Done in 64 seconds
Max abs differnce: 1.0835776720341528E-13



#nodes= 1024
java apache matrix=array of array SVD Done in 41.5 seconds
optimised java SVD Done in 5.4 seconds
native - large TINY Done in 3.9 seconds
EigenSystemBEAST Done in 24.4 seconds
EigenSystemApache Done in 22.1 seconds (version0) 12.3 seconds (version2) 11.6 seconds (version3)

#nodes= 1600
java apache matrix=array of array SVD 193.6 seconds
optimised java SVD Done in 19.4 seconds
native SVD in 16.8 seconds


#nodes= 576
AStep 1 1.7 seconds
AStep 2 0.3 seconds
Step 3 0.2 seconds
Done in 2.3 seconds
Max abs differnce: 4.608292913932388E-14

#nodes= 1024
AStep 1 16.7 seconds
AStep 2 1.6 seconds
Step 3 1.6 seconds
Done in 20.0 seconds (version 0) 12.1 -14.1 (version 4) 12.5 - 14.6 (version 5)

#nodes= 1600
ST 1 Done in 10.3 seconds
ST 1 Done in 2.4 seconds
ST 2 Done in 12.6 seconds
AStep 1 25.4 seconds
AStep 2 7.9 seconds
Step 3 6.7 seconds
Done in 40.0 seconds
Max abs differnce: 1.046766839873925E-13

#nodes= 2500
AStep 1 88.9 seconds
AStep 2 27.7 seconds
Step 3 27.8 seconds
Done in 144.5 148.8seconds version5 141.5 - 148.8 version4
*/