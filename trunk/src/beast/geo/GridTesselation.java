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
			singularValueDecomposition();
			//SingularValueDecomposition svd = new SingularValueDecomposition(m);
//			RealMatrix U = svd.getU();
//			RealMatrix S = svd.getS();
//			RealMatrix V = svd.getV();
			long end = System.currentTimeMillis();
			
			System.out.println("Done in " + (end-start)/1000+" seconds");
//
//			System.out.println(U);
//			System.out.println(S);
//			System.out.println(V);
			
			
			
			System.exit(0);
			
			
			
		}
	
		
	    /** Relative threshold for small singular values. */
	    private static final double EPS = 0x1.0p-52;
	    /** Absolute threshold for small singular values. */
	    private static final double TINY = 0x1.0p-966;
	    /** Computed singular values. */
	    private double[] singularValues;
	    /** row dimension = column dimension */
	    private int n;
	    /** Indicator for transposed matrix. */
	    private boolean transposed;
	    /**
	     * Tolerance value for small singular values, calculated once we have
	     * populated "singularValues".
	     **/
	    private double tol;


		/* adapted from org.apache.commons.math3.linear.SingularValueDecomposition */
	    public void singularValueDecomposition() {
			int n = nodes.size();
	        final double[][] A = new double[n][n];
			for (GraphNode t : nodes) {
				for (GraphNode node : t.neighbours) {
					A[t.id][node.id] = 1; 
				}
				A[t.id][t.id] = -t.neighbours.length; 
			}


	         // "m" is always the largest dimension.
            transposed = false;

	        singularValues = new double[n];
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
	                singularValues[k] = 0;
	                for (int i = k; i < n; i++) {
	                    singularValues[k] = FastMath.hypot(singularValues[k], A[i][k]);
	                }
	                if (singularValues[k] != 0) {
	                    if (A[k][k] < 0) {
	                        singularValues[k] = -singularValues[k];
	                    }
	                    for (int i = k; i < n; i++) {
	                        A[i][k] /= singularValues[k];
	                    }
	                    A[k][k] += 1;
	                }
	                singularValues[k] = -singularValues[k];
	            }
	            for (int j = k + 1; j < n; j++) {
	                if (k < nct &&
	                    singularValues[k] != 0) {
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
	        
	        
	        // Set up the final bidiagonal matrix or order p.
	        int p = n;
	        if (nct < n) {
	            singularValues[nct] = A[nct][nct];
	        }
	        if (n < p) {
	            singularValues[p - 1] = 0;
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
	            if (singularValues[k] != 0) {
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
	                    = TINY + EPS * (FastMath.abs(singularValues[k]) +
	                                    FastMath.abs(singularValues[k + 1]));

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
	                    if (FastMath.abs(singularValues[ks]) <= TINY + EPS * t) {
	                        singularValues[ks] = 0;
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
	                        double t = FastMath.hypot(singularValues[j], f);
	                        final double cs = singularValues[j] / t;
	                        final double sn = f / t;
	                        singularValues[j] = t;
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
	                        double t = FastMath.hypot(singularValues[j], f);
	                        final double cs = singularValues[j] / t;
	                        final double sn = f / t;
	                        singularValues[j] = t;
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
	                    final double maxPm1Pm2 = FastMath.max(FastMath.abs(singularValues[p - 1]),
	                                                          FastMath.abs(singularValues[p - 2]));
	                    final double scale = FastMath.max(FastMath.max(FastMath.max(maxPm1Pm2,
	                                                                                FastMath.abs(e[p - 2])),
	                                                                   FastMath.abs(singularValues[k])),
	                                                      FastMath.abs(e[k]));
	                    final double sp = singularValues[p - 1] / scale;
	                    final double spm1 = singularValues[p - 2] / scale;
	                    final double epm1 = e[p - 2] / scale;
	                    final double sk = singularValues[k] / scale;
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
	                        f = cs * singularValues[j] + sn * e[j];
	                        e[j] = cs * e[j] - sn * singularValues[j];
	                        g = sn * singularValues[j + 1];
	                        singularValues[j + 1] = cs * singularValues[j + 1];

	                        for (int i = 0; i < n; i++) {
	                            t = cs * V[i][j] + sn * V[i][j + 1];
	                            V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
	                            V[i][j] = t;
	                        }
	                        t = FastMath.hypot(f, g);
	                        cs = f / t;
	                        sn = g / t;
	                        singularValues[j] = t;
	                        f = cs * e[j] + sn * singularValues[j + 1];
	                        singularValues[j + 1] = -sn * e[j] + cs * singularValues[j + 1];
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
	                    if (singularValues[k] <= 0) {
	                        singularValues[k] = singularValues[k] < 0 ? -singularValues[k] : 0;

	                        for (int i = 0; i <= pp; i++) {
	                            V[i][k] = -V[i][k];
	                        }
	                    }
	                    // Order the singular values.
	                    while (k < pp) {
	                        if (singularValues[k] >= singularValues[k + 1]) {
	                            break;
	                        }
	                        double t = singularValues[k];
	                        singularValues[k] = singularValues[k + 1];
	                        singularValues[k + 1] = t;
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
	        tol = FastMath.max(n * singularValues[0] * EPS,
	                           FastMath.sqrt(Precision.SAFE_MIN));

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
		
		public static void main(String[] args) throws Exception {
			final GridTesselation tessel = new GridTesselation();
			tessel.loadLibrary();
			tessel.initByName("depth", 4, "bbox", "5 120 47 154", "reddistance", "1.0", 
					"allNeighbors", false, "useGreatCircle", true);
			
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
