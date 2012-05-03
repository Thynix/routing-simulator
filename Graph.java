import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Class to represent and generate graphs of small world networks.
 * At present only limited Kleinberg graphs are generated.
 * Functions to evaluate the graph topology are also provided.
 */
public class Graph {
	private ArrayList<SimpleNode> nodes;
	private double[] locations;

	/**
	 * Private constructor; call one of the generator functions instead.
	 *
	 * @param nNodes Initial internal capacity.  The actual graph may have
	 * more or fewer nodes.
	 */
	private Graph(int nNodes) {
		if (nNodes <= 0)
			throw new IllegalArgumentException("Must have positive nodes.");
		nodes = new ArrayList<SimpleNode>(nNodes);
		locations = null;
	}

	private void generateNodes(GraphParam param, Random rand) {
		locations = new double[param.n];
		if (param.evenSpacing) {
			for (int i = 0; i < param.n; i++) locations[i] = (1.0 * i) / param.n;
		} else {
			for (int i = 0; i < param.n; i++) locations[i] = rand.nextDouble();
		}

		Arrays.sort(locations);
		for (int i = 0; i < param.n; i++) {
			SimpleNode node = new SimpleNode(locations[i], rand.nextDouble() < param.pLowUptime, param.pInstantReject, rand);
			node.index = i;
			nodes.add(node);
		}

	}

	//TODO: Using GraphParam as an argument is beginning to smell: this takes additional arguments and ignores the number of close connections.
	/**
	 * Generates a graph with only long connections and a peer count distribution as described in the given file.
	 * @param param graph generation parameters. Local and remote connections irrelevant as given distribution is followed.
	 * @param rand used for random numbers
	 * @param filename Desired occurrences in format of "[number of peers] [occurrences]" on each line.
	 * @param forceSize If true, the size in param is used. If not, the sum of occurrences in the distribution file.
	 * @return specified graph
	 */
	public static Graph generatePeerDistGraph(GraphParam param, Random rand, String filename, boolean forceSize) {
		WeightedDistribution distribution = new WeightedDistribution(filename, new Random(rand.nextLong()));
		param = new GraphParam(forceSize ? param.n : distribution.totalOccurances, param.p, param.q, param.pLowUptime, param.pInstantReject, param.evenSpacing, param.fastGeneration);
		Graph g = new Graph(param.n);
		g.generateNodes(param, rand);
		for (SimpleNode node : g.nodes) {
			g.nodes.set(node.index, new WeightedDegreeNode(node.getLocation(), node.lowUptime(), param.pInstantReject, rand, distribution));
			g.nodes.get(node.index).index = node.index;
		}

		//Probability of not making a connection with a peer which has its desired degree.
		final double rejectProbability = 0.98;
		//TODO: Some way to get desired peer distributions cleanly? This is copy-paste from below because it needs to drop out mid-loop.
		//make far links
		double[] sumProb = new double[param.n];
		for (int i = 0; i < param.n; i++) {
			WeightedDegreeNode src = (WeightedDegreeNode) g.nodes.get(i);
			if (src.atDegree()) continue;
			WeightedDegreeNode dest;
			if (param.fastGeneration) {
				//Continuous approximation to 1/d distribution; accurate in the large n case.
				//Treats spacing as even, whether or not that is accurate.
				//Assumes nodes are sorted in location order.
				double maxSteps = param.n / 2.0;
				while (!src.atDegree()) {
					int steps = (int)Math.round(Math.pow(maxSteps, rand.nextDouble()));
					assert steps >= 0 && steps <= param.n / 2;
					int idx = rand.nextBoolean() ? i + steps : i - steps;
					if (idx < 0) idx += param.n;
					if (idx >= param.n) idx -= param.n;
					dest = (WeightedDegreeNode)g.nodes.get(idx);
					/* Reject if already connected, connection to self, if with high probability if
					 * candidate peer has its desired degree
					 */
					if (idx == i || src.isConnected(dest) ||
					    (dest.atDegree() && rand.nextDouble() < rejectProbability)) continue;
					src.connect(dest);
				}
			} else {
				//Slow generation operates on exact node locations (even or otherwise).
				//Does not require sorted node order.
				//Precisely accurate even with uneven locations and small graph size.

				//Find normalizing constant for this node
				double norm = 0.0;
				for (int j = 0; j < param.n; j++) {
					if (i != j) {
						norm += 1.0 / Location.distance(src.getLocation(), g.nodes.get(j).getLocation());
					}
					sumProb[j] = norm;
					//TODO: That norm increased? Wouldn't it only decrease if distance is negative?
					if (j > 0) assert sumProb[j] >= sumProb[j-1];
				}

				//Make q distant connections
				while (!src.atDegree()) {
					double x = rand.nextDouble() * norm;
					int idx = Arrays.binarySearch(sumProb, x);
					if (idx < 0) idx = -1 - idx;
					dest = (WeightedDegreeNode)g.nodes.get(idx);
					if (src == dest || src.isConnected(dest) ||
					    (dest.atDegree() && rand.nextDouble() < rejectProbability)) continue;
					src.connect(dest);
				}
			}
		}

		return g;
	}

	/**
	 * Generate a one-dimensional Kleinberg Graph with given parameters.
	 * See The Small-World Phenomenon: An Algorithmic Perspective
	 * Jon Kleinberg, 1999
	 * http://www.cs.cornell.edu/home/kleinber/swn.pdf
	 * We use a modified version where edges are not directed.
	 * Note that q specifies outgoing links, and p is distance not link
	 * count.  Average node degree will be 2 * (p + q), minimum node
	 * degree 2 * p + q.  There is no maximum node degree.
	 * TODO: Adjacent link support not yet tested.
	 *
	 * @param param Contains graph parameters such as size and number of various-distance connections.
	 * @param rand Random number source used for initialization: locations and probabilities.
	 * @return A Graph with the desired structure
	 */
	public static Graph generate1dKleinbergGraph(GraphParam param, Random rand) {
		//TODO: Arguments list is cleaner, but this is a mess.
		final int n = param.n;
		final int q = param.q;
		final int p = param.p;
		final double pLowUptime = param.pLowUptime;
		final double pInstantReject = param.pInstantReject;
		final boolean evenSpacing = param.evenSpacing;
		final boolean fastGeneration = param.fastGeneration;

		Graph g = new Graph(n);

		//make nodes
		g.generateNodes(param, rand);

		//make adjacent links
		for (int i = 0; i < n; i++) {
			for (int j = 1; j <= p; j++) {
				//Each node connects to the p nodes after it,
				//and will be connected to by the p nodes before it.
				int target = (i + j) % n;
				g.nodes.get(i).connect(g.nodes.get(target));
			}
		}

		//make far links
		double[] sumProb = new double[n];
		for (int i = 0; i < n; i++) {
			SimpleNode src = g.nodes.get(i);
			SimpleNode dest = null;
			if (fastGeneration) {
				//Continuous approximation to 1/d distribution; accurate in the large n case.
				//Treats spacing as even, whether or not that is accurate.
				//Assumes nodes are sorted in location order.
				double maxSteps = n / 2.0;
				for (int j = 0; j < q; j++) {
					int steps = (int)Math.round(Math.pow(maxSteps, rand.nextDouble()));
					assert steps >= 0 && steps <= n / 2;
					int idx = rand.nextBoolean() ? i + steps : i - steps;
					if (idx < 0) idx += n;
					if (idx >= n) idx -= n;
					dest = g.nodes.get(idx);
					if (idx == i || src.isConnected(dest)) {
						j--;
						continue;
					}
					src.connect(dest);
				}
			} else {
				//Slow generation operates on exact node locations (even or otherwise).
				//Does not require sorted node order.
				//Precisely accurate even with uneven locations and small graph size.

				//Find normalizing constant for this node
				double norm = 0.0;
				for (int j = 0; j < n; j++) {
					if (i != j) {
						norm += 1.0 / Location.distance(src.getLocation(), g.nodes.get(j).getLocation());
					}
					sumProb[j] = norm;
					if (j > 0) assert sumProb[j] >= sumProb[j-1];
				}

				//Make q distant connections
				for (int j = 0; j < q; j++) {
					double x = rand.nextDouble() * norm;
					int idx = Arrays.binarySearch(sumProb, x);
					if (idx < 0) idx = -1 - idx;
					dest = g.nodes.get(idx);
					if (src == dest || src.isConnected(dest)) {
						j--;
						continue;
					}
					src.connect(dest);
				}
			}
		}

		return g;
	}

	/**
	 * Get a node by index.
	 *
	 * @param i Index of node to get
	 * @return Node at index i
	 */
	public SimpleNode getNode(int i) {
		return nodes.get(i);
	}

	/**
	 * Print some topology statistics.
	 */
	public void printGraphStats(boolean verbose) {
		if (verbose) {
			int nEdges = nEdges();
			double meanDegree = ((double)(2 * nEdges)) / nodes.size();
			System.out.println(	"Graph stats:");
			System.out.println(	"Size:					" + size());
			System.out.println(	"Edges:					" + nEdges);
			System.out.println(	"Min degree:				" + minDegree());
			System.out.println(	"Max degree:				" + maxDegree());
			System.out.println(	"Mean degree:				" + meanDegree);
			System.out.println(	"Degree stddev:				" + Math.sqrt(degreeVariance()));
			System.out.println(	"Mean local clustering coefficient:	" + meanLocalClusterCoeff());
			System.out.println(	"Global clustering coefficient:		" + globalClusterCoeff());
			System.out.println();
		} else {
			double[] cc = localClusterCoeff();
			int[] deg = degrees();
			ArrayStats ccStats = new ArrayStats(cc);
			ArrayStats degStats = new ArrayStats(deg);
			System.out.print(size() + "\t" + nEdges() + "\t" + minDegree() + "\t" + maxDegree() + "\t" + globalClusterCoeff() + "\t");
			System.out.print(ccStats.mean() + "\t" + ccStats.stdDev() + "\t" + ccStats.skewness() + "\t" + ccStats.kurtosis() + "\t");
			System.out.print(degStats.mean() + "\t" + degStats.stdDev() + "\t" + degStats.skewness() + "\t" + degStats.kurtosis() + "\t");
		}
	}

	/**Print column headers for printGraphStats(false).*/
	public static void printGraphStatsHeader() {
		System.out.print("nNodes\tnEdges\tminDegree\tmaxDegree\tglobalClusterCoeff\tlocalCCMean\tlocalCCStdDev\tlocalCCSkew\tlocalCCKurtosis\t");
		System.out.print("degreeMean\tdegreeStdDev\tdegreeSkew\tdegreeKurtosis\t");
	}

	/**Get the topology stats as an array.*/
	public double[] graphStats() {
		double[] cc = localClusterCoeff();
		int[] deg = degrees();
		ArrayStats ccStats = new ArrayStats(cc);
		ArrayStats degStats = new ArrayStats(deg);
		return new double[] {size(), nEdges(), minDegree(), maxDegree(), globalClusterCoeff(),
			ccStats.mean(), ccStats.stdDev(), ccStats.skewness(), ccStats.kurtosis(),
			degStats.mean(), degStats.stdDev(), degStats.skewness(), degStats.kurtosis(),
		};
	}

	/**Edge length distribution*/
	public double[] edgeLengths() {
		int nEdges = nEdges();
		double[] lengths = new double[nEdges];
		int e = 0;
		for (int i = 0; i < nodes.size(); i++) {
			SimpleNode node = nodes.get(i);
			assert node.index == i;
			ArrayList<SimpleNode> conn = node.getConnections();
			for (int j = 0; j < conn.size(); j++) {
				if (conn.get(j).index < i) continue;
				assert conn.get(j).index != i;
				lengths[e++] = node.distanceToLoc(conn.get(j).getLocation());
			}
		}
		assert e == nEdges;
		return lengths;
	}

	/**
	 * Get the number of nodes in this graph.
	 *
	 * @return Size of the graph
	 */
	public int size() {
		return nodes.size();
	}

	/**
	 * Count edges in this graph.
	 *
	 * @return Total number of edges
	 */
	public int nEdges() {
		int nEdges = 0;
		for (int i = 0; i < nodes.size(); i++) {
			nEdges += nodes.get(i).degree();
		}

		//edges are undirected and will be double counted
		assert nEdges % 2 == 0;
		return nEdges / 2;
	}

	/**
	 * Find the minimum degree of any node in the graph.
	 *
	 * @return Minimum node degree
	 */
	public int minDegree() {
		int n = size();
		if (n == 0) return 0;
		int min = nodes.get(0).degree();
		for (int i = 1; i < n; i++) {
			min = Math.min(min, nodes.get(i).degree());
		}
		return min;
	}

	/**
	 * Find the maximum degree of any node in the graph.
	 *
	 * @return Maximum node degree
	 */
	public int maxDegree() {
		int n = size();
		if (n == 0) return 0;
		int max = nodes.get(0).degree();
		for (int i = 1; i < n; i++) {
			max = Math.max(max, nodes.get(i).degree());
		}
		return max;
	}

	/**
	 * Compute the variance in the degree of the nodes.
	 *
	 * @return Variance of node degree
	 */
	public double degreeVariance() {
		long sumDegrees = 0;
		long sumSquareDegrees = 0;
		long n = nodes.size();
		if (n == 0) return 0;
		for (int i = 0; i < n; i++) {
			int d = nodes.get(i).degree();
			sumDegrees += d;
			sumSquareDegrees += d * d;
		}

		double var = ((double)sumSquareDegrees)/((double)n) - ((double)(sumDegrees * sumDegrees))/((double)(n * n));
		return var;
	}


	/**
	 * Calculate the mean clustering coefficients of nodes in the graph.
	 * See http://en.wikipedia.org/wiki/Clustering_coefficient
	 * This is *not* the same as the global clustering coefficient
	 * described there; this is the unweighted mean of the local
	 * coefficients, which gives undue weight to low-degree nodes.
	 *
	 * @return Mean local clustering coefficient
	 */
	public double meanLocalClusterCoeff() {
		double sumCoeff = 0.0;
		int n = nodes.size();
		if (n == 0) return 0;
		for (int i = 0; i < n; i++) {
			sumCoeff += nodes.get(i).localClusterCoeff();
		}
		double mean = sumCoeff / n;
		assert mean >= 0.0 && mean <= 1.0;
		return mean;
	}

	public double[] localClusterCoeff() {
		int n = nodes.size();
		double[] cc = new double[n];
		for (int i = 0; i < n; i++) cc[i] = nodes.get(i).localClusterCoeff();
		return cc;
	}

	public int[] degrees() {
		int n = nodes.size();
		int[] d = new int[n];
		for (int i = 0; i < n; i++) d[i] = nodes.get(i).degree();
		return d;
	}

	/**
	 * Calculate the global clustering coefficient.
	 * See http://en.wikipedia.org/wiki/Clustering_coefficient
	 *
	 * @return Global clustering coefficient
	 */
	public double globalClusterCoeff() {
		int nClosed = 0;
		int nTotal = 0;

		for (int i = 0; i < nodes.size(); i++) {
			SimpleNode n = nodes.get(i);
			int degree = n.degree();
			nClosed += n.closedTriplets();
			nTotal += (degree * (degree - 1)) / 2;
		}

		return ((double)(nClosed)) / ((double)(nTotal));
	}

	/**
	 * Determine whether a request was routed to its optimal location.
	 *
	 * @param r Request to check for precise routing
	 * @return Whether the request was precisely routed
	 */
	public boolean preciseRoute(Request r) {
		//Even spacing only version:
		//return r.minDist() <= (1.0 / (nodes.size() * 2.0));
		int idx = Arrays.binarySearch(locations, r.getLocation());
		if (idx < 0) idx = -1 - idx;
		int pre = idx - 1;
		if (pre == -1) pre = locations.length - 1;
		if (idx == locations.length) idx = 0;
		return r.minDist() <= Location.distance(r.getLocation(), locations[pre])
			&& r.minDist() <= Location.distance(r.getLocation(), locations[idx]);
	}

	/**
	 * Randomize the locations of the nodes prior to swapping etc.
	 *
	 * @param rand The source of randomness to use
	 * @param evenSpacing Whether the locations should be evenly spaced
	 */
	public void randomizeLocations(Random rand, boolean evenSpacing) {
		double[] newLocations = new double[size()];
		if (evenSpacing) {
			for (int i = 0; i < newLocations.length; i++) {
				newLocations[i] = ((double)(i)) / newLocations.length;
			}
			shuffle(newLocations, rand);
		} else {
			for (int i = 0; i < newLocations.length; i++) {
				newLocations[i] = rand.nextDouble();
			}
		}
		setLocations(newLocations);
	}

	/**
	 * Set locations for the nodes in this Graph.
	 *
	 * @param newLocations The set of new locations
	 */
	public void setLocations(double[] newLocations) {
		if (newLocations.length != size())
			throw new IllegalArgumentException("Wrong location count.");
		for (int i = 0; i < newLocations.length; i++) {
			nodes.get(i).setLocation(newLocations[i]);
		}
	}

	/**
	 * Knuth shuffle an array.
	 *
	 * @param a The array to shuffle
	 * @param rand The randomness source to use
	 */
	public static void shuffle(double[] a, Random rand) {
		for (int i = 0; i < a.length; i++) {
			int j = i + rand.nextInt(a.length - i);
			if (i == j) continue;
			double t = a[i];
			a[i] = a[j];
			a[j] = t;
		}
	}

	/**
	 * Perform the darknet location swapping algorithm repeatedly.
	 *
	 * @param nAttempts Number of swaps to attempt
	 * @param uniform Whether to use centralized uniform probabilities or decentralized walks
	 * @param walkDist Number of hops to walk if using decentralized walks
	 * @param uniformWalk Whether to walk uniformly or attempt to correct for high degree node bias
	 * @param rand Randomness source
	 * @return Number of swap requests accepted
	 */
	public int darknetSwap(int nAttempts, boolean uniform, int walkDist, boolean uniformWalk, Random rand) {
		int nAccepted = 0;
		for (int i = 0; i < nAttempts; i++) {
			SimpleNode origin;
			SimpleNode target;
			origin = nodes.get(rand.nextInt(nodes.size()));
			if (uniform) {
				//centralized uniform swapping -- choose both nodes from a flat distribution
				target = nodes.get(rand.nextInt(nodes.size()));
			} else {
				//decentralized random walk
				int hops = walkDist;
				if (!uniformWalk) hops *= 2;	//correct for fact that some hops get rejected
				target = origin.randomWalk(hops, uniformWalk, rand);
			}
			if (origin.attemptSwap(target)) nAccepted++;
		}
		return nAccepted;
	}

	public int[] randomWalkDistTest(int nWalks, int hopsPerWalk, boolean uniform, Random rand) {
		int[] choiceFreq = new int[size()];
		int dupCount = 0;
		for (int i = 0; i < nWalks; i++) {
			SimpleNode origin = nodes.get(rand.nextInt(size()));
			SimpleNode dest = origin.randomWalk(hopsPerWalk, uniform, rand);
			choiceFreq[dest.index]++;
			if (origin == dest) dupCount++;
		}
		System.out.println("Origin selected as dest on " + dupCount + " walks out of " + nWalks);
		return choiceFreq;
	}

	/**
	 * Create some graphs; test them for statistical properties of interest.
	 */
	public static void main(String[] args) {
		int nNodes = 4000;
		int p = 0;
		int q = 6;

		int nWalks = 10 * 1000 * 1000;
		int nBuckets = 400;
		int hopsPerWalkUniform = 20;
		int hopsPerWalkCorrected = 40;

		int nTrials = 4;
		int[][][] pdfs = new int[nTrials][3][nBuckets];

		for (int trial = 0; trial < nTrials; trial++) {
			System.out.println("Creating test graph...");
			Random rand = new MersenneTwister(trial);
			Graph g = generate1dKleinbergGraph(new GraphParam(nNodes, p, q, 0.0, 0.0, true, false), rand);
			g.printGraphStats(true);
			int[] uniformWalkDist;
			int[] weightedWalkDist;
			int[] refDist = new int[nNodes];
			System.out.println("Computing reference distribution...");
			for (int i = 0; i < nWalks; i++) refDist[rand.nextInt(nNodes)]++;
			System.out.println("Computing uniform walks...");
			uniformWalkDist = g.randomWalkDistTest(nWalks, hopsPerWalkUniform, true, rand);
			System.out.println("Computing weighted walks...");
			weightedWalkDist = g.randomWalkDistTest(nWalks, hopsPerWalkCorrected, false, rand);

			Arrays.sort(refDist);
			Arrays.sort(uniformWalkDist);
			Arrays.sort(weightedWalkDist);
			int nodesPerBucket = nNodes / nBuckets;
			assert nBuckets * nodesPerBucket == nNodes;

			for (int i = 0; i < nNodes; i++) {
				pdfs[trial][0][i / nodesPerBucket] += refDist[i];
				pdfs[trial][1][i / nodesPerBucket] += uniformWalkDist[i];
				pdfs[trial][2][i / nodesPerBucket] += weightedWalkDist[i];
			}
		}

		System.out.println("Distribution PDFs:");
		for (int i = 0; i < nTrials; i++) {
			System.out.print("Reference\tUniform\tWeighted\t");
		}
		System.out.println();
		for (int i = 0; i < nBuckets; i++) {
			for (int trial = 0; trial < nTrials; trial++) {
				System.out.print(pdfs[trial][0][i] + "\t" + pdfs[trial][1][i] + "\t" + pdfs[trial][2][i] + "\t");
			}
			System.out.println();
		}
	}
}
