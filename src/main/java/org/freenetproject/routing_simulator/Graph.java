package org.freenetproject.routing_simulator;

import java.io.*;
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

	public interface LinkLengthSource {
		/**
		 * @return a desired link length for a connection determined by the link length distribution scheme.
		 * This will be attempted to be matched as closely as location distribution allows.
		 */
		public double getLinkLength(Random random);
	}

	/**
	 * Generates link lengths with uniform / flat probability. Terrible distribution.
	 */
	public static class UniformLinkSource implements LinkLengthSource {
		@Override
		public double getLinkLength(Random random) {
			return random.nextDouble() * 0.5;
		}
	}

	public static class ConformingLinkSource implements LinkLengthSource {
		private final ArrayList<Double> lengths;
		public ConformingLinkSource(String filename) {
			lengths = new ArrayList<Double>();
			try {
				BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
				//TODO: Read all, put into ArrayList, get a link length selects from that.
				String line;
				//TODO: This seems like a C++ way of doing things. What's the Java way?
				while ( (line = reader.readLine()) != null) {
					//File format has link length as first value, separated by a space.
					lengths.add(Double.valueOf(line.split(" ")[0]));
				}
			} catch (FileNotFoundException e) {
				System.out.println(e);
				System.out.println("Unable to open file \"" + filename + "\".");
				System.exit(1);
			} catch (IOException e) {
				System.out.println(e);
				System.exit(2);
			}
		}

		@Override
		public double getLinkLength(Random random) {
			return lengths.get(random.nextInt(lengths.size()));
		}
	}

	public interface DegreeSource {
		/**
		 * @return degree conforming to the distribution.
		 */
		public int getDegree();
	}

	public static class FixedDegreeSource implements DegreeSource {
		private final int degree;

		public FixedDegreeSource(int degree) {
			this.degree = degree;
		}

		@Override
		public int getDegree() {
			return degree;
		}
	}

	public static class ConformingDegreeSource implements DegreeSource {
		private final WeightedDistribution distribution;

		public ConformingDegreeSource(String filename, Random random) {
			this.distribution  = new WeightedDistribution(filename, random);
		}

		@Override
		public int getDegree() {
			return distribution.randomValue();
		}
	}

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
		if (param.fastGeneration) {
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

	private static class DistanceEntry implements Comparable<DistanceEntry> {
		public final double distance;
		public final int index;

		public DistanceEntry(double distance, int index) {
			this.distance = distance;
			this.index = index;
		}

		@Override
		public int compareTo(DistanceEntry other) {
			return Double.compare(this.distance, other.distance);
		}
	}

	//TODO: Using GraphParam as an argument is beginning to smell: this takes additional arguments and ignores the number of close connections.
	/**
	 * Generates a graph with link length distribution and peer count distribution as described in the given sources.
	 * @param param graph generation parameters. Local and remote connections irrelevant as given distribution is followed.
	 * @param rand used for random numbers
	 * @return specified graph
	 */
	public static Graph generateGraph(GraphParam param, Random rand, DegreeSource degreeSource, LinkLengthSource linkLengthSource) {
		Graph g = new Graph(param.n);
		g.generateNodes(param, rand);
		for (SimpleNode node : g.nodes) {
			g.nodes.set(node.index, new WeightedDegreeNode(node.getLocation(), node.lowUptime(), param.pInstantReject, rand, degreeSource.getDegree()));
			g.nodes.get(node.index).index = node.index;
		}

		//Probability of not making a connection with a peer which has its desired degree.
		final double rejectProbability = 0.98;
		//TODO: Some way to get desired peer distributions cleanly? This is copy-paste from generate1dKleinbergGraph because it needs to drop out mid-loop.
		//make far links
		DistanceEntry[] distances = new DistanceEntry[param.n];
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
					distances[j] = new DistanceEntry(Location.distance(src.getLocation(), g.nodes.get(j).getLocation()), j);
				}

				//System.out.println("distances size is " + );
				Arrays.sort(distances);

				//Make q distant connections
				while (!src.atDegree()) {
					double length = linkLengthSource.getLinkLength(rand);
					int idx = Arrays.binarySearch(distances, new DistanceEntry(length, -1));
					if (idx < 0) idx = -1 - idx;
					if (idx >= param.n) idx = param.n - 1;
					dest = (WeightedDegreeNode)g.nodes.get(distances[idx].index);
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
					/*
					 * The array is sorted by location and evenly spaced, so a change in index goes
					 * a consistent distance away.
					 */
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

				/*
				 * Find normalizing constant for this node - sum distance probabilities so that they
				 * they are in increasing order and may be searched through to find the closest link.
				 * Note that this means here the probability is proportional to 1/distance.
				 * sumProb is a non-normalized CDF of probabilities by node index.
				 */
				double norm = 0.0;
				for (int j = 0; j < n; j++) {
					if (i != j) {
						norm += 1.0 / Location.distance(src.getLocation(), g.nodes.get(j).getLocation());
					}
					sumProb[j] = norm;
					//CDF must be non-decreasing
					if (j > 0) assert sumProb[j] >= sumProb[j-1];
				}

				/* Make q distant connections - */
				for (int j = 0; j < q; j++) {
					/*
					 * sumProb is a CDF, so to weight by it pick a "Y value" and find closest index.
					 * norm is now the highest (and last) value in the CDF, so this is picking
					 * a distance probability sum and finding the closest node for that distance.
					 * Because there are more nodes which match values in highly represented domains
					 * (steeper in the CDF) a random value is more likely to be in those areas.
					 */
					double x = rand.nextDouble() * norm;
					assert x <= norm;
					int idx = Arrays.binarySearch(sumProb, x);
					/*
					 * If such value is not actually present, as it might not be due to being
					 * floating point, use the index where it would be inserted:
					 * idx = -insertion point - 1
					 * insertion point = -1 - idx
					 * The insertion point would be the length of the array and thus out of bounds
					 * if all elements were less than it, but this will not happen as norm is the
					 * greatest element and nextDouble() is [0, 1). This does not mean it will not
					 * choose the greatest element as insertion point is the index of the first
					 * greater element.
					 * TODO: Does this result in choosing the greater value when the lesser is closer? Looks like it yes.
					 */
					if (idx < 0) idx = -1 - idx;
					//idx is index of the first greater element, but use the lesser if it is closer.
					if (idx > 0 && Math.abs(x - sumProb[idx - 1]) < Math.abs(x - sumProb[idx])) idx--;
					//Assert that this actually is the closest.
					if (idx > 0) assert Math.abs(x - sumProb[idx]) < Math.abs(x - sumProb[idx - 1]);
					if (idx < sumProb.length - 1) assert Math.abs(x - sumProb[idx]) < Math.abs(x - sumProb[idx + 1]);
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
	 * Writes graph to a file. Format:
	 * <ul>
	 *      <li>Number of nodes.</li>
	 *      <li>Serialized SimpleNodes.</li>
	 * </ul>
	 * @param destination file to write graph to.
	 */
	public void write(File destination) {
		try {
			final FileOutputStream outputStream = new FileOutputStream(destination);
			final ObjectOutputStream output = new ObjectOutputStream(outputStream);

			// Number of nodes.
			output.writeInt(nodes.size());

			// Nodes.
			for (SimpleNode node : nodes) output.writeObject(node);

			/*
			 * Write connections starting from zero index to higher index nodes. This way each connection
			 * will only be written once: the other end at a higher index will not write the connection to
			 * the lower node, which has already been written.
			 *
			 * Add to intermediate list first in order to be able to write the number of connections for the
			 * purposes of reading more easily.
			 * TODO: Better to read pairs until error? Less extensible.
			 */
			final ArrayList<Integer> connectionIndexes = new ArrayList<Integer>();
			int writtenConnections = 0;
			for (int i = 0; i < nodes.size(); i++) {
				for (SimpleNode node: nodes.get(i).getConnections()) {
					if (node.index > i) {
						writtenConnections++;
						connectionIndexes.add(i);
						connectionIndexes.add(node.index);
					}
				}
			}

			output.writeInt(writtenConnections);
			System.out.println("Writing " + writtenConnections + " connections.");
			assert connectionIndexes.size() == writtenConnections * 2;
			for (Integer index : connectionIndexes) {
				output.writeInt(index);
			}

			output.flush();
			outputStream.close();
		} catch (IOException e) {
			System.err.println("Could not write to " + destination.getAbsolutePath() + ":");
			e.printStackTrace();
			System.exit(3);
		}
	}

	/**
	 * Constructs the graph from a file.
	 * @param source file to read the graph from.
	 * @param random Randomness source to give to nodes.
	 * @return graph defined by the file.
	 */
	public static Graph read(File source, Random random) {
		try {
			final FileInputStream inputStream = new FileInputStream(source);
			final ObjectInputStream input = new ObjectInputStream(inputStream);

			// Number of nodes.
			final int networkSize = input.readInt();
			final Graph graph = new Graph(networkSize);
			graph.locations = new double[networkSize];

			// Nodes.
			for (int i = 0; i < networkSize; i++) {
				SimpleNode node = (SimpleNode)input.readObject();
				node.setRand(random);
				node.index = i;
				graph.locations[i] = node.getLocation();
				graph.nodes.add(node);
			}

			final int writtenConnections = input.readInt();
			System.out.println("Reading " + writtenConnections + " connections.");
			//Each connection consists of two indexes in a pair.
			for (int i = 0; i < writtenConnections; i++) {
				final int from = input.readInt();
				final int to = input.readInt();
				//System.out.println(from + " " + to);
				graph.nodes.get(from).connect(graph.nodes.get(to));
			}

			return graph;
		} catch (IOException e) {
			System.err.println("Could not read from " + source.getAbsolutePath() + ":");
			e.printStackTrace();
			System.exit(4);
			return null;
		} catch (ClassNotFoundException e) {
			System.err.println("Unexpected class in graph:");
			e.printStackTrace();
			System.exit(5);
			return null;
		}
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
			Graph g = generate1dKleinbergGraph(new GraphParam(nNodes, p, q, 0.0, 0.0, true), rand);
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
