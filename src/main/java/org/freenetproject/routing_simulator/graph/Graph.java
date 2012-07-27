package org.freenetproject.routing_simulator.graph;

import org.freenetproject.routing_simulator.graph.degree.PoissonDegreeSource;
import org.freenetproject.routing_simulator.graph.linklength.KleinbergLinkSource;
import org.freenetproject.routing_simulator.graph.linklength.LinkLengthSource;
import org.freenetproject.routing_simulator.graph.degree.DegreeSource;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;
import org.freenetproject.routing_simulator.util.ArrayStats;
import org.freenetproject.routing_simulator.util.MersenneTwister;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
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

	/**
	 * Probability of not making a connection with a peer which has its desired degree.
	 */
	private static final double rejectProbability = 0.98;

	/**
	 * Private constructor; call one of the generator functions instead.
	 *
	 * @param nodes The nodes which make up the network.
	 * @see
	 */
	private Graph(ArrayList<SimpleNode> nodes) {
		this.nodes = nodes;
	}

	public static ArrayList<SimpleNode> generateNodes(final int nNodes, final Random rand, boolean fastGeneration, DegreeSource source) {
		double[] locations = new double[nNodes];
		if (fastGeneration) {
			for (int i = 0; i < nNodes; i++) locations[i] = (1.0 * i) / nNodes;
		} else {
			for (int i = 0; i < nNodes; i++) locations[i] = rand.nextDouble();
		}

		// Increasing index should also mean increasing location.
		Arrays.sort(locations);

		final ArrayList<SimpleNode> nodes = new ArrayList<SimpleNode>(nNodes);

		for (int i = 0; i < nNodes; i++) {
			//TODO: index in constructor
			SimpleNode node = new SimpleNode(locations[i], rand, source.getDegree());
			node.index = i;
			nodes.add(node);
		}

		return nodes;
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

	/**
	 * Connects a directed graph with lattice links between X and X - 1 mod N.
	 * Each node has a single shortcut edge with an endpoint determined by the link length source.
	 * @param linkLengthSource Provides shortcut endpoints.
	 */
	public static Graph connectSandberg(ArrayList<SimpleNode> nodes, LinkLengthSource linkLengthSource) {
		Graph g = new Graph(nodes);

		// Base graph: Edge from X to X - 1 mod N for all nodes 0 to N - 1.
		for (int i = 0; i < nodes.size(); i++) {
			// Modulo of negative not defined: manually wrap.
			int wrapped = i - 1;
			if (wrapped < 0) wrapped += nodes.size();
			g.getNode(i).connectOutgoing(g.getNode(wrapped));
		}

		// Shortcuts: Edges from each node to an endpoint.
		for (int i = 0; i < nodes.size(); i++) {
			g.getNode(i).connectOutgoing(linkLengthSource.getPeer(g.getNode(i)));
		}

		return g;
	}

	/**
	 * Generates a graph with link length distribution and peer count distribution as described in the given sources.
	 * @param rand used for random numbers
	 * @return specified graph
	 */
	public static Graph connectGraph(ArrayList<SimpleNode> nodes, Random rand, LinkLengthSource linkLengthSource) {
		Graph g = new Graph(nodes);

		DistanceEntry[] distances = new DistanceEntry[nodes.size()];
		for (int i = 0; i < nodes.size(); i++) {
			SimpleNode src = g.nodes.get(i);
			if (src.atDegree()) continue;
			SimpleNode dest;

				// Fill distance entry array.
				for (int j = 0; j < nodes.size(); j++) {
					distances[j] = new DistanceEntry(Location.distance(src.getLocation(), g.nodes.get(j).getLocation()), j);
				}

				//System.out.println("distances size is " + );
				Arrays.sort(distances);

				// Make connections until at desired degree.
				while (!src.atDegree()) {
					dest = linkLengthSource.getPeer(src);
					if (src == dest || src.isConnected(dest) ||
					    (dest.atDegree() && rand.nextDouble() < rejectProbability)) continue;
					src.connect(dest);
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
			 * TODO: Assumes undirected connections.
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
			final Graph graph = new Graph(new ArrayList<SimpleNode>(networkSize));

			// Nodes.
			for (int i = 0; i < networkSize; i++) {
				SimpleNode node = (SimpleNode)input.readObject();
				node.setRand(random);
				node.index = i;
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

	/**
	 * Edge length distribution. Treats edges as directed.
	 *
	 * @param includeLatticeLinks If true, links from a node at index X to X - 1 mod N will be included. If false
	 *                            they will not.
	 */
	public ArrayList<Double> edgeLengths(final boolean includeLatticeLinks) {
		ArrayList<Double> lengths = new ArrayList<Double>();
		for (SimpleNode node : nodes) {
			for (SimpleNode peer : node.getConnections()) {
				if (!includeLatticeLinks && node.index == (peer.index + 1) % size()) continue;
				lengths.add(node.distanceToLoc(peer.getLocation()));
			}
		}
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

		//TODO: There are directed links now.
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

		int nWalks = 10 * 1000 * 1000;
		int nBuckets = 400;
		int hopsPerWalkUniform = 20;
		int hopsPerWalkCorrected = 40;

		int nTrials = 4;
		int[][][] pdfs = new int[nTrials][3][nBuckets];

		for (int trial = 0; trial < nTrials; trial++) {
			System.out.println("Creating test graph...");
			Random rand = new MersenneTwister(trial);
			final ArrayList<SimpleNode> nodes = Graph.generateNodes(nNodes, rand, true, new PoissonDegreeSource(12));
			Graph g = connectGraph(nodes, rand, new KleinbergLinkSource(rand, nodes));
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
