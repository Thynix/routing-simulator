package org.freenetproject.routing_simulator.graph;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Pair;
import org.freenetproject.routing_simulator.graph.degree.DegreeSource;
import org.freenetproject.routing_simulator.graph.degree.PoissonDegreeSource;
import org.freenetproject.routing_simulator.graph.linklength.KleinbergLinkSource;
import org.freenetproject.routing_simulator.graph.linklength.LinkLengthSource;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;
import org.freenetproject.routing_simulator.util.ArrayStats;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

/**
 * Class to represent and generate graphs of small world networks.
 * At present only limited Kleinberg graphs are generated.
 * Functions to evaluate the graph topology are also provided.
 */
public class Graph {
	private final ArrayList<SimpleNode> nodes;

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

	public static ArrayList<SimpleNode> generateNodes(final int nNodes, final RandomGenerator rand, boolean fastGeneration, DegreeSource source) {
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
			nodes.add(new SimpleNode(locations[i], rand, source.getDegree(), i));
		}

		return nodes;
	}

	/**
	 * Connects a directed graph with lattice links between X and X - 1 mod N.
	 * Each node has a single shortcut edge with an endpoint determined by the link length source.
	 * @param nodes Nodes which make up the network.
	 * @param linkLengthSource Provides shortcut endpoints.
	 * @param shortcuts Number of shortcut edges. Must be non-negative.
	 */
	public static Graph connectSandberg(ArrayList<SimpleNode> nodes, int shortcuts, LinkLengthSource linkLengthSource) {
		Graph g = new Graph(nodes);

		// Base graph of lattice edges: Edge from X to X - 1 mod N for all nodes 0 to N - 1.
		g.addLatticeLinks(true);

		// Shortcuts: Edges from each node to an endpoint.
		for (SimpleNode origin : g.nodes) {
			SimpleNode endpoint;
			// -1 to account for the single lattice edge.
			while (origin.degree() - 1 < shortcuts) {
				do {
					endpoint = linkLengthSource.getPeer(origin);
				} while (origin.isConnected(endpoint));
				origin.connectOutgoing(endpoint);
			}
		}

		return g;
	}

	/**
	 * Adds lattice links. Should be the first thing to add edges to a network.
	 *
	 * Connects node X to X - 1 mod N for all X, where X is node index and N is the network size.
	 *
	 * @param directed If true, the lattice links will be directed. If false, undirected.
	 */
	private void addLatticeLinks(final boolean directed) {
		assert nEdges() == 0;

		for (int i = 0; i < size(); i++) {
			/*
			 * From X to X - 1, wrapped to the network size. Java implements modulus such that it produces
			 * -1 for -1 % N, not N - 1 as it does in the definition of the lattice links. Going from
			 * X + 1 to X is equivalent.
			 */
			final SimpleNode from = getNode((i + 1) % size());
			final SimpleNode to = getNode(i);
			if (from.isConnected(to) || to.isConnected(from)) {
				throw new IllegalStateException("Connection already existed.");
			}

			if (directed) from.connectOutgoing(to);
			else from.connect(to);
		}

		// There is an edge between every node in a circle.
		assert nEdges() == size();
	}

	/**
	 * Gives the node in question random connections until it meets its desired degree. Does not add disconnected
	 * peers.
	 *
	 * @param node node to add connections to.
	 * @param random source of entropy for selecting which nodes to connect to.
	 */
	public List<SimpleNode> bootstrap(final SimpleNode node, final RandomGenerator random) {
		List<SimpleNode> disconnectedNodes = new ArrayList<SimpleNode>();
		SimpleNode peer;
		do {
			peer = getNode(random.nextInt(size()));

			/*
			 * Do not connect to self - reference comparison should be sufficient, or make a duplicate
			 * connection. Avoid connecting to disconnected nodes lest it fragment the network.
			 */
			if (node == peer || node.isConnected(peer) || peer.degree() == 0) continue;

			// Reference comparison should be sufficient.
			assert !node.equals(peer);

			/*
			 * If the peer is already at its degree, connect only if not rejected.
			 * Drop a random connection to keep the connection count invariant.
			 */
			if (!peer.atDegree() || peer.atDegree() && random.nextDouble() > rejectProbability) {
				SimpleNode disconnected = peer.swapConnections(node);
				if (disconnected.degree() == 0) disconnectedNodes.add(disconnected);
			}
		} while (!node.atDegree());

		return disconnectedNodes;
	}

	/**
	 * Connects a graph such that all nodes have a single (non-lattice, if possible) undirected connection to a
	 * single super node. Ignores nodes' desired degree.
	 *
	 * @param nodes Nodes which make up the network.
	 * @param lattice If true, adds lattice edges. If false, does not.
	 *
	 * @return Graph with the specified edges added.
	 */
	public static Graph connectSuperNode(ArrayList<SimpleNode> nodes, boolean lattice) {
		Graph graph = new Graph(nodes);
		assert nodes.size() > 1;
		if (lattice) graph.addLatticeLinks(false);

		final SimpleNode superNode = graph.getNode(0);
		for (int i = 1; i < graph.size(); i++) {
			SimpleNode peer = graph.getNode(i);
			if (!superNode.isConnected(peer)) superNode.connect(peer);
		}

		return graph;
	}

	/**
	 * Adds links to a graph which conform to the link length distribution and peer count distribution given.
	 *
	 * @param g Graph to add edges to.
	 * @param rand Provides probabilities for whether to connect to a node which is already at its desired degree.
	 * @param linkLengthSource Provides peers which give conforming connections.
	 *
	 * @return Graph with specified edges added.
	 */
	public static Graph connectGraph(Graph g, RandomGenerator rand, LinkLengthSource linkLengthSource) {
		SimpleNode destination;
		for (SimpleNode src : g.nodes) {
			if (src.atDegree()) continue;

			// Make connections until at desired degree.
			while (!src.atDegree()) {
				destination = linkLengthSource.getPeer(src);
				if (src == destination || src.isConnected(destination) ||
				    (destination.atDegree() && rand.nextDouble() < rejectProbability)) continue;
				src.connect(destination);
			}
		}

		return g;
	}

	public static Graph connectGraph(ArrayList<SimpleNode> nodes, RandomGenerator rand, LinkLengthSource linkLengthSource) {
		return connectGraph(new Graph(nodes), rand, linkLengthSource);
	}

	/**
	 * Generates a graph with undirected lattice connections, and link length distribution and peer count
	 * distribution as described in the given sources.
	 *
	 * @param nodes Nodes which make up the network.
	 *
	 * @see Graph#connectGraph(Graph, org.apache.commons.math3.random.RandomGenerator, org.freenetproject.routing_simulator.graph.linklength.LinkLengthSource)
	 */
	public static Graph connectGraphLattice(ArrayList<SimpleNode> nodes, RandomGenerator rand, LinkLengthSource linkLengthSource) {
		Graph graph = new Graph(nodes);
		graph.addLatticeLinks(false);
		return Graph.connectGraph(graph, rand, linkLengthSource);
	}

	/**
	 * Writes graph to a file. Format:
	 * <ul>
	 *      <li>Number of nodes.</li>
	 *      <li>SimpleNodes.</li>
	 *      <li>Connections: index from, index to</li>
	 * </ul>
	 * @param output stream to write graph to.
	 */
	public void write(DataOutputStream output) {
		try {

			// Number of nodes.
			output.writeInt(nodes.size());

			// Nodes.
			for (SimpleNode node : nodes) node.write(output);

			/*
			 * Write every connection; undirected edges are two directed edges.
			 */
			final ArrayList<Integer> connectionIndexes = new ArrayList<Integer>();
			int writtenConnections = 0;
			for (SimpleNode from : nodes) {
				for (SimpleNode to : from.getConnections()) {
					writtenConnections++;
					connectionIndexes.add(from.index);
					connectionIndexes.add(to.index);
				}
			}

			output.writeInt(writtenConnections);
			System.out.println("Writing " + writtenConnections + " connections.");
			assert connectionIndexes.size() == writtenConnections * 2;
			for (Integer index : connectionIndexes) {
				output.writeInt(index);
			}

			output.flush();
			output.close();
		} catch (IOException e) {
			System.err.println("Could not write to output stream:");
			e.printStackTrace();
			System.exit(3);
		}
	}

	/**
	 * Constructs the graph from a file.
	 * @param input stream to read the graph from.
	 * @param random Randomness source to give to nodes.
	 * @return graph defined by the file.
	 */
	public static Graph read(DataInputStream input, RandomGenerator random) {
		try {
			// Number of nodes.
			final int networkSize = input.readInt();
			final Graph graph = new Graph(new ArrayList<SimpleNode>(networkSize));

			// Nodes.
			for (int i = 0; i < networkSize; i++) {
				graph.nodes.add(new SimpleNode(input, i, random));
			}

			final int writtenConnections = input.readInt();
			System.out.println("Reading " + writtenConnections + " connections.");
			//Each connection consists of two indexes in a pair.
			for (int i = 0; i < writtenConnections; i++) {
				final int from = input.readInt();
				final int to = input.readInt();
				//System.out.println(from + " " + to);
				graph.nodes.get(from).connectOutgoing(graph.nodes.get(to));
			}

			return graph;
		} catch (IOException e) {
			System.err.println("Could not read from input stream:");
			e.printStackTrace();
			System.exit(4);
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
			System.out.print(ccStats.mean() + "\t" + ccStats.stdDev() + "\t" + ccStats.skewness() + "\t");
			System.out.print(degStats.mean() + "\t" + degStats.stdDev() + "\t" + degStats.skewness() + "\t");
		}
	}

	/**Print column headers for printGraphStats(false).*/
	public static void printGraphStatsHeader() {
		System.out.print("nNodes\tnEdges\tminDegree\tmaxDegree\tglobalClusterCoeff\tlocalCCMean\tlocalCCStdDev\tlocalCCSkew\t");
		System.out.print("degreeMean\tdegreeStdDev\tdegreeSkew\t");
	}

	/**Get the topology stats as an array.*/
	public double[] graphStats() {
		double[] cc = localClusterCoeff();
		int[] deg = degrees();
		ArrayStats ccStats = new ArrayStats(cc);
		ArrayStats degStats = new ArrayStats(deg);
		return new double[] {size(), nEdges(), minDegree(), maxDegree(), globalClusterCoeff(),
			ccStats.mean(), ccStats.stdDev(), ccStats.skewness(),
			degStats.mean(), degStats.stdDev(), degStats.skewness(),
		};
	}

	/**
	 * Edge length distribution. Treats edges as directed.
	 *
	 * @param includeLatticeLinks If true, links from nodes with adjacent indexes will be included. If false
	 *                            they will not.
	 */
	public ArrayList<Double> edgeLengths(final boolean includeLatticeLinks) {
		ArrayList<Double> lengths = new ArrayList<Double>();
		for (SimpleNode node : nodes) {
			for (SimpleNode peer : node.getConnections()) {
				if (!includeLatticeLinks) {
					if ( node.index == (peer.index + 1) % size() ||
					     peer.index == (node.index + 1) % size()) continue;
				}
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
		// Indexes in <lesser, greater> order of a connection.
		HashSet<Pair<Integer, Integer>> connections = new HashSet<Pair<Integer, Integer>>();
		int undirected = 0;

		for (SimpleNode origin : nodes) {
			for (SimpleNode peer : origin.getConnections()) {
				/*
				 * If the set already contained the element the connection is two mutual directed edges,
				 * which in this case is considered one directed edge.
				 */
				if (origin.index < peer.index) {
					if (!connections.add(new Pair<Integer, Integer>(origin.index, peer.index))) undirected++;
				} else {
					if (!connections.add(new Pair<Integer, Integer>(peer.index, origin.index))) undirected++;
				}
			}
		}

		System.out.println("Out of the edges " + undirected + " are undirected.");
		return connections.size();
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
		for (SimpleNode node : nodes) {
			int d = node.degree();
			sumDegrees += d;
			sumSquareDegrees += d * d;
		}

		return ((double)sumSquareDegrees)/((double)n) - ((double)(sumDegrees * sumDegrees))/((double)(n * n));
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
		for (SimpleNode node : nodes) {
			sumCoeff += node.localClusterCoeff();
		}
		double mean = sumCoeff / n;
		assert mean >= 0.0 && mean <= 1.0;
		return mean;
	}

	private double[] localClusterCoeff() {
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
	private double globalClusterCoeff() {
		int nClosed = 0;
		int nTotal = 0;

		for (SimpleNode n : nodes) {
			int degree = n.degree();
			nClosed += n.closedTriplets();
			nTotal += (degree * (degree - 1)) / 2;
		}

		return ((double)(nClosed)) / ((double)(nTotal));
	}

	private int[] randomWalkDistTest(int nWalks, int hopsPerWalk, boolean uniform, RandomGenerator rand) {
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
			RandomGenerator rand = new MersenneTwister(trial);
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
