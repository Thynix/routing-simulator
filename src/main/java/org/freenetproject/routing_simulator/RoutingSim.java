package org.freenetproject.routing_simulator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.freenetproject.routing_simulator.graph.Graph;
import org.freenetproject.routing_simulator.graph.degree.ConformingDegreeSource;
import org.freenetproject.routing_simulator.graph.degree.DegreeSource;
import org.freenetproject.routing_simulator.graph.degree.FixedDegreeSource;
import org.freenetproject.routing_simulator.graph.degree.PoissonDegreeSource;
import org.freenetproject.routing_simulator.graph.linklength.ConformingLinkSource;
import org.freenetproject.routing_simulator.graph.linklength.KleinbergLinkSource;
import org.freenetproject.routing_simulator.graph.linklength.LinkLengthSource;
import org.freenetproject.routing_simulator.graph.linklength.UniformLinkSource;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;
import org.freenetproject.routing_simulator.util.ArrayStats;
import org.freenetproject.routing_simulator.util.MersenneTwister;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Class to perform routing simulations on Graphs.
 */
public class RoutingSim {
	//generic output format
	static final DecimalFormat outputFormat = new DecimalFormat("0.000000");

	/**
	 * Checks that a path is a directory which can be written to, and attempts to create it if it does not exist.
	 * Outputs descriptive messages if given anything other than an existing writable directory.
	 * @param path path to check
	 * @return True if the directory (now) exists and is writable; false otherwise.
	 */
	private static boolean writableDirectory(String path) {
		File file = new File(path);
		if (!file.exists()) {
			if (!file.mkdirs()) {
				System.out.println("Unable to create degree output directory \"" + file.getAbsolutePath() + "\".");
				return false;
			} else {
				System.out.println("Degree output directory \"" + file.getAbsolutePath() + "\" did not exist, so it was created.");
			}
		} else if (!file.isDirectory()) {
			System.out.println("Degree output path \"" + file.getAbsolutePath() + "\" is not a directory as it should be.");
			return false;
		} else if (!file.canWrite()) {
			System.out.println("No write access to degree output directory \"" + file.getAbsolutePath() + "\".");
			return false;
		}
		return true;
	}

	/**
	 * Attempts to open the given path for writing.
	 * @param path File to open.
	 * @return If successful, an output stream for the file. If not, outputs a message and returns null.
	 */
	private static FileOutputStream writableFile(String path) {
		try {
			return new FileOutputStream(new File(path));
		} catch (FileNotFoundException e) {
			System.out.println("Unable to open \"" + path + "\" for output:");
			e.printStackTrace();
			return null;
		}
	}

	private static boolean readableFile(String path) {
		File file = new File(path);
		if (!file.exists() || !file.canRead() || !file.isFile()) {
			System.out.println("Cannot read \"" + file.getAbsolutePath() + "\" as a file.");
			return false;
		}
		return true;
	}

	/**
	 * Main simulator program.  Generate a set of graphs of different
	 * parameters, run a set of requests on them, and print assorted
	 * stats.
	 *
	 * @param args Command-line arguments; not used.
	 */
	public static void main(String[] args) throws ParseException {
		//TODO: All this options stuff is an ugly, verbose mess. Is there some way to clean it up or at least move it somewhere else? Maybe main only does options parsing? (Then calls other stuff.)
		Options options = new Options();
		//TODO: Default values for arguments.
		//TODO: Line lengths.
		//Overall
		options.addOption("D", "output-degree", true, "Output file for degree distribution.");
		options.addOption("L", "output-link", true, "Output file for link length distribution.");
		options.addOption("I", "include-lattice", false, "Include links from index X to X - 1 mod N when outputting link lengths. If the graph does not actually have lattice connections this is recommended.");
		options.addOption("q", "quiet", false, "No simulation output to stdout. Messages about arguments are still output.");
		options.addOption("v", "verbose", false, "Progress updates.");
		options.addOption("h", "help", false, "Display this message.");
		options.addOption("v", "version", false, "Display software version.");
		options.addOption("S", "seed", true, "Seed used by psuedorandom number generator.");

		//Graphs: General generation options
		//TODO: Scale degree distribution (Also results?) to arbitrary network size - attempt to avoid distortion.
		options.addOption("s", "size", true, "Number of nodes in the network. Currently ignored when using --degree unless --force-size is specified.");
		//TODO: Does it make sense to use fastGeneration without evenSpacing? Assuming it doesn't.
		options.addOption("f", "fast-generation", false, "If present, the simulator will assign locations with even spacing and, when using --ideal-link, take shortcuts to speed up graph generation.");
		options.addOption("G", "load-graph", true, "Path to load a saved graph from.");
		options.addOption("g", "save-graph", true, "Path to save a graph after simulation is run on it.");
		options.addOption("s", "sandberg-graph", false, "Generate a directed graph with an edge from x to x -1 mod N for all x = 0 ... N - 1 as described in early section 2.2.1 in \"Searching in a Small World\".");

		//Graphs: link length distribution
		options.addOption("l", "ideal-link", false, "Kleinberg's ideal distribution: proportional to 1/d.");
		options.addOption("f", "flat-link", false, "Intentionally terrible distribution: uniformly random.");
		options.addOption("c", "conforming-link", true, "Distribution conforming to a file. Takes a path to a degree distribution file of the format \"[degree] [number of occurrences]\\n\"\"");

		//Graphs: degree distribution
		options.addOption("F", "fixed-degree", true, "All nodes are as close to the specified degree as practical.");
		options.addOption("C", "conforming-degree", true, "Distribution conforming to a file. Takes a path to a degree distribution file of the format \"[degree] [number of occurrences]\\n\"");
		options.addOption("i", "poisson-degree", true, "Distribution conforming to a Poisson distribution with the given mean.");

		//Simulations: Routing policies
		options.addOption("R", "route", true, "Simulate routing the given number of requests. Requires that --output-route, --fold-policy, and --output-hops be specified.");;
		options.addOption("o", "output-route", true, "File to which routing information is output.");
		StringBuilder description = new StringBuilder("Path folding policy:");
		for (PathFolding policy : PathFolding.values()) description.append(" ").append(policy);
		options.addOption("P", "fold-policy", true,  description.toString());
		options.addOption("H", "output-hops", true, "Base filename to output hop histograms for each sink policy. Appended with -<policy-num> for each.");

		//Simulations: Probe distribution
		options.addOption("p", "probe", true, "Simulate running probes from random locations for the specified number of maximum hops. Requires that --output-probe be specified.");
		options.addOption("m", "metropolis-hastings", false, "If present, probes will be routed with Metropolis-Hastings correction. If not, peers will be selected entirely at random.");
		options.addOption("O", "output-probe", true, "Directory to which probe distribution is output as \"[node ID] [times seen]\\n\" for a reference of random selection from the whole and at each hop up to the specified maximum hops.");

		CommandLineParser parser = new GnuParser();
		CommandLine cmd = parser.parse(options, args);

		if (cmd.hasOption("help")) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "java -jar simulator.jar", options );
			return;
		}

		if (cmd.hasOption("version")) {
			System.out.println("Freenet Routing Simulator v 0.0.1-dev");
			return;
		}

		//Check that required arguments are specified and that combinations make sense.
		if (cmd.hasOption("quiet") && cmd.hasOption("verbose")) {
			System.out.println("Quiet with verbose does not make sense.");
			return;
		}
		if (cmd.hasOption("quiet") && !((cmd.hasOption("probe") && cmd.hasOption("output-probe")) || cmd.hasOption("output-degree") || cmd.hasOption("output-link") || cmd.hasOption("save-graph") || cmd.hasOption("route"))) {
			System.out.println("Simulation will produce no output: --quiet is specified, but not any option which outputs to a file.");
		}

		int degreeOptions = 0;
		for (String option : new String[] { "fixed-degree", "conforming-degree" }) {
			if (cmd.hasOption(option)) degreeOptions++;
		}

		int linkOptions = 0;
		for (String option : new String[] { "ideal-link", "flat-link", "conforming-link" }) {
			if (cmd.hasOption(option)) linkOptions++;
		}

		if (degreeOptions > 1 || linkOptions > 1) {
			System.out.println("Graph cannot be generated with multiple methods at once.");
			return;
		}
		if (degreeOptions == 0 && linkOptions == 0 && !cmd.hasOption("load-graph") && !cmd.hasOption("sandberg-graph")) {
			System.out.println("No graph generation method specified.");
			return;
		}

		if (cmd.hasOption("route") && !(cmd.hasOption("output-route") && cmd.hasOption("fold-policy") && cmd.hasOption("output-hops"))) {
			System.out.println("--route was specified, but not all of its required parameters: --output-route, --fold-policy, --output-hops.");
			return;
		}
		if (cmd.hasOption("probe") && !cmd.hasOption("output-probe")) {
			System.out.println("--probe was specified, but not --output-probe.");
			return;
		}

		if (!cmd.hasOption("size") && !cmd.hasOption("load-graph")) {
			System.out.println("Network size not specified. (--size)");
			return;
		}

		final PathFolding pathFolding;
		if (cmd.hasOption("fold-policy")) {
			try {
				pathFolding = PathFolding.valueOf(cmd.getOptionValue("fold-policy"));
			} catch (IllegalArgumentException e) {
				System.out.println("The folding policy \"" + cmd.getOptionValue("fold-policy") + "\" is invalid.");
				System.out.println("Possible values are:");
				for (PathFolding policy : PathFolding.values()) {
					System.out.println(policy.toString());
				}
				e.printStackTrace();
				System.exit(15);
				/*
				 * Data flow analysis does not appear to take into account that System.exit() ends the
				 * program, so the return is here to prevent pathFolding from being possibly
				 * uninitialized.
				 */
				return;
			}
		} else {
			pathFolding = PathFolding.NONE;
		}

		//Check for problems with specified paths.
		//Check if input files can be read.
		if (cmd.hasOption("conforming-degree") && ! readableFile(cmd.getOptionValue("conforming-degree"))) System.exit(9);
		if (cmd.hasOption("conforming-link") && ! readableFile(cmd.getOptionValue("conforming-link"))) System.exit(10);
		if (cmd.hasOption("load-graph") && !readableFile(cmd.getOptionValue("load-graph"))) System.exit(11);

		//Check if output paths are directories that can be written to, and create them if they do not exist.
		if (cmd.hasOption("output-probe") && !writableDirectory(cmd.getOptionValue("output-probe"))) System.exit(12);

		//Check that output files exist and are writable or can be created.
		FileOutputStream degreeOutput = null, linkOutput = null;
		final File graphOutput;
		if (cmd.hasOption("output-degree") && (degreeOutput = writableFile(cmd.getOptionValue("output-degree"))) == null ) System.exit(13);
		if (cmd.hasOption("output-link") && (linkOutput = writableFile(cmd.getOptionValue("output-link"))) == null) System.exit(14);
		if (cmd.hasOption("save-graph")) {
			graphOutput = new File(cmd.getOptionValue("save-graph"));
			//Just check for this one; saving takes a File. Better to check here than after simulation runs.
			try {
				FileOutputStream outputStream = new FileOutputStream(graphOutput);
				outputStream.close();
			} catch (FileNotFoundException e) {
				System.err.println("Could not open saved graph:");
				e.printStackTrace();
				System.exit(7);
			} catch (IOException e) {
				System.err.println("Unexpected IOException:");
				e.printStackTrace();
				System.exit(8);
			}
		} else {
			graphOutput = null;
		}

		//TODO: What is a sink policy? Looks like 2 is used as a hard-coded second dimension in Request, so this isn't currently a good configuration target.
		int[] sinkPolsUsed = {0, 1};
		PrintStream[] histogramOutput = new PrintStream[sinkPolsUsed.length];
		for (int i = 0; i < sinkPolsUsed.length; i++) {
			histogramOutput[i] = new PrintStream(writableFile(cmd.getOptionValue("output-hops")+ "-" + sinkPolsUsed[i]));
		}

		int nRequests = cmd.hasOption("route") ? Integer.valueOf(cmd.getOptionValue("route")) : 4000;
		int maxHops = cmd.hasOption("probe") ? Integer.valueOf(cmd.getOptionValue("probe")) : 50;

		Random rand;

		//TODO: What is this used for?
		//TODO: Hard-coded dimensions bad - clean this up.
		double[][][] avgStats = new double[1][13][1];

		final boolean quiet = cmd.hasOption("quiet");
		final boolean verbose = cmd.hasOption("verbose");
		final int seed = cmd.hasOption("seed") ? Integer.valueOf(cmd.getOptionValue("seed")) : 0;

		if (!verbose) {
			System.out.println();
			System.out.print("p\tq\tpLow\tpInst\tevenSpacing\tfastGeneration\t");
			Graph.printGraphStatsHeader();
			System.out.println();
		}

		rand = new MersenneTwister(seed);

		//Time tracking: report time taken for each graph setting if verbose; upon completion otherwise.
		long startTime = System.currentTimeMillis();
		long lastTime = startTime;

		// Load the graph; otherwise generate.
		Graph g;
		if (cmd.hasOption("load-graph")) {
			g = Graph.read(new File(cmd.getOptionValue("load-graph")), rand);
		} else {
			final int networkSize = Integer.valueOf(cmd.getOptionValue("size"));

			final DegreeSource degreeSource;
			if (cmd.hasOption("conforming-degree")) degreeSource = new ConformingDegreeSource(cmd.getOptionValue("conforming-degree"), rand);
			else if (cmd.hasOption("poisson-degree")) degreeSource = new PoissonDegreeSource(Integer.valueOf(cmd.getOptionValue("poisson-degree")));
			else if (cmd.hasOption("fixed-degree")) degreeSource = new FixedDegreeSource(Integer.valueOf(cmd.getOptionValue("fixed-degree")));
			else /* if (cmd.hasOption("sandberg-graph") */ degreeSource = new FixedDegreeSource(0);

			final ArrayList<SimpleNode> nodes = Graph.generateNodes(networkSize, rand, true, degreeSource);

			final LinkLengthSource linkLengthSource;
			if (cmd.hasOption("conforming-link")) linkLengthSource = new ConformingLinkSource(cmd.getOptionValue("conforming-link"), rand, nodes);
			else if (cmd.hasOption("ideal-link")) linkLengthSource = new KleinbergLinkSource(rand, nodes);
			else if (cmd.hasOption("flat-link")) linkLengthSource = new UniformLinkSource(rand, nodes);
			else throw new IllegalStateException("Link length distribution undefined.");

				if (cmd.hasOption("sandberg-graph")) {
				g = Graph.connectSandberg(nodes, linkLengthSource);
				} else {
				g = Graph.connectGraph(nodes, rand, linkLengthSource);
			}
		}

		if (!quiet) {
			g.printGraphStats(verbose);
			System.out.println("Generation took (ms): " + (System.currentTimeMillis() - lastTime));
			lastTime = System.currentTimeMillis();
		}

		double[] indivStats = g.graphStats();
		//TODO: 13 is defined because this doubles array has 13 different types of information. This should be a class.
		for (int i = 0; i < 13; i++) {
			avgStats[0][i][0] = indivStats[i];
		}

		if (cmd.hasOption("probe")) {
			rand = new MersenneTwister(seed);
			//Uniform probes if --metropolis-hastings is not specified.
			probeDistribution(g, rand, maxHops, quiet, verbose, cmd.getOptionValue("output-probe"), !cmd.hasOption("metropolis-hastings"));
		}

		if (cmd.hasOption("route")) {
			rand = new MersenneTwister(seed);
			simulate(g, rand, nRequests, cmd.getOptionValue("output-route"), pathFolding, histogramOutput);
		}

		if (cmd.hasOption("output-degree")) {
			int[] degrees = new int[g.maxDegree() + 1];
			for (int degree : g.degrees()) {
				degrees[degree]++;
			}
			for (int i = 0; i < g.maxDegree(); i++) {
				try {
					degreeOutput.write((i + " " + degrees[i] + "\n").getBytes());
				} catch (IOException e) {
					System.out.println("Unexpected error encoding string for degree output:");
					System.out.println(e);
					e.printStackTrace();
					return;
				}
			}
		}

		if (cmd.hasOption("output-link")) {
			ArrayList<Double> lengths = g.edgeLengths(cmd.hasOption("include-lattice"));
			//Output is intended for gnuplot CDF - second value is Y and should sum to 1.
			double normalized = 1.0/lengths.size();
			for (double length : lengths) {
				try {
					linkOutput.write((length + " " + normalized + "\n").getBytes());
				} catch (IOException e) {
					System.out.println("Unexpected error encoding string for link length output:");
					System.out.println(e);
					e.printStackTrace();
					return;
				}
			}
		}

		if (cmd.hasOption("save-graph")) g.write(graphOutput);

		if (verbose) System.out.println();
		if (!quiet) {
			System.out.println("Time taken (ms): " + (System.currentTimeMillis() - lastTime));
		}

		if (!quiet) {
			System.out.println("Average stats:");
			Graph.printGraphStatsHeader();
			//TODO: Why is this hardcoded to 13?
			for (int j = 0; j < 13; j++) {
				//TODO: Hard-coded dimensions bad.
				ArrayStats s = new ArrayStats(avgStats[0][j]);
				System.out.print("(" + outputFormat.format(s.mean()) + " " + outputFormat.format(s.stdDev()) + ")\t");
			}
			System.out.println();
			System.out.println("Total time taken (ms): " + (System.currentTimeMillis() - startTime));
		}
	}

	/**
	 * @param array to convert
	 * @return string in which "[index] [value]" pairs are newline-delimited.
	 */
	private static String stringArray(int[] array) {
		String s = "";
		for (int i = 0; i < array.length; i++) s += i + " " + array[i] + "\n";
		return s;
	}

	private static void writeArray(int[] array, File target) {
		try {
			FileOutputStream outputStream = new FileOutputStream(target);
			outputStream.write(stringArray(array).getBytes());
		} catch (FileNotFoundException e) {
			System.out.println("Cannot open file \"" + target.getAbsolutePath() + "\" for writing.");
			System.out.println(e);
			System.exit(1);
		} catch (IOException e) {
			System.out.println(e);
			System.exit(1);
		}
	}

	public static void probeDistribution(Graph g, Random rand, int maxHops, boolean quiet, boolean verbose, final String containingPath, boolean uniform) {
		File output = new File(containingPath);
		assert output.isDirectory();
		if (!output.exists()) {
			if (!output.mkdirs()) {
				System.out.println("Unable to create requested output directory \"" + containingPath + "\".");
				System.exit(1);
			}
		}

		//TODO: nTrials and nProbes configurable on command line.
		final int nTrials = 100;
		final int nProbes = g.size() * 30;
		System.out.println("Determining baseline");
		int[] baselineOccurrences = new int[g.size()];
		/*
		 * Find baseline for visibility by selecting the same number of nodes from the entire network at random
		 * as endpoints at each HTL. Sort occurrences each run, then add to final result array to represent
		 * actual spread from each run and avoid node index influence.
		 */
		for (int i = 0; i < nTrials; i++) {
			int[] trialOccurrences = new int[g.size()];
			for (int walk = 0; walk < nProbes; walk++) {
				trialOccurrences[g.getNode(rand.nextInt(g.size())).index]++;
			}
			Arrays.sort(trialOccurrences);
			assert baselineOccurrences.length == trialOccurrences.length;
			for (int j = 0; j < trialOccurrences.length; j++) {
				baselineOccurrences[j] += trialOccurrences[j];
			}
		}

		output = new File(containingPath + "reference.dat");
		writeArray(baselineOccurrences, output);

		System.out.println("Simulating HTL");
		//Find distribution of nodes reached with random walk for increasing hops from all nodes.
		//maxHops + 1 is because the starting node is at zero hops.
		int[][] hopOccurrences = new int[maxHops + 1][g.size()];
		ArrayList<SimpleNode> trace;
		for (int nodeIndex = 0; nodeIndex < nTrials; nodeIndex++) {
			SimpleNode source = g.getNode(rand.nextInt(g.size()));
			SimpleNode alongTrace;
			int[][] trialOccurrences = new int[maxHops + 1][g.size()];
			for (int walk = 0; walk < nProbes; walk++) {
				trace = source.randomWalkList(maxHops, uniform, rand);
				//Traces: starting point (zero hops), then maxHops hops from there.
				assert trace.size() == maxHops + 1;
				for (int fromEnd = 0; fromEnd <= maxHops; fromEnd++) {
					//fromEnd of trace: hops along. 0 is starting node.
					alongTrace = trace.get(fromEnd);
					trialOccurrences[fromEnd][alongTrace.index]++;
				}
			}
			assert hopOccurrences.length == trialOccurrences.length;
			for (int i = 0; i < trialOccurrences.length; i++) {
				assert hopOccurrences[i].length == trialOccurrences[i].length;
				Arrays.sort(trialOccurrences[i]);
				for (int j = 0; j < trialOccurrences[i].length; j++) {
					hopOccurrences[i][j] += trialOccurrences[i][j];
				}
			}
		}

		System.out.println("Sorting results.");
		for (int hops = 0; hops <= maxHops; hops++) {
			output = new File(containingPath + "probe-" + hops + ".dat");
			writeArray(hopOccurrences[hops], output);
		}
	}

	public static void simulate(Graph g, Random rand, int nRequests, final String outputPath,
	                            final PathFolding policy, final PrintStream[] histogramOutput) {
		File outputFile = new File(outputPath);
		PrintStream stream = null;
		try {
			stream = new PrintStream(outputFile);
		} catch (IOException e) {
			System.err.println("Cannot open file \"" + outputPath + ": " + e);
			System.exit(1);
		}
		stream.println("Routing " + nRequests + " requests on network of size " + g.size() + ".");
		long startTime = System.currentTimeMillis();


		for (int i = 0; i < nRequests; i++) {
			final SimpleNode origin = g.getNode(rand.nextInt(g.size()));
			origin.greedyRoute(rand.nextDouble(), 50, policy);
		}
		//TODO: This is only used in 2D capacity in one line: can make 1D?
		//Request[][] requests = new Request[nRequests][nIntersectTests];

		//Separated into two loops so that the same set
		//of requests are generated for each routing
		//policy.
		//This is why the instances of Random must be the same each run: comparability between routing schemes.
		//TODO: Intersect means what?
		//int[][] requestOrigins = new int[nRequests][nIntersectTests];
		/*for (int i = 0; i < nRequests; i++) {
			double l = rand.nextDouble();
			//create multiple requests with same target location but different origin
			for (int j = 0; j < nIntersectTests; j++) {
				Request r = new Request(l, rand, routePolicy, false, g);
				requests[i][j] = r;
				requestOrigins[i][j] = rand.nextInt(g.size());
			}


		//Route all requests.
		//TODO: What does route do?
		for (int i = 0; i < nRequests; i++) {
			for (int j = 0; j < nIntersectTests; j++) {
				Request r = requests[i][j];
				SimpleNode origin = g.getNode(requestOrigins[i][j]);
				origin.route(r, null, policy);
			}
		}

		//Count how many requests were routed to their exact destination.
		int nPreciseRouted = 0;
		for (int i = 0; i < nRequests; i++) {
			for (int j = 0; j < nIntersectTests; j++) {
				if (g.preciseRoute(requests[i][j])) nPreciseRouted++;
			}
		}
		stream.println("Requests precisely routed: " + nPreciseRouted);
		//assert nPreciseRouted == nRequests * nIntersectTests;

		//TODO: What is a sink policy? Is "pol" policy?
		//For each sink policy, compute and print stats.
		//Sink policy has no effect on routing decisions,
		//so we compute all sink policies during one
		//routing trial.
		//TODO: Is there ever more than one sink policy?
		for (int sp = 0; sp < sinkPolsUsed.length; sp++) {
			int sinkPolicy = sinkPolsUsed[sp];
			//TODO: "Hist" is history? Histogram?
			int sinkHistSize = 13;
			int[] sinkHist = new int[sinkHistSize];
			int[] lowUptimeSinkHist = new int[sinkHistSize];
			int[] highUptimeSinkHist = new int[sinkHistSize];
			ArrayList<Integer> hopsHist = new ArrayList<Integer>(50);

			double meanSinkCount = 0.0;
			double meanLowUptimeSinkCount = 0.0;
			double meanHighUptimeSinkCount = 0.0;

			for (int i = 0; i < nRequests; i++) {
				for (int j = 0; j < nIntersectTests; j++) {
					Request r = requests[i][j];

					int nLowUptimeSinks = r.nLowUptimeSinks(sinkPolicy);
					int nHighUptimeSinks = r.nHighUptimeSinks(sinkPolicy);
					int nSinks = nLowUptimeSinks + nHighUptimeSinks;

					meanSinkCount += (double)nSinks;
					meanLowUptimeSinkCount += (double)nLowUptimeSinks;
					meanHighUptimeSinkCount += (double)nHighUptimeSinks;

					nSinks = Math.min(sinkHistSize - 1, nSinks);
					nLowUptimeSinks = Math.min(sinkHistSize - 1, nLowUptimeSinks);
					nHighUptimeSinks = Math.min(sinkHistSize - 1, nHighUptimeSinks);

					sinkHist[nSinks]++;
					lowUptimeSinkHist[nLowUptimeSinks]++;
					highUptimeSinkHist[nHighUptimeSinks]++;

					int nHops = r.hopsTaken();
					// Expand histogram to make room for current value instead of clamping.
					hopsHist.ensureCapacity(nHops + 1);
					while (hopsHist.size() < nHops + 1) hopsHist.add(0);
					hopsHist.set(nHops, hopsHist.get(nHops) + 1);
				}
			}

			stream.println("Hops histogram for policy " + sinkPolicy + ":");

			for (int i = 0; i < hopsHist.size(); i++) {
				histogramOutput[sp].println(i + " " + hopsHist.get(i));
			}
			histogramOutput[sp].println();
			histogramOutput[sp].flush();

			//Before this means are just totals; divide by number of iterations above.
			meanSinkCount /= nRequests * nIntersectTests;
			meanLowUptimeSinkCount /= nRequests * nIntersectTests;
			meanHighUptimeSinkCount /= nRequests * nIntersectTests;

			//TODO: Would these be more readable as actual histogram plots?
			stream.println("Sink count histograms for policy " + sinkPolicy + ":");
			stream.println("n\tTotal\tlowU\thighU");
			for (int i = 0; i < sinkHistSize; i++) {
				stream.println(i + "\t" + sinkHist[i] + "\t" + lowUptimeSinkHist[i] + "\t" + highUptimeSinkHist[i]);
			}
			stream.println("Mean sinks: " + meanSinkCount);
			stream.println("Mean low uptime sinks: " + meanLowUptimeSinkCount);
			stream.println("Mean high uptime sinks: " + meanHighUptimeSinkCount);
			stream.println();
		}

		//TODO: What does decrement mean? TODO: This is 2D (requests) to 1D (decrements)... 1D requests would be cleaner here.
		//Populate decrements array with htlDecrements() from each request.
		int[] decrements = new int[nRequests * nIntersectTests];
		for (int i = 0; i < nRequests; i++) {
			for (int j = 0; j < nIntersectTests; j++) {
				int d = requests[i][j].htlDecrements();
				decrements[i*nIntersectTests + j] = d;
			}
		}
		//TODO: Why not have the printArraySummary() sort it?
		Arrays.sort(decrements);
		stream.println("HTL Decrements:");
		stream.print(printArraySummary(decrements, true));

		//TODO: What is 19 for? The length is used for comparison with number of hops?
		double[] logDistance = new double[19];
		int[] logDistCount = new int[19];
		boolean byDecrements = false;

		//TODO: What is this doing? Filling out a histogram?
		for (int i = 0; i < nRequests; i++) {
			for (int j = 0; j < nIntersectTests; j++) {
				Request r = requests[i][j];
				int[] hops;
				if (byDecrements) {
					hops = r.getPathDecrements();
				} else {
					hops = r.getPathHTL();
				}
				double[] pathDistance = r.getPathDistance();
				for (int k = 0; k < hops.length; k++) {
					int dec = hops[k];
					if (dec >= logDistance.length) break;
					logDistCount[dec]++;
					double logDist = Math.log(pathDistance[k]) / Math.log(2);
					logDistance[dec] += logDist;
				}
			}
		}

		stream.println("Routing distance by " + (byDecrements ? "decrements:" : "HTL:"));
		stream.println((byDecrements ? "Dec" : "HTL") + "\tCount\tLogDist");
		for (int i = 0; i < logDistance.length; i++) {
			double meanLogDist = logDistance[i];
			if (logDistCount[i] > 0) {
				meanLogDist /= logDistCount[i];
			} else {
				assert meanLogDist == 0.0;
			}
			stream.println(i + "\t" + logDistCount[i] + "\t" + meanLogDist);
		}
		stream.println();

		*/
		/*
		int[] maxHopsToIntersect = null;
		int[] pairedMaxHTI;
		int[][] hopsToSink;
		if (nIntersectTests > 1) {
			maxHopsToIntersect = new int[nRequests];
			pairedMaxHTI = new int[nRequests];
			hopsToSink = new int[sinkPolsUsed.length][nRequests];
			for (int i = 0; i < nRequests; i++) {
				maxHopsToIntersect[i] = -1;
				pairedMaxHTI[i] = -1;
				//for (int j = 0; j < nIntersectTests; j++) {
				//TODO: What?! This only runs once!
				for (int j = 0; j < 1; j++) {
					for (int k = 0; k < nIntersectTests; k++) {
						//TODO: So... k starts at 1 then. No self-routing?
						if (j == k) continue;
						int hti = requests[i][j].hopsToIntersect(requests[i][k]);
						maxHopsToIntersect[i] = Math.max(maxHopsToIntersect[i], hti);
						//TODO: Is zero selected to just have a sampling? What is the meaning of being at zero?
						if (j == 0 || k == 0) pairedMaxHTI[i] = Math.max(pairedMaxHTI[i], hti);
						if (j == 0) {
							for (int s = 0; s < sinkPolsUsed.length; s++) {
								hopsToSink[s][i] = Math.max(hopsToSink[s][i], requests[i][k].hopsToIntersect(requests[i][j], true, sinkPolsUsed[s]));
							}
						}
					}
				}*/
				/*Only valid if doing precise routing
				assert maxHopsToIntersect[i] >= 0;
				assert pairedMaxHTI[i] >= 0;
				*/
		/*
			}
			Arrays.sort(maxHopsToIntersect);
			stream.println("Max hops to intersection:");
			stream.print(printArraySummary(maxHopsToIntersect, true));
			if (printPairedMaxHTI) {
				Arrays.sort(pairedMaxHTI);
				stream.println("Paired max hops to intersection:");
				stream.print(printArraySummary(pairedMaxHTI, true));
			}
			stream.println();
			for (int s = 0; s < sinkPolsUsed.length; s++) {
				stream.println("Max hops to sink, policy " + sinkPolsUsed[s] + ":");
				Arrays.sort(hopsToSink[s]);
				stream.print(printArraySummary(hopsToSink[s], true));
				stream.println();
			}
		}*/
		stream.println("Time taken (ms): " + (System.currentTimeMillis() - startTime));
		stream.print("Summary:\t" + g.size() + "\t" + g.nEdges() + "\t" + g.minDegree() + "\t");
		stream.print(g.maxDegree() + "\t" + Math.sqrt(g.degreeVariance()) + "\t" + g.meanLocalClusterCoeff() + "\t");
		/*stream.print(g.globalClusterCoeff() + "\t" + nRequests * nIntersectTests + "\t" + nPreciseRouted + "\t");
		stream.print(printArraySummary(decrements, false));
		if (nIntersectTests > 1) {
			System.out.print(printArraySummary(maxHopsToIntersect, false));
		}
		stream.println();
		stream.println();*/
		stream.println();
	}

	/**
	 * Utility method to summarize a distribution.  Assumes array is
	 * sorted, as per Arrays.sort.
	 *
	 * @param a The array to summarize
	 * @param verbose Whether to print the long or short version
	 * @return A String summarizing the array
	 */
	public static String printArraySummary(int[] a, boolean verbose) {
		if (!isSorted(a)) throw new IllegalArgumentException("Array must be sorted.");
		double m = mean(a);
		double stdDev = stdDev(a);
		int pct50 = a[((int)(a.length * 0.5))];
		int pct90 = a[((int)(a.length * 0.9))];
		int pct97 = a[((int)(a.length * 0.97))];
		int pct99 = a[((int)(a.length * 0.99))];
		final StringBuilder summary = new StringBuilder();
		if (verbose) {
			summary.append("Mean:").append(m).append("\n");
			summary.append("Std Dev:").append(stdDev).append("\n");
			summary.append("50th percentile:").append(pct50).append("\n");
			summary.append("50th percentile:").append(pct50).append("\n");
			summary.append("90th percentile:").append(pct90).append("\n");
			summary.append("97th percentile:").append(pct97).append("\n");
			summary.append("99th percentile:").append(pct99).append("\n");
		} else {
			summary.append(m).append("\t").append(stdDev).append("\t").append(pct50).append("\t").append(pct90).append("\t").append(pct97).append("\t").append(pct99).append("\t");
		}
		return summary.toString();
	}

	public static double mean(int[] a) {
		double m = 0.0;
		for (int i = 0; i < a.length; i++) m += a[i];
		m /= a.length;
		return m;
	}

	public static double sumSquares(int[] a) {
		double ss = 0.0;
		for (int i = 0; i < a.length; i++) ss += ((double)a[i])*((double)a[i]);
		return ss;
	}

	public static double variance(int[] a) {
		double ss = sumSquares(a);
		double m = mean(a);
		return (ss / a.length) - m * m;
	}

	public static double stdDev(int[] a) {
		return Math.sqrt(variance(a));
	}

	/**
	 * @param a
	 * @return Whether the array is sorted in non-decreasing order.
	 */
	public static boolean isSorted(int[] a) {
		for (int i = 0; i < a.length - 1; i++) if (a[i] > a[i+1]) return false;
		return true;
	}

	/* TODO: Is there an alternative to copy-pasting to get array summary for both doubles and ints?
	 * Using Generics with <? extends Number> doesn't seem helpful as it wouldn't support operators and there are no
	 * common methods in Number.
	 */

	/**
	 * Utility method to summarize a distribution.  Assumes array is
	 * sorted, as per Arrays.sort.
	 *
	 * @param a The array to summarize
	 * @param verbose Whether to print the long or short version
	 * @return A String summarizing the array
	 */
	public static String printArraySummary(double[] a, boolean verbose) {
		if (!isSorted(a)) throw new IllegalArgumentException("Array must be sorted.");
		double m = mean(a);
		double stdDev = stdDev(a);
		double pct50 = a[((int)(a.length * 0.5))];
		double pct80 = a[((int)(a.length * 0.8))];
		double pct90 = a[((int)(a.length * 0.9))];
		double pct97 = a[((int)(a.length * 0.97))];
		double pct99 = a[((int)(a.length * 0.99))];
		String s;
		if (verbose) {
			s =		"Mean:			" + m + "\n";
			s = s +		"Std Dev:		" + stdDev + "\n";
			s = s +		"50th percentile:	" + pct50 + "\n";
			s = s +		"90th percentile:	" + pct90 + "\n";
			s = s +		"97th percentile:	" + pct97 + "\n";
			s = s +		"99th percentile:	" + pct99 + "\n";
		} else {
			s = m + "\t" + stdDev + "\t" + pct50 + "\t" + pct90 + "\t" + pct97 + "\t" + pct99 + "\t";
		}
		return s;
	}

	public static double mean(double[] a) {
		double m = 0.0;
		for (int i = 0; i < a.length; i++) m += a[i];
		m /= a.length;
		return m;
	}

	public static double sumSquares(double[] a) {
		double ss = 0.0;
		for (int i = 0; i < a.length; i++) ss += a[i] * a[i];
		return ss;
	}

	public static double variance(double[] a) {
		double ss = sumSquares(a);
		double m = mean(a);
		return (ss / a.length) - m * m;
	}

	public static double stdDev(double[] a) {
		return Math.sqrt(variance(a));
	}

	public static boolean isSorted(double[] a) {
		for (int i = 0; i < a.length - 1; i++) if (a[i] > a[i+1]) return false;
		return true;
	}
}
