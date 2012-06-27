package org.freenetproject.routing_simulator;

import org.apache.commons.cli.*;

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
	//E-series preferred numbers
	//TODO: What are these for? Don't appear to be used.
	public static final int[] Esix =    {	10,	15,	22,	33,	47,	68	};
	public static final int[] Etwelve = {	10, 12, 15, 18, 22, 27, 33, 39, 47, 56, 68, 82	};

	//generic output format
	static DecimalFormat outputFormat = new DecimalFormat("0.000000");

	private enum GraphGenerator {
		DEGREE,
		IDEAL
	}

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
			System.out.println(e);
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
		//TODO: Output routing simulation results to file
		options.addOption("q", "quiet", false, "No simulation output to stdout. Messages about arguments are still output.");
		options.addOption("v", "verbose", false, "Progress updates.");
		options.addOption("h", "help", false, "Display this message.");
		options.addOption("v", "version", false, "Display software version.");
		options.addOption("S", "seed", true, "Seed used by psuedorandom number generator.");

		//Graphs: General generation options
		//TODO: Scale degree distribution (Also results?) to arbitrary network size - attempt to avoid distortion.
		options.addOption("s", "size", true, "Number of nodes in the network. Currently ignored when using --degree unless --force-size is specified.");
		//TODO: Does it make sense to use fastGeneration without evenSpacing? Assuming it doesn't.
		options.addOption("f", "fast-generation", false, "If present, the simulator will assign locations as per --evenspacing and take shortcuts to speed up graph generation.");
		options.addOption("e", "even-spacing", false, "If present, the simulator will space nodes evenly throughout locations.");

		//Graphs: 1D Kleinberg
		options.addOption("i", "ideal", false, "Use an ideal 1D Kleinberg graph. Requires that --size, --local, --remote, --instant-reject, and --low-uptime be specified.");
		options.addOption("l", "local", true, "Number of local connections per node.");
		options.addOption("r", "remote", true, "Number of remote connections per node.");
		options.addOption("I", "instant-reject", true, "Probability between 0.0 and 1.0 that a connection is instantly rejected.");
		options.addOption("u", "low-uptime", true, "Probability between 0.0 and 1.0 that a node has low uptime.");

		//Graphs: From degree distribution
		options.addOption("d", "degree", true, "Use a graph following a given degree distribution. Takes a path to a degree distribution file of the format \"[degree] [number of occurrences]\\n\"");
		options.addOption("F", "force-size", false, "When using --degree force generation of --size nodes. May cause severe distortion.");
		options.addOption(null, "link", true, "Use a graph following a given link length distribution. Takes a path to a link distribution file containing a file for the format \"[link length] [arbitrary]\n");

		//Simulations: Routing policies
		//TODO: But what do the various numbers actually mean?
		options.addOption("R", "route", true, "Simulate routing policy of the specified number; possible policies are 1 through 6. Requires that --instant-reject, --low-uptime, --requests, --output-route, and --intersect-tests be specified.");
		//TODO: Explain more on these - what are their effects?
		options.addOption("q", "requests", true, "Number of requests to run.");
		options.addOption("n", "intersect-tests", true, "Number of intersect tests per request: same target but random origin.");
		options.addOption("o", "output-route", true, "File to which routing information is output.");

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
		if (cmd.hasOption("quiet") && !((cmd.hasOption("probe") && cmd.hasOption("output-probe")) || cmd.hasOption("output-degree") || cmd.hasOption("output-link"))) {
			System.out.println("Simulation will produce no output: --quiet is specified, but not any option which outputs to a file.");
		}
		if (cmd.hasOption("ideal") && cmd.hasOption("degree")) {
			System.out.println("Graph cannot be generated with multiple methods at once.");
			return;
		}
		if (!cmd.hasOption("ideal") && !cmd.hasOption("degree")) {
			System.out.println("No graph generation method specified.");
			return;
		}
		if (cmd.hasOption("ideal") && (!cmd.hasOption("size") || !cmd.hasOption("local") || !cmd.hasOption("remote"))) {
			System.out.println("--ideal was specified, but not one or more of its required parameters: --size, --local, --remote.");
			return;
		}
		if (cmd.hasOption("degree") && cmd.hasOption("force-size") && !cmd.hasOption("size")) {
			System.out.println("--degree and --force-size were specified but not --size.");
			return;
		}
		if (cmd.hasOption("route") && (!cmd.hasOption("requests") || !cmd.hasOption("intersect-tests") || !cmd.hasOption("instant-reject") || !cmd.hasOption("low-uptime") || !cmd.hasOption("output-route"))) {
			System.out.println("--route was specified, but not one or more of its required parameters: --requests, --intersect-tests, --instant-reject, --low-uptime.");
			return;
		}
		if (cmd.hasOption("probe") && !cmd.hasOption("output-probe")) {
			System.out.println("--probe was specified, but not --output-probe.");
			return;
		}

		//Check for problems with specified paths.
		//Check if input files can be read.
		if (cmd.hasOption("degree") && ! readableFile(cmd.getOptionValue("degree"))) return;
		if (cmd.hasOption("link") && ! readableFile(cmd.getOptionValue("link"))) return;

		//Check if output paths are directories that can be written to, and create them if they do not exist.
		if (cmd.hasOption("output-probe") && !writableDirectory(cmd.getOptionValue("output-probe"))) return;

		//Check that output files exist and are writable or can be created.
		FileOutputStream degreeOutput = null, linkOutput = null;
		if (cmd.hasOption("output-degree") && (degreeOutput = writableFile(cmd.getOptionValue("output-degree"))) == null ) return;
		if (cmd.hasOption("output-link") && (linkOutput = writableFile(cmd.getOptionValue("output-link"))) == null) return;


		final GraphGenerator graphType;
		if (cmd.hasOption("ideal")) graphType = GraphGenerator.IDEAL;
		else /*if (cmd.hasOption(""))*/ graphType = GraphGenerator.DEGREE;

		int nRequests = cmd.hasOption("requests") ? Integer.valueOf(cmd.getOptionValue("requests")) : 4000;
		int nIntersectTests = cmd.hasOption("intersect-tests") ? Integer.valueOf(cmd.getOptionValue("intersect-tests")) : 2;
		//TODO: What is a sink policy? Looks like 2 is used as a hard-coded second dimension in Request, so this isn't currently a good configuration target.
		int[] sinkPolsUsed = {0, 1};
		int routePol = cmd.hasOption("route") ? Integer.valueOf(cmd.getOptionValue("route")) : 3;
		int maxHops = cmd.hasOption("probe") ? Integer.valueOf(cmd.getOptionValue("probe")) : 50;

		Random rand;

		//TODO: What is this used for?
		//TODO: Hard-coded dimensions bad - clean this up.
		double[][][] avgStats = new double[1][13][1];

		//Time tracking: report time taken for each graph setting if verbose; upon completion otherwise.
		long startTime = System.currentTimeMillis();
		long lastTime = startTime;
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
				GraphParam gp = new GraphParam(
					cmd.hasOption("size") ? Integer.valueOf(cmd.getOptionValue("size")) : 2,
					cmd.hasOption("local") ? Integer.valueOf(cmd.getOptionValue("local")) : 0,
					cmd.hasOption("remote") ? Integer.valueOf(cmd.getOptionValue("remote")) : 0,
					cmd.hasOption("low-uptime") ? Double.valueOf(cmd.getOptionValue("low-uptime")) : 0,
					cmd.hasOption("instant-reject") ? Double.valueOf(cmd.getOptionValue("instant-reject")) : 0,
					cmd.hasOption("evenspacing") || cmd.hasOption("fast-generation"),
					cmd.hasOption("fast-generation"));

				if (verbose) {
					System.out.print("Generating " + graphType.name() + " graph of " + gp.n + " nodes, with ");
					System.out.println("parameters p = " + gp.p + ", q = " + gp.q + ".");
					System.out.print("pLowUptime = " + gp.pLowUptime + ", pInstantReject = " + gp.pInstantReject);
					System.out.println(", evenSpacing = " + gp.evenSpacing + ", fastGeneration = " + gp.fastGeneration);
				}

				Graph g;
				if (graphType == GraphGenerator.IDEAL) g = Graph.generate1dKleinbergGraph(gp, rand);
				else /*if (graphType == GraphGenerator.DEGREE)*/ {
					Graph.LinkLengthSource source;
					if (cmd.hasOption("link")) source = new Graph.ConformingLinkSource(cmd.getOptionValue("link"));
					else source = new Graph.UniformLinkSource();
					g = Graph.generatePeerDistGraph(gp, rand, cmd.getOptionValue("degree"), cmd.hasOption("force-size"), source);

				if (!quiet) g.printGraphStats(verbose);
				if (verbose) {
					System.out.println("Time taken (ms): " + (System.currentTimeMillis() - lastTime));
					lastTime = System.currentTimeMillis();
				}

				double[] indivStats = g.graphStats();
				//TODO: 13 is defined because this doubles array has 13 different types of information. This should be a class.
				for (int i = 0; i < 13; i++) {
					avgStats[0][i][0] = indivStats[i];
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
					double[] lengths = g.edgeLengths();
					//Output is intended for gnuplot CDF - second value is Y and should sum to 1.
					double normalized = 1.0/lengths.length;
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

				if (cmd.hasOption("probe")) {
					rand = new MersenneTwister(seed);
					//Uniform probes if --metropolis-hastings is not specified.
					probeDistribution(g, rand, maxHops, quiet, verbose, cmd.getOptionValue("output-probe"), !cmd.hasOption("metropolis-hastings"));
				}

				if (cmd.hasOption("route")) {
					rand = new MersenneTwister(seed);
					simulate(g, rand, nRequests, nIntersectTests, routePol, sinkPolsUsed, verbose, cmd.getOptionValue("output-route"));
				}
				if (verbose) System.out.println();
			if (!quiet) {
				System.out.println("Time taken (ms): " + (System.currentTimeMillis() - lastTime));
			}

		if (!quiet) {
			System.out.println("Average stats:");
			System.out.print("p\tq\tpLow\tpInst\tevenSpacing\tfastGeneration\t");
			Graph.printGraphStatsHeader();
			System.out.println();
				System.out.print(gp.p + "\t" + gp.q + "\t" + gp.pLowUptime + "\t" +
						gp.pInstantReject + "\t" + gp.evenSpacing + "\t" + gp.fastGeneration + "\t");
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

		/* 13:05:40    evanbd | As a general rule, for well-behaved (meaning: basically normal) data, your percentiles will be valid up until you have
		                      | around 30 points not included. So if you want a valid 99th percentile, you need to take 3000 data points.
		   13:05:46    evanbd | Plus or minus a bunch.
		   13:06:24    evanbd | This is a heuristic, not a rule; if you want a rule, you need to go put a confidence interval on your percentile, and that
		                      | gets all complicated.
		   13:07:02    evanbd | Note that the above number doesn't change much with network size.
		   13:07:13    evanbd | 30*(1/1-0.99)
		 */
		//TODO: 0.99 is desired percentile - any concerns for validity of occurrence distribution?
		final int nTrials =(int)(30*(1/(1-0.99)));
		final int nProbes = g.size();
		System.out.println("Determining baseline");
		int[] occurrences = new int[g.size()];
		/*
		 * Find baseline for visibility by selecting the same number of nodes from the entire network at random
		 * as endpoints at each HTL.
		 */
		for (int i = 0; i < nTrials; i++) {
			for (int walk = 0; walk < nProbes; walk++) {
				occurrences[g.getNode(rand.nextInt(g.size())).index]++;
			}
		}

		output = new File(containingPath + "reference.dat");
		writeArray(occurrences, output);

		System.out.println("Simulating HTL");
		//Find distribution of nodes reached with random walk for increasing hops from all nodes.
		//maxHops + 1 is because the starting node is at zero hops.
		int[][] hopOccurrences = new int[maxHops + 1][g.size()];
		ArrayList<SimpleNode> trace;
		for (int nodeIndex = 0; nodeIndex < nTrials; nodeIndex++) {
			SimpleNode source = g.getNode(rand.nextInt(g.size()));
			SimpleNode alongTrace;
			for (int walk = 0; walk < nProbes; walk++) {
				trace = source.randomWalkList(maxHops, uniform, rand);
				//Traces: starting point (zero hops), then maxHops hops from there.
				assert trace.size() == maxHops + 1;
				for (int fromEnd = 0; fromEnd <= maxHops; fromEnd++) {
					//fromEnd of trace: hops along. 0 is starting node.
					alongTrace = trace.get(fromEnd);
					hopOccurrences[fromEnd][alongTrace.index]++;
				}
			}
		}

		System.out.println("Sorting results.");
		for (int hops = 0; hops <= maxHops; hops++) {
			output = new File(containingPath + "probe-" + hops + ".dat");
			Arrays.sort(hopOccurrences[hops]);
			writeArray(hopOccurrences[hops], output);
		}
	}

	public static void simulate(Graph g, Random rand, int nRequests, int nIntersectTests,
			int routePolicy, int[] sinkPolsUsed, boolean printPairedMaxHTI, final String outputPath) {
		File outputFile = new File(outputPath);
		PrintStream stream = null;
		try {
			stream = new PrintStream(outputFile);
		} catch (IOException e) {
			System.err.println("Cannot open file \"" + outputPath + ": " + e);
			System.exit(1);
		}
		stream.println("Routing " + nRequests * nIntersectTests + " requests, policy " + routePolicy + " on network of size " + g.size() + ".");
		long startTime = System.currentTimeMillis();

		//TODO: This is only used in 2D capacity in one line: can make 1D?
		Request[][] requests = new Request[nRequests][nIntersectTests];

		//Separated into two loops so that the same set
		//of requests are generated for each routing
		//policy.
		//This is why the instances of Random must be the same each run: comparability between routing schemes.
		//TODO: Intersect means what?
		int[][] requestOrigins = new int[nRequests][nIntersectTests];
		for (int i = 0; i < nRequests; i++) {
			double l = rand.nextDouble();
			//create multiple requests with same target location but different origin
			for (int j = 0; j < nIntersectTests; j++) {
				Request r = new Request(l, rand, routePolicy, false, g);
				requests[i][j] = r;
				requestOrigins[i][j] = rand.nextInt(g.size());
			}
		}

		//Route all requests.
		//TODO: What does route do?
		for (int i = 0; i < nRequests; i++) {
			for (int j = 0; j < nIntersectTests; j++) {
				Request r = requests[i][j];
				SimpleNode origin = g.getNode(requestOrigins[i][j]);
				origin.route(r, null);
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
		for (int sp = 0; sp < sinkPolsUsed.length; sp++) {
			int sinkPolicy = sinkPolsUsed[sp];
			//TODO: "Hist" is history? Histogram?
			int sinkHistSize = 13;
			int[] sinkHist = new int[sinkHistSize];
			int[] lowUptimeSinkHist = new int[sinkHistSize];
			int[] highUptimeSinkHist = new int[sinkHistSize];
			int[] hopsHist = new int[50];

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

					int nHops = Math.min(hopsHist.length - 1, r.hopsTaken());
					hopsHist[nHops]++;
				}
			}

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

		//
		int[] maxHopsToIntersect = null;
		int[] pairedMaxHTI = null;
		int[][] hopsToSink = null;
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
				}
				/*Only valid if doing precise routing
				assert maxHopsToIntersect[i] >= 0;
				assert pairedMaxHTI[i] >= 0;
				*/
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
		}
		stream.println("Time taken (ms): " + (System.currentTimeMillis() - startTime));
		stream.print("Summary:\t" + g.size() + "\t" + g.nEdges() + "\t" + g.minDegree() + "\t");
		stream.print(g.maxDegree() + "\t" + Math.sqrt(g.degreeVariance()) + "\t" + g.meanLocalClusterCoeff() + "\t");
		stream.print(g.globalClusterCoeff() + "\t" + nRequests * nIntersectTests + "\t" + nPreciseRouted + "\t");
		stream.print(printArraySummary(decrements, false));
		if (nIntersectTests > 1) {
			System.out.print(printArraySummary(maxHopsToIntersect, false));
		}
		stream.println();
		stream.println();
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
		int pct80 = a[((int)(a.length * 0.8))];
		int pct90 = a[((int)(a.length * 0.9))];
		int pct97 = a[((int)(a.length * 0.97))];
		int pct99 = a[((int)(a.length * 0.99))];
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
		for (int i = 0; i < a.length; i++) ss += ((double)a[i])*((double)a[i]);
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
