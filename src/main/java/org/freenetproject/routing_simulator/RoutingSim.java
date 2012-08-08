package org.freenetproject.routing_simulator;

import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.graph.Graph;
import org.freenetproject.routing_simulator.graph.linklength.LinkLengthSource;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;
import org.freenetproject.routing_simulator.util.ArrayStats;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Class to perform routing simulations on Graphs.
 */
public class RoutingSim {
	//generic output format
	private static final DecimalFormat outputFormat = new DecimalFormat("0.000000");

	/**
	 * Main simulator program.  Generate a set of graphs of different
	 * parameters, run a set of requests on them, and print assorted
	 * stats.
	 *
	 * @param args Command-line arguments; not used.
	 */
	public static void main(String[] args) throws ParseException {
		Arguments arguments = Arguments.parse(args);
		if (arguments == null) System.exit(1);

		//TODO: What is a sink policy? Looks like 2 is used as a hard-coded second dimension in Request, so this isn't currently a good configuration target.
		int[] sinkPolsUsed = {0, 1};
		PrintStream[] histogramOutput = new PrintStream[sinkPolsUsed.length];
		for (int aSinkPolsUsed : sinkPolsUsed) {
			//histogramOutput[i] = new PrintStream(writableFile(cmd.getOptionValue("output-hops")+ "-" + sinkPolsUsed[i]));
		}

		//TODO: What is this used for?
		//TODO: Hard-coded dimensions bad - clean this up.
		double[][][] avgStats = new double[1][13][1];

		if (!arguments.verbose) {
			System.out.println();
			System.out.print("p\tq\tpLow\tpInst\tevenSpacing\tfastGeneration\t");
			Graph.printGraphStatsHeader();
			System.out.println();
		}

		RandomGenerator rand = new MersenneTwister(arguments.seed);

		//Time tracking: report time taken for each graph setting if verbose; upon completion otherwise.
		long startTime = System.currentTimeMillis();
		long lastTime = startTime;

		// Load the graph; otherwise generate.
		Graph g;
		if (arguments.graphGenerator == GraphGenerator.LOAD) {
			g = Graph.read(arguments.graphInput, rand);
		} else {
			final ArrayList<SimpleNode> nodes = Graph.generateNodes(arguments.networkSize, rand, arguments.fastGeneration, arguments.getDegreeSource(rand));
			final LinkLengthSource linkLengthSource = arguments.getLinkLengthSource(rand, nodes);

			switch (arguments.graphGenerator) {
				case SANDBERG: g = Graph.connectSandberg(nodes, arguments.shortcuts, linkLengthSource); break;
				case SUPER_NODE: g = Graph.connectSuperNode(nodes, arguments.lattice); break;
				case STANDARD:
					if (arguments.lattice) g = Graph.connectGraphLattice(nodes, rand, linkLengthSource);
					else g = Graph.connectGraph(nodes, rand, linkLengthSource);
					break;
				default: throw new IllegalStateException("Missing implementation piece for graph generation method " + arguments.graphGenerator.name());
			}
		}

		if (!arguments.quiet) {
			g.printGraphStats(arguments.verbose);
			System.out.println("Generation took (ms): " + (System.currentTimeMillis() - lastTime));
			lastTime = System.currentTimeMillis();
		}

		double[] indivStats = g.graphStats();
		//TODO: 13 is defined because this doubles array has 13 different types of information. This should be a class.
		for (int i = 0; i < 11; i++) {
			avgStats[0][i][0] = indivStats[i];
		}

		if (arguments.runProbe) {
			// Re-initialize random number source so behavior here does not depend on previous usage, only the seed.
			rand = new MersenneTwister(arguments.seed);
			//Uniform probes if --metropolis-hastings is not specified.
			//TODO: Pass in checked directory.
			probeDistribution(g, rand, arguments.maxHops, arguments.quiet, arguments.verbose, arguments.outputProbe, !arguments.metropolisHastings);
		}

		if (arguments.runRoute) {
			rand = new MersenneTwister(arguments.seed);
			simulate(g, rand, arguments.nRequests, arguments.outputRoute, arguments.foldingPolicy, arguments.routingPolicy, histogramOutput);
		}

		if (arguments.degreeOutput != null) {
			int[] degrees = new int[g.maxDegree() + 1];
			for (int degree : g.degrees()) {
				degrees[degree]++;
			}
			for (int i = 0; i < g.maxDegree(); i++) {
				try {
					arguments.degreeOutput.write((i + " " + degrees[i] + "\n").getBytes());
				} catch (IOException e) {
					System.out.println("Unexpected error encoding string for degree output:");
					System.out.println(e);
					e.printStackTrace();
					return;
				}
			}
		}

		if (arguments.linkOutput != null) {
			ArrayList<Double> lengths = g.edgeLengths(arguments.includeLattice);
			//Output is intended for gnuplot CDF - second value is Y and should sum to 1.
			double normalized = 1.0/lengths.size();
			for (double length : lengths) {
				try {
					arguments.linkOutput.write((length + " " + normalized + "\n").getBytes());
				} catch (IOException e) {
					System.out.println("Unexpected error encoding string for link length output:");
					System.out.println(e);
					e.printStackTrace();
					return;
				}
			}
		}

		if (arguments.graphOutput != null) g.write(arguments.graphOutput);

		if (arguments.verbose) System.out.println();
		if (!arguments.quiet) {
			System.out.println("Time taken (ms): " + (System.currentTimeMillis() - lastTime));
		}

		if (!arguments.quiet) {
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

	private static void probeDistribution(Graph g, RandomGenerator rand, int maxHops, boolean quiet, boolean verbose, final String containingPath, boolean uniform) {
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

	private static void simulate(Graph g, RandomGenerator rand, int nRequests, final String outputPath,
	                             final FoldingPolicy foldingPolicy, final RoutingPolicy routingPolicy, final PrintStream[] histogramOutput) {
		/*File outputFile = new File(outputPath);
		PrintStream stream = null;
		try {
			stream = new PrintStream(outputFile);
		} catch (IOException e) {
			System.err.println("Cannot open file \"" + outputPath + ": " + e);
			System.exit(1);
		}
		stream.println("Routing " + nRequests + " requests on network of size " + g.size() + ".");*/
		long startTime = System.currentTimeMillis();


		for (int i = 0; i < nRequests; i++) {
			final SimpleNode origin = g.getNode(rand.nextInt(g.size()));
			/*
			 * It causes distortion to select among node locations for destinations as they may be less
			 * evenly distributed, but it allows determining if a request was routed exactly based on
			 * whether the target location is equal.
			 */
			final SimpleNode destination = g.getNode(rand.nextInt(g.size()));
			origin.route(destination.getLocation(), 50, routingPolicy, foldingPolicy);
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
		//stream.println("Time taken (ms): " + (System.currentTimeMillis() - startTime));
		//stream.print("Summary:\t" + g.size() + "\t" + g.nEdges() + "\t" + g.minDegree() + "\t");
		//stream.print(g.maxDegree() + "\t" + Math.sqrt(g.degreeVariance()) + "\t" + g.meanLocalClusterCoeff() + "\t");
		/*stream.print(g.globalClusterCoeff() + "\t" + nRequests * nIntersectTests + "\t" + nPreciseRouted + "\t");
		stream.print(printArraySummary(decrements, false));
		if (nIntersectTests > 1) {
			System.out.print(printArraySummary(maxHopsToIntersect, false));
		}
		stream.println();
		stream.println();*/
		//stream.println();
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

	private static double mean(int[] a) {
		double m = 0.0;
		for (int anA : a) m += anA;
		m /= a.length;
		return m;
	}

	private static double sumSquares(int[] a) {
		double ss = 0.0;
		for (int anA : a) ss += ((double) anA) * ((double) anA);
		return ss;
	}

	private static double variance(int[] a) {
		double ss = sumSquares(a);
		double m = mean(a);
		return (ss / a.length) - m * m;
	}

	private static double stdDev(int[] a) {
		return Math.sqrt(variance(a));
	}

	/**
	 * @param a The array to check.
	 * @return True if the array is sorted in non-decreasing order.
	 */
	private static boolean isSorted(int[] a) {
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

	private static double mean(double[] a) {
		double m = 0.0;
		for (double anA : a) m += anA;
		m /= a.length;
		return m;
	}

	private static double sumSquares(double[] a) {
		double ss = 0.0;
		for (double anA : a) ss += anA * anA;
		return ss;
	}

	private static double variance(double[] a) {
		double ss = sumSquares(a);
		double m = mean(a);
		return (ss / a.length) - m * m;
	}

	private static double stdDev(double[] a) {
		return Math.sqrt(variance(a));
	}

	private static boolean isSorted(double[] a) {
		for (int i = 0; i < a.length - 1; i++) if (a[i] > a[i+1]) return false;
		return true;
	}
}
