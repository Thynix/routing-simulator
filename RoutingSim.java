import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

/**
 * Class to perform routing simulations on Graphs.
 */
public class RoutingSim {
	//E-series preferred numbers
	public static final int[] Esix =    {	10,	15,	22,	33,	47,	68	};
	public static final int[] Etwelve = {	10, 12, 15, 18, 22, 27, 33, 39, 47, 56, 68, 82	};

	//generic output format
	static DecimalFormat outputFormat = new DecimalFormat("0.000000");
	
	/**
	 * Main simulator program.  Generate a set of graphs of different
	 * parameters, run a set of requests on them, and print assorted
	 * stats.
	 *
	 * @param args Command-line arguments; not used.
	 */
	public static void main(String[] args) {
		int nRequests = 25000;
		int nIntersectTests = 2;
		int[] sinkPolsUsed = {0, 1};
		//TODO: What are these? What is a pol?
		int[] routePolsUsed = {3};
		int nTrials = 1;

		boolean printGraphStats = true;
		boolean printPairedMaxHTI = false;
		boolean verbose = false;
		boolean printIndivStats = true;
		boolean printAvgStats = false;

		Random rand;

		int[] graph_n;
		int[] graph_p;
		int[] graph_q;
		double[] graph_pLow;
		double[] graph_pInstant;
		boolean[][] graph_genMode;

		/*
		graph_n = {100, 330, 1000, 1500, 2200, 3300, 4700, 6800, 10000};
		graph_p = {0, 1};
		graph_q = {4, 5, 6, 7, 8};
		graph_pLow = {0.1};
		graph_pInstant = {0.};
		graph_genMode = {{true, true}, {true, false}, {false, false}};	//{evenSpacing, fastGeneration}
		*/

		/* Initialize each combination of graphParams.
		 * n: Network size.
		 * p: Local connections per node.
		 * q: Remote connections per node.
		 * pLow: Probability of low uptime. High/low uptime is used when tabulating sinks to categorize as high
		 *       or low uptime.
		 * pInstant: Probability that a request to route is instantly rejected.
		 * evenSpacing: True: locations are distributed uniformly throughout the network.
		 *              False: locations are uniformly random.
		 * fastGeneration: True: assumes even distribution for simpler connection generation.
		 *                 False: precise connection generation but involves more computation and may be slower.
		 */
		//Network size.
		graph_n = new int[] {10000};
		//Local connections per node.
		graph_p = new int[] {0};
		//Remote connections per node.
		graph_q = new int[] {6};

		int maxHops = 50;
		graph_pLow = new double[] {0.1};
		graph_pInstant = new double[] {0.0};
		graph_genMode = new boolean[][] {{false, false}};	//{evenSpacing, fastGeneration}

		GraphParam[] graphParam;
		graphParam = new GraphParam[graph_n.length * graph_p.length *
			graph_q.length * graph_pLow.length * graph_pInstant.length *
			graph_genMode.length];

		int graph_idx = 0;
		for (int n : graph_n) {
			for (int p : graph_p) {
				for (int q : graph_q) {
					for (double pLow : graph_pLow) {
						for (double pInstant : graph_pInstant) {
							for (boolean[] genMode : graph_genMode) {
								graphParam[graph_idx++] = new GraphParam(n, p, q, pLow, pInstant, genMode[0], genMode[1]);
							}
						}
					}
				}
			}
		}

		//TODO: What is this used for?
		double[][][] avgStats = new double[graphParam.length][13][nTrials];

		//Time tracking: report time taken for each graph setting if verbose; upon completion otherwise.
		long startTime = System.currentTimeMillis();
		long lastTime = startTime;

		System.out.println(	"Simulating " + graphParam.length + " distinct graph parameter sets, " + nTrials + " trials each.");
		System.out.println();

		if (!verbose && printIndivStats) {
			System.out.print("p\tq\tpLow\tpInst\tevenSpacing\tfastGeneration\t");
			Graph.printGraphStatsHeader();
			System.out.println();
		}

		//Run nTrials trials.
		for (int trial = 0; trial < nTrials; trial++) {
			System.out.print("Trial " + trial + "... ");
			if (verbose || printIndivStats) System.out.println();
			for (int graphIter = 0; graphIter < graphParam.length; graphIter++) {
				rand = new MersenneTwister(trial);
				GraphParam gp = graphParam[graphIter];

				if (printGraphStats) {
					if (verbose) {
						System.out.print("Generating 1d Kleinberg graph of " + gp.n + " nodes, with ");
						System.out.println("parameters p = " + gp.p + ", q = " + gp.q + ".");
						System.out.print("pLowUptime = " + gp.pLowUptime + ", pInstantReject = " + gp.pInstantReject);
						System.out.println(", evenSpacing = " + gp.evenSpacing + ", fastGeneration = " + gp.fastGeneration);
					} else if (printIndivStats) {
						System.out.print(gp.p + "\t" + gp.q + "\t" + gp.pLowUptime + "\t" +
								gp.pInstantReject + "\t" + gp.evenSpacing + "\t" + gp.fastGeneration + "\t");
					}
				}

				//TODO: Graph constructor should just take GramParam instead of all its parts. (As well as MersanneTwister)
				Graph g = Graph.generatePeerDistGraph(gp, rand, "../../stats/peerDist_1407.dat");
				if (printGraphStats && printIndivStats) g.printGraphStats(verbose);
				if (verbose) {
					System.out.println("Time taken (ms): " + (System.currentTimeMillis() - lastTime));
					lastTime = System.currentTimeMillis();
				}

				double[] indivStats = g.graphStats();
				for (int i = 0; i < 13; i++) {
					avgStats[graphIter][i][trial] = indivStats[i];
				}

				//TODO: foreach
				for (int rp = 0; rp < routePolsUsed.length; rp++) {
					rand = new MersenneTwister(trial);
					simulate(g, rand, nRequests, nIntersectTests, routePolsUsed[rp], sinkPolsUsed, printPairedMaxHTI, maxHops);
				}
				if (verbose || printIndivStats) System.out.println();
			}
			if (!verbose) {
				System.out.println("Time taken (ms): " + (System.currentTimeMillis() - lastTime));
				lastTime = System.currentTimeMillis();
			}
		}

		if (printAvgStats) {
			System.out.println("Average stats:");
			System.out.print("p\tq\tpLow\tpInst\tevenSpacing\tfastGeneration\t");
			Graph.printGraphStatsHeader();
			System.out.println();
			for (int i = 0; i < graphParam.length; i++) {
				GraphParam gp = graphParam[i];
				System.out.print(gp.p + "\t" + gp.q + "\t" + gp.pLowUptime + "\t" +
						gp.pInstantReject + "\t" + gp.evenSpacing + "\t" + gp.fastGeneration + "\t");
				//TODO: Why is this hardcoded to 13?
				for (int j = 0; j < 13; j++) {
					ArrayStats s = new ArrayStats(avgStats[i][j]);
					System.out.print("(" + outputFormat.format(s.mean()) + " " + outputFormat.format(s.stdDev()) + ")\t");
				}
				System.out.println();
			}
		}
		System.out.println("Total time taken (ms): " + (System.currentTimeMillis() - startTime));
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

	public static void probeDistribution(Graph g, Random rand, int maxHops) {
		final String containingPath = "occurenceDistribution/";
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
		//Find baseline for visibility by selecting nodes from the entire network at random.
		for (int i = 0; i < nTrials; i++) {
			for (int walk = 0; walk < g.size(); walk++) {
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
				//False: not using uniform routing: MH correction.
				trace = source.randomWalkList(maxHops, false, rand);
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
			int routePolicy, int[] sinkPolsUsed, boolean printPairedMaxHTI) {

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
		System.out.println("Requests precisely routed: " + nPreciseRouted);
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
			System.out.println("Sink count histograms for policy " + sinkPolicy + ":");
			System.out.println("n\tTotal\tlowU\thighU");
			for (int i = 0; i < sinkHistSize; i++) {
				System.out.println(i + "\t" + sinkHist[i] + "\t" + lowUptimeSinkHist[i] + "\t" + highUptimeSinkHist[i]);
			}
			System.out.println("Mean sinks: " + meanSinkCount);
			System.out.println("Mean low uptime sinks: " + meanLowUptimeSinkCount);
			System.out.println("Mean high uptime sinks: " + meanHighUptimeSinkCount);
			System.out.println();
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
		System.out.println("HTL Decrements:");
		System.out.print(printArraySummary(decrements, true));

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

		System.out.println("Routing distance by " + (byDecrements ? "decrements:" : "HTL:"));
		System.out.println((byDecrements ? "Dec" : "HTL") + "\tCount\tLogDist");
		for (int i = 0; i < logDistance.length; i++) {
			double meanLogDist = logDistance[i];
			if (logDistCount[i] > 0) {
				meanLogDist /= logDistCount[i];
			} else {
				assert meanLogDist == 0.0;
			}
			System.out.println(i + "\t" + logDistCount[i] + "\t" + meanLogDist);
		}
		System.out.println();

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
			System.out.println("Max hops to intersection:");
			System.out.print(printArraySummary(maxHopsToIntersect, true));
			if (printPairedMaxHTI) {
				Arrays.sort(pairedMaxHTI);
				System.out.println("Paired max hops to intersection:");
				System.out.print(printArraySummary(pairedMaxHTI, true));
			}
			System.out.println();
			for (int s = 0; s < sinkPolsUsed.length; s++) {
				System.out.println("Max hops to sink, policy " + sinkPolsUsed[s] + ":");
				Arrays.sort(hopsToSink[s]);
				System.out.print(printArraySummary(hopsToSink[s], true));
				System.out.println();
			}
		}
		System.out.println("Time taken (ms): " + (System.currentTimeMillis() - startTime));
		System.out.print("Summary:\t" + g.size() + "\t" + g.nEdges() + "\t" + g.minDegree() + "\t");
		System.out.print(g.maxDegree() + "\t" + Math.sqrt(g.degreeVariance()) + "\t" + g.meanLocalClusterCoeff() + "\t");
		System.out.print(g.globalClusterCoeff() + "\t" + nRequests * nIntersectTests + "\t" + nPreciseRouted + "\t");
		System.out.print(printArraySummary(decrements, false));
		if (nIntersectTests > 1) {
			System.out.print(printArraySummary(maxHopsToIntersect, false));
		}
		System.out.println();
		System.out.println();
		System.out.println();
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
