package org.freenetproject.routing_simulator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.graph.degree.ConformingDegreeSource;
import org.freenetproject.routing_simulator.graph.degree.DegreeSource;
import org.freenetproject.routing_simulator.graph.degree.FixedDegreeSource;
import org.freenetproject.routing_simulator.graph.degree.PoissonDegreeSource;
import org.freenetproject.routing_simulator.graph.linklength.ConformingLinkSource;
import org.freenetproject.routing_simulator.graph.linklength.KleinbergLinkSource;
import org.freenetproject.routing_simulator.graph.linklength.LinkLengthSource;
import org.freenetproject.routing_simulator.graph.linklength.UniformLinkSource;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;

import static org.freenetproject.routing_simulator.util.File.readableFile;
import static org.freenetproject.routing_simulator.util.File.writableDirectory;
import static org.freenetproject.routing_simulator.util.File.writableFile;

/**
 * Parses arguments, determines any errors, and provides access to the values.
 *
 * Static methods parse and validate arguments; an instance stores the values on its fields.
 */
public class Arguments {

	public final boolean quiet, verbose, lattice, fastGeneration, runProbe, metropolisHastings, runRoute, includeLattice, bootstrap;
	public final int seed, networkSize, shortcuts, maxHops, nRequests;
	public final GraphGenerator graphGenerator;
	public final DataInputStream degreeInput, linkInput, graphInput;
	public final DataOutputStream degreeOutput, linkOutput, graphOutput;
	public final String outputProbe, outputRoute;
	public final FoldingPolicy foldingPolicy;
	public final RoutingPolicy routingPolicy;

	private final CommandLine cmd;

	private static final RoutingPolicy ROUTING_DEFAULT = RoutingPolicy.GREEDY;

	private Arguments(boolean quiet, boolean verbose, boolean lattice, boolean fastGeneration, boolean runProbe, boolean metropolisHastings, boolean runRoute, boolean includeLattice, boolean bootstrap,
	                  int seed, int networkSize, int shortcuts, int maxHops, int nRequests,
	                  GraphGenerator graphGenerator,
	                  DataInputStream degreeInput, DataInputStream linkInput, DataInputStream graphInput,
	                  DataOutputStream degreeOutput, DataOutputStream linkOutput, DataOutputStream graphOutput,
	                  String outputProbe, String outputRoute,
	                  FoldingPolicy foldingPolicy,
	                  RoutingPolicy routingPolicy,
	                  CommandLine cmd) {

		this.quiet = quiet;
		this.verbose = verbose;
		this.lattice = lattice;
		this.fastGeneration = fastGeneration;
		this.runProbe = runProbe;
		this.metropolisHastings = metropolisHastings;
		this.runRoute = runRoute;
		this.bootstrap = bootstrap;
		this.seed = seed;
		this.networkSize = networkSize;
		this.shortcuts = shortcuts;
		this.maxHops = maxHops;
		this.nRequests = nRequests;
		this.graphGenerator = graphGenerator;
		this.degreeInput = degreeInput;
		this.linkInput = linkInput;
		this.graphInput = graphInput;
		this.degreeOutput = degreeOutput;
		this.linkOutput = linkOutput;
		this.graphOutput = graphOutput;
		this.outputProbe = outputProbe;
		this.outputRoute = outputRoute;
		this.foldingPolicy = foldingPolicy;
		this.routingPolicy = routingPolicy;
		this.includeLattice = includeLattice;
		this.cmd = cmd;
	}

	public DegreeSource getDegreeSource(RandomGenerator random) {
		final DegreeSource degreeSource;

		if (cmd.hasOption("conforming-degree")) degreeSource = new ConformingDegreeSource(degreeInput, random);
		else if (cmd.hasOption("poisson-degree")) degreeSource = new PoissonDegreeSource(Integer.valueOf(cmd.getOptionValue("poisson-degree")));
		else if (cmd.hasOption("fixed-degree")) degreeSource = new FixedDegreeSource(Integer.valueOf(cmd.getOptionValue("fixed-degree")));
		else /* if (cmd.hasOption("sandberg-graph" || cmd.hasOption("supernode-graph") */ degreeSource = new FixedDegreeSource(0);

		return degreeSource;
	}

	public LinkLengthSource getLinkLengthSource(RandomGenerator random, ArrayList<SimpleNode> nodes) {
		final LinkLengthSource linkLengthSource;

		if (cmd.hasOption("conforming-link")) linkLengthSource = new ConformingLinkSource(linkInput, random, nodes);
		else if (cmd.hasOption("ideal-link")) linkLengthSource = new KleinbergLinkSource(random, nodes);
		else if (cmd.hasOption("flat-link")) linkLengthSource = new UniformLinkSource(random, nodes);
		else /* if cmd.hasOption("supernode-graph") */ linkLengthSource = null;

		return linkLengthSource;
	}

	private static Options generateOptions() {
		Options options = new Options();
		//TODO: Default values for arguments.
		//TODO: Can omit single-letter somehow? There are enough letters but not really a way to assign them in a way which makes sense.
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
		options.addOption("s", "size", true, "Number of nodes in the network. Currently ignored when using --degree unless --force-size is specified.");
		// TODO: Reinstate or possibly remove.
		options.addOption("f", "fast-generation", false, "If present, the simulator will assign locations with even spacing and, when using --ideal-link, take shortcuts to speed up graph generation.");
		options.addOption("G", "load-graph", true, "Path to load a saved graph from.");
		options.addOption("g", "save-graph", true, "Path to save a graph after simulation is run on it.");
		options.addOption("s", "sandberg-graph", true, "Generate a directed graph with an edge from x to x -1 mod N for all x = 0 ... N - 1 as described in early section 2.2.1 in \"Searching in a Small World.\" Takes the number of shortcuts to make; the paper specifies 1 shortcut.");
		options.addOption("l", "lattice", false, "Generate a graph with undirected lattice links, the given degree distribution, and the given link length distribution");
		options.addOption("s", "supernode-graph", false, "Generate a graph with all nodes connected to a single supernode.");

		//Graphs: link length distribution
		options.addOption("l", "ideal-link", false, "Kleinberg's ideal distribution: proportional to 1/d.");
		options.addOption("f", "flat-link", false, "Intentionally terrible distribution: uniformly random.");
		options.addOption("c", "conforming-link", true, "Distribution conforming to a file. Takes a path to a degree distribution file of the format \"[degree] [number of occurrences]\\n\"\"");

		//Graphs: degree distribution
		options.addOption("F", "fixed-degree", true, "All nodes are as close to the specified degree as practical.");
		options.addOption("C", "conforming-degree", true, "Distribution conforming to a file. Takes a path to a degree distribution file of the format \"[degree] [number of occurrences]\\n\"");
		options.addOption("i", "poisson-degree", true, "Distribution conforming to a Poisson distribution with the given mean.");

		//Simulations: Routing policies
		options.addOption("R", "route", true, "Simulate routing the given number of requests. Requires that --output-route, --fold-policy, and --output-hops be specified.");
		options.addOption("o", "output-route", true, "File to which routing information is output.");
		StringBuilder description = new StringBuilder("Path folding policy:");
		for (FoldingPolicy policy : FoldingPolicy.values()) description.append(" ").append(policy);
		options.addOption("P", "fold-policy", true,  description.toString());

		description = new StringBuilder("Routing policy used. Default is " + ROUTING_DEFAULT.name() +". Possible policies:");
		for (RoutingPolicy policy : RoutingPolicy.values()) description.append(" ").append(policy.name());
		options.addOption("r", "route-policy", true, description.toString());

		options.addOption("H", "output-hops", true, "Base filename to output hop histograms for each sink policy. Appended with -<policy-num> for each.");
		options.addOption("b", "bootstrap", false, "If specified, nodes which lose all their connections due to path folding will be connected to random nodes.");

		//Simulations: Probe distribution
		options.addOption("p", "probe", true, "Simulate running probes from random locations for the specified number of maximum hops. Requires that --output-probe be specified.");
		options.addOption("m", "metropolis-hastings", false, "If present, probes will be routed with Metropolis-Hastings correction. If not, peers will be selected entirely at random.");
		options.addOption("O", "output-probe", true, "Directory to which probe distribution is output as \"[node ID] [times seen]\\n\" for a reference of random selection from the whole and at each hop up to the specified maximum hops.");

		return options;
	}

	/**
	 * Parses command line arguments and validates them.
	 *
	 * @param args Arguments to parse.
	 * @return An Arguments instance with the parsed arguments, or null in the case of an error.
	 * @throws ParseException Failed to parse command line.
	 */
	public static Arguments parse(String[] args) throws ParseException {
		final Options options = generateOptions();
		final CommandLineParser parser = new GnuParser();
		final CommandLine cmd = parser.parse(options, args);

		if (cmd.hasOption("help")) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "java -jar simulator.jar", options );
			return null;
		}

		if (cmd.hasOption("version")) {
			System.out.println("Freenet Routing Simulator v 0.0.1-dev");
			return null;
		}

		//Check that required arguments are specified and that combinations make sense.
		if (cmd.hasOption("quiet") && cmd.hasOption("verbose")) {
			System.out.println("Quiet with verbose does not make sense.");
			return null;
		}
		if (cmd.hasOption("quiet") && !((cmd.hasOption("probe") && cmd.hasOption("output-probe")) || cmd.hasOption("output-degree") || cmd.hasOption("output-link") || cmd.hasOption("save-graph") || cmd.hasOption("route"))) {
			System.out.println("Simulation will produce no output: --quiet is specified, but not any option which outputs to a file.");
			return null;
		}

		int degreeOptions = 0;
		for (String option : new String[] { "fixed-degree", "conforming-degree", "poisson-degree" }) {
			if (cmd.hasOption(option)) degreeOptions++;
		}

		int linkOptions = 0;
		for (String option : new String[] { "ideal-link", "flat-link", "conforming-link" }) {
			if (cmd.hasOption(option)) linkOptions++;
		}

		if (degreeOptions > 1 || linkOptions > 1) {
			System.out.println("Graph cannot be generated with multiple methods at once.");
			return null;
		}

		/*
		 * Determine graph generation method to use, and if exactly one was specified correctly. Allowing more
		 * than one generation method to be specified would make it ambiguous which one was used.
		 */
		int validGeneration = 0;
		GraphGenerator graphGenerator = null;
		if (cmd.hasOption("load-graph")) {
			validGeneration++;
			graphGenerator = GraphGenerator.LOAD;
		}

		if (cmd.hasOption("supernode-graph")) {
			validGeneration++;
			graphGenerator = GraphGenerator.SUPER_NODE;
		}

		if (degreeOptions == 0 && linkOptions == 1 && cmd.hasOption("sandberg-graph")) {
			validGeneration++;
			graphGenerator = GraphGenerator.SANDBERG;
		}

		if (degreeOptions == 1 && linkOptions == 1) {
			validGeneration++;
			graphGenerator = GraphGenerator.STANDARD;
		}

		if (!cmd.hasOption("size") && graphGenerator != GraphGenerator.LOAD) {
			System.out.println("Network size not specified. (--size)");
			return null;
		}

		/*
		 * - Sandberg-graph requires link.
		 * - Load-graph does not require link or degree.
		 * - Without these two, degree and link must be specified; optionally lattice too.
		 */
		if (validGeneration == 0 || validGeneration > 1) {
			System.out.println("Either zero or too many graph generation methods specified.");
			System.out.println("Valid graph generators are:");
			System.out.println(" * --load-graph");
			System.out.println(" * --sandberg-graph with --*-link and --size");
			System.out.println(" * --*-degree, --*-link, --size, and optionally --lattice");
			System.out.println(" * --supernode-graph with --size and optionally --lattice");
			return null;
		}

		// By this point a single valid graph generator should be specified.
		assert graphGenerator != null;

		if (cmd.hasOption("route") && !cmd.hasOption("fold-policy")) {
			System.out.println("--route was specified, but not --fold-policy.");
			return null;
		}
		if (cmd.hasOption("probe") && !cmd.hasOption("output-probe")) {
			System.out.println("--probe was specified, but not --output-probe.");
			return null;
		}

		final FoldingPolicy foldingPolicy;
		if (cmd.hasOption("fold-policy")) {
			try {
				foldingPolicy = FoldingPolicy.valueOf(cmd.getOptionValue("fold-policy"));
			} catch (IllegalArgumentException e) {
				System.out.println("The folding policy \"" + cmd.getOptionValue("fold-policy") + "\" is invalid.");
				System.out.println("Possible values are:");
				for (FoldingPolicy policy : FoldingPolicy.values()) {
					System.out.println(policy.toString());
				}
				e.printStackTrace();
				return null;
			}
		} else {
			foldingPolicy = FoldingPolicy.NONE;
		}

		final RoutingPolicy routingPolicy;
		if (cmd.hasOption("route-policy")) {
			final String policy = cmd.getOptionValue("route-policy");
			try {
				routingPolicy = RoutingPolicy.valueOf(policy);
			} catch (IllegalArgumentException e) {
				System.out.println("The routing policy \"" + policy + "\" is invalid.");
				System.out.println("Possible values are:");
				for (RoutingPolicy policyName : RoutingPolicy.values()) {
					System.out.println(policyName.toString());
				}
				return null;
			}
		} else {
			routingPolicy = ROUTING_DEFAULT;
		}

		//Check for problems with specified paths.
		//Check if input files can be read.
		final DataInputStream degreeInput, linkInput, graphInput;
		try {
			degreeInput = readableFile("conforming-degree", cmd);
			linkInput = readableFile("conforming-link", cmd);
			graphInput = readableFile("load-graph", cmd);
		} catch (FileNotFoundException e) {
			return null;
		}

		//Check if output paths are directories that can be written to, and create them if they do not exist.
		final File probeOutputDirectory;
		if (cmd.hasOption("output-probe") && (probeOutputDirectory = writableDirectory(cmd.getOptionValue("output-probe"))) == null) return null;

		//Check that output files exist and are writable or can be created.
		final DataOutputStream degreeOutput, linkOutput, graphOutput;
		try {
			degreeOutput = writableFile("output-degree", cmd);
			linkOutput = writableFile("output-link", cmd);
			graphOutput = writableFile("save-graph", cmd);
		} catch (FileNotFoundException e) {
			return null;
		}

		final boolean quiet = cmd.hasOption("quiet");
		final boolean verbose = cmd.hasOption("verbose");
		final boolean lattice = cmd.hasOption("lattice");
		final boolean fastGeneration = cmd.hasOption("fast-generation");
		final int seed = cmd.hasOption("seed") ? Integer.valueOf(cmd.getOptionValue("seed")) : 0;
		final int networkSize = cmd.hasOption("size") ? Integer.valueOf(cmd.getOptionValue("size")) : 0;
		final int nRequests = cmd.hasOption("route") ? Integer.valueOf(cmd.getOptionValue("route")) : 0;
		final int maxHops = cmd.hasOption("probe") ? Integer.valueOf(cmd.getOptionValue("probe")) : 0;
		final int shortcuts = cmd.hasOption("sandberg-graph") ? Integer.valueOf(cmd.getOptionValue("sandberg-graph")) : 0;

		return new Arguments(quiet, verbose, lattice, fastGeneration, cmd.hasOption("probe"), cmd.hasOption("metropolis-hastings"), cmd.hasOption("route"), cmd.hasOption("include-lattice"), cmd.hasOption("bootstrap"),
		                     seed, networkSize, shortcuts, maxHops, nRequests,
		                     graphGenerator,
		                     degreeInput, linkInput, graphInput,
		                     degreeOutput, linkOutput, graphOutput,
		                     cmd.getOptionValue("output-probe"), cmd.getOptionValue("output-route"),
		                     foldingPolicy,
                                     routingPolicy,
		                     cmd);
	}
}
