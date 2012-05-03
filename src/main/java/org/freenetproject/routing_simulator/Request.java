package org.freenetproject.routing_simulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Class to hold simple routable requests.
 * Has a location and an HTL, and tracks stats.
 */
public class Request {
	/** The maximum htl of a request.  All requests are started at this htl. */
	public static final int MAX_HTL = 18;

	/** The probability of decrementing the htl when the request is at max htl. */
	public static final double DEC_AT_MAX = 0.5;

	/** The probability of dropping a request when it is at htl == 1. */
	public static final double DROP_AT_MIN = 0.25;
	
	/**
	 * No storing above this.  Freenet caches requests at htl 16 and below,
	 * and caches / stores inserts at htl 15 and below.
	 */
	public static final int MAX_STORE_HTL = 15;

	/** Number of sink policies in use. */
	public static final int SINK_POLICIES = 2;

	private int htl;
	private double location;

	private boolean done;
	private boolean rnf;

	private ArrayList<SimpleNode> path;
	private boolean[][] sinks;
	private int[] pathHTL;
	private int[] pathDecrements;

	private Random rand;

	/**
	 * The policy that will be used to route this request.
	 * 
	 * @see SimpleNode.route
	 */
	public int routePolicy;

	/**
	 * If true, this request will route until it reaches its optimum location.
	 */
	public boolean routeFully;
	private Graph parentGraph;	//The graph this request is being routed on
	private int htlDecrements;

	/**
	 * Main constructor.
	 *
	 * @param location The location to route toward
	 * @param rand The java.util.Random to use for creating and routing
	 * this request
	 * @param routePolicy The routing policy to use for this request
	 */
	public Request(double location, Random rand, int routePolicy, boolean routeFully, Graph g) {
		if (location < 0.0 || location >= 1.0)
			throw new IllegalArgumentException("Invalid location: "+location);

		htl = MAX_HTL;
		this.location = location;
		done = false;
		rnf = false;
		path = new ArrayList<SimpleNode>();
		if (rand == null) {
			this.rand = new Random();
		} else {
			this.rand = rand;
		}
		sinks = new boolean[100][2];
		pathHTL = new int[100];
		pathDecrements = new int[100];
		for (int i = 0; i < pathHTL.length; i++) pathHTL[i] = -1;

		this.routePolicy = routePolicy;
		this.routeFully = routeFully;
		parentGraph = g;
		htlDecrements = 0;
	}

	/**
	 * Construct a request with randomized location.
	 *
	 * @param rand The java.util.Random to use for creating and routing
	 * this request
	 * @param routePolicy The routing policy to use for this request
	 */
	public Request(Random rand, int routePolicy, boolean routeFully, Graph g) {
		this(rand.nextDouble(), rand, routePolicy, routeFully, g);
	}

	/**
	 * Get the HTL of this request.
	 *
	 * @return HTL
	 */
	public int getHTL() {
		return htl;
	}

	/**
	 * Get the location being routed to
	 *
	 * @return The location to route to
	 */
	public double getLocation() {
		return location;
	}

	/** Probabilistically decrement the HTL.*/
	public void decrementHTL() {
		htlDecrements++;
		if (htl == MAX_HTL) {
			if (rand.nextDouble() < DEC_AT_MAX) {
				htl--;
			}
		} else if (htl == 1) {
			if (rand.nextDouble() < DROP_AT_MIN) {
				htl--;
			}
		} else {
			htl--;
		}
		//If routing to precise location, and not there yet, don't drop.
		if (routeFully) {
			if (htl == 0) {
				if (!parentGraph.preciseRoute(this))
					htl = 1;
			} else if (htl == 1) {
				if (parentGraph.preciseRoute(this))
					htl = 0;
			}
		}
	}

	/**
	 * Handle statistics on arrival at a new node.
	 *
	 * @param node The node this request has arrived at
	 */
	public void handleArrival(SimpleNode node) {
		if (path.indexOf(node) != -1) {
			throw new IllegalArgumentException("Illegal loop.");
		}
		if (node == null) {
			throw new NullPointerException();
		}
		path.add(node);
		int idx = path.size() - 1;
		if (idx < pathHTL.length) {
			pathHTL[idx] = htl;
		}
		if (idx < pathDecrements.length) {
			pathDecrements[idx] = htlDecrements;
		}
	}

	/**
	 * Record sink information when a node sinks a request.
	 *
	 * @param node The node sinking this request
	 * @param policies Array telling which policies will sink this request
	 */
	public void sink(SimpleNode node, boolean[] policies) {
		if (node == null)
			throw new NullPointerException();
		if (policies.length != SINK_POLICIES)
			throw new IllegalArgumentException();

		int idx = path.indexOf(node);
		if (idx == -1)
			throw new IllegalArgumentException("Node not in path.");
		if (idx >= sinks.length)
			return;
		for (int i = 0; i < SINK_POLICIES; i++) {
			if (policies[i]) sinks[idx][i] = true;
		}
	}

	/**
	 * The minimum routing distance achieved while routing this request.
	 *
	 * @return Minimum routing distance achieved
	 */
	public double minDist() {
		double d = 1.0;
		for (int i = 0; i < path.size(); i++) {
			SimpleNode n = path.get(i);
			d = Math.min(d, n.distanceToLoc(location));
		}
		assert d >= 0.0 && d <= 0.5;
		return d;
	}

	/**
	 * How many routing hops it took this node to reach its minimum
	 * routing distance.
	 *
	 * @return Hops to reach minimum routing distance
	 */
	public int hopsToMinDist() {
		double dist = 1.0;
		int hops = -1;
		for (int i = 0; i < path.size(); i++) {
			SimpleNode n = path.get(i);
			double d = n.distanceToLoc(location);
			if (d < dist) {
				dist = d;
				hops = i;
			}
		}
		assert dist >= 0.0 && dist <= 0.5;
		assert hops >= 0;
		return hops;
	}

	/**
	 * The HTL this request had when it reached its minimum routing
	 * distance.
	 *
	 * @return HTL of this request at its minimum routing distance
	 */
	public int htlAtMinDist() {
		int hops = hopsToMinDist();
		if (hops < pathHTL.length) {
			return pathHTL[hops];
		} else {
			return pathHTL[pathHTL.length - 1];
		}
	}

	/**
	 * The number of low-uptime nodes that sunk this request.
	 *
	 * @param policy The sink policy to check
	 * @return Number of low-uptime sinks
	 */
	public int nLowUptimeSinks(int policy) {
		int n = 0;
		for (int i = 0; i < path.size(); i++) {
			if (sinks[i][policy] && path.get(i).lowUptime()) {
				n++;
			}
		}
		return n;
	}

	/**
	 * The number of high-uptime nodes that sunk this request.
	 *
	 * @param policy The sink policy to check
	 * @return Number of high-uptime sinks
	 */
	public int nHighUptimeSinks(int policy) {
		int n = 0;
		for (int i = 0; i < path.size(); i++) {
			if (sinks[i][policy] && !path.get(i).lowUptime()) {
				n++;
			}
		}
		return n;
	}

	/**
	 * The number of nodes this request visited during routing.
	 *
	 * @return The number of routing hops taken
	 */
	public int hopsTaken() {
		return path.size();
	}

	/**
	 * How many times the HTL was decremented.  This is not quite the
	 * same as the number of hops taken.  This is not probabilistic.
	 *
	 * @return The number of times this request's HTL has been decremented.
	 */
	public int htlDecrements() {
		return htlDecrements;
	}

	/**
	 * The number of hops this request took before it intersected the path
	 * of the specified request.
	 *
	 * @param r The Request to check for intersection with
	 * @param sinksOnly Whether to check only nodes which were sinks for r
	 * @param sinkPolicy If sinksOnly is true, controls which sink policy to check
	 * @return Number of hops until intersection; -1 if no intersection
	 */
	public int hopsToIntersect(Request r, boolean sinksOnly, int sinkPolicy) {
		if (path.size() < 1) return -1;
		if (r.path.size() < 1) return -1;

		if (sinksOnly && (sinkPolicy < 0 || sinkPolicy >= SINK_POLICIES))
			throw new IllegalArgumentException("Invalid sink policy.");

		for (int i = 0; i < path.size(); i++) {
			if (i >= pathDecrements.length) return pathDecrements[pathDecrements.length - 1];
			SimpleNode pathNode = path.get(i);
			for (int j = 0; j < r.path.size(); j++) {
				if (pathNode == r.path.get(j)) {
					if (!sinksOnly || r.sinks[j][sinkPolicy])
						return pathDecrements[i];
				}
			}
		}
		return -1;
	}

	public int hopsToIntersect(Request r) {
		return hopsToIntersect(r, false, -1);
	}

	public int[] getPathDecrements() {
		return Arrays.copyOf(pathDecrements, Math.min(pathDecrements.length, path.size()));
	}

	public int[] getPathHTL() {
		return Arrays.copyOf(pathHTL, Math.min(pathHTL.length, path.size()));
	}

	public double[] getPathDistance() {
		double[] d = new double[path.size()];
		for (int i = 0; i < d.length; i++) {
			d[i] = path.get(i).distanceToLoc(location);
		}
		return d;
	}
}
