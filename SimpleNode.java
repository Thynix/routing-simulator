import java.util.ArrayList;
import java.util.Random;

/**
 * A simple node model.  Has a location and a set of connections.
 */
public class SimpleNode {
	public static final int ROUTE_SUCCESS = 0;
	public static final int ROUTE_RNF = 1;

	private double location;
	private ArrayList<SimpleNode> connections;
	public Request lastRequest;
	public boolean lowUptime;

	private double pInstantReject;
	private Random rand;

	/**Index of this node in the graph; purely for convenience, not used in any decision making.*/
	public int index;

	/**
	 * Public constructor.
	 *
	 * @param location The routing location of this node.
	 * @param maxConnections The maximum number of connections this node can have.
	 * @param id The unique id of this node.
	 * @param lowUptime Whether this is a low-uptime node
	 */
	public SimpleNode(double location, boolean lowUptime, double pInstantReject, Random rand) {
		if (location < 0.0 || location >= 1.0)
			throw new IllegalArgumentException("Location must be in [0,1).");

		this.location = location;
		connections = new ArrayList<SimpleNode>();
		lastRequest = null;
		this.lowUptime = lowUptime;
		this.pInstantReject = pInstantReject;
		this.rand = rand;
		index = -1;
	}

	/**
	 * Get the routing distance to a location.
	 *
	 * @param l Location to compute distance of
	 * @return The circular routing distance
	 * @see Location.distance
	 */
	public double distanceToLoc(double l) {
		return Location.distance(location, l);
	}

	/**
	 * Get the routing distance to a request location.
	 *
	 * @param r Request to compute distance of
	 * @return The circular routing distance
	 * @see Location.distance
	 */
	public double distanceToLoc(Request r) {
		return distanceToLoc(r.getLocation());
	}

	/**
	 * Compute the minimum FOAF distance to a location.
	 * That is, the minimum distance of this node, or any of its peers, to
	 * that location; used for FOAF (Friend of a Friend) routing modes.
	 *
	 * @param l The location to compute the distance to
	 * @param ignoreLoc Any peers with locations exactly equal to these
	 * are ignored
	 * @return The minimum distance of the location to us or a peer
	 */
	public double minFOAFDist(double l, double[] ignoreLoc) {
		double d = distanceToLoc(l);
distloop:
		for (int i = 0; i < connections.size(); i++) {
			if (ignoreLoc != null) {
				double peerLoc = connections.get(i).location;
				for (int j = 0; j < ignoreLoc.length; j++) {
					if (peerLoc == ignoreLoc[j]) {
						continue distloop;
					}
				}
			}
			d = Math.min(d, connections.get(i).distanceToLoc(l));
		}
		return d;
	}

	/**
	 * Compute the minimum FOAF distance to a request.
	 *
	 * @param r The request to compute distance for
	 * @return Minimum FOAF distance
	 */
	public double minFOAFDist(Request r) {
		return minFOAFDist(r.getLocation(), null);
	}

	/**
	 * Compute the minimum FOAF distance to a request, with ignore
	 * locations.
	 *
	 * @param r The request to compute distance for
	 * @param ignoreLoc Set of locations to ignore
	 * @return Minimum FOAF distance
	 */
	public double minFOAFDist(Request r, double[] ignoreLoc) {
		return minFOAFDist(r.getLocation(), ignoreLoc);
	}

	/**
	 * Get the location of this node.
	 *
	 * @return Node location
	 */
	public double getLocation() {
		return location;
	}

	/** Set the location of this node.*/
	public void setLocation(double location) {
		if (location < 0.0 || location >= 1.0)
			throw new IllegalArgumentException("Invalid location");

		this.location = location;
	}

	/**
	 * Get whether this is a low uptime node.  This is used in sink
	 * decisions, but not routing decisions.
	 *
	 * @return True iff this is a low uptime node
	 */
	public boolean lowUptime() {
		return lowUptime;
	}

	/**
	 * Determine which sink policies this node will sink the data under.
	 * All policies only sink at low HTL, to protect privacy of inserters.
	 * Policy 0: Current Freenet (build 1234) policy.  Sink if we have a
	 * better location than all our high-uptime peers, ignoring routing
	 * decisions.
	 * Policy 1: evand's proposed change: sink if we have a better
	 * location than both the origin and destination, where low uptime
	 * (or null) nodes are considered to have bad locations.
	 *
	 * @param r The request in question
	 * @param routable Which of this node's peers are routable
	 * @param origin The node that this request came to us from; null for
	 * locally originated requests
	 * @param dest The node that this request will be routed to; null if
	 * the request will not be routed onward (HTL == 0, RNF, etc)
	 * @return An array indicating which policies sink the data
	 */
	protected boolean[] isSink(Request r, boolean[] routable, SimpleNode origin, SimpleNode dest) {
		boolean[] result = new boolean[Request.SINK_POLICIES];

		if (r.getHTL() > Request.MAX_STORE_HTL) {
			for (int i = 0; i < result.length; i++) result[i] = false;
			return result;
		}

		//current Freenet sink policy (build 1233)
		double myDist = distanceToLoc(r.getLocation());
		double minDist = 1.0;
		for (int i = 0; i < connections.size(); i++) {
			if (connections.get(i) == null) continue;
			double d = connections.get(i).distanceToLoc(r.getLocation());
			if (!connections.get(i).lowUptime() && d < minDist)
				minDist = d;
		}
		result[0] = myDist <= minDist;

		//proposed new sink policy
		minDist = 1.0;
		if (origin != null && !origin.lowUptime())
			minDist = Math.min(minDist, origin.distanceToLoc(r));
		if (dest != null && !dest.lowUptime())
			minDist = Math.min(minDist, dest.distanceToLoc(r));
		result[1] = myDist <= minDist;

		return result;
	}

	/**
	 * Whether this node has already handled this request.  Only the most
	 * recently routed request is tracked, so loop detection will not work
	 * if multiple requests are on the network at the same time.
	 *
	 * @param r Request to check for routing loops
	 * @return Whether this node has already handled the request
	 */
	public boolean isLoop(Request r) {
		return (r == lastRequest);
	}

	/** Return code indicating that routing succeeded. */
	public static final int RESULT_SUCCESS = 0;
	/** Return code indicating that routing was unable to find a destination. */
	public static final int RESULT_RNF = 1;
	public static final int RESULT_INSTANT_REJECT = 2;

	/**
	 * Route a request.  Routing policy is chosen by
	 * <code>r.routePolicy</code>.
	 *
	 * @param r The request to route
	 * @return Code indicating result
	 */
	public int route(Request r, SimpleNode origin) {
		assert r.getHTL() <= Request.MAX_HTL;
		assert r.getHTL() > 0;
		assert !isLoop(r);
		//the originating node didn't instant reject the request
		if (origin != null && rand.nextDouble() < pInstantReject) {
			return RESULT_INSTANT_REJECT;
		}
		lastRequest = r;
		r.handleArrival(this);
		r.decrementHTL();
		if (r.getHTL() == 0) {
			//ran out of HTL; successfully routed.
			//r.sink(this, isSink(r, routable, origin, destNode));
			return RESULT_SUCCESS;
		}

		boolean[] routable = new boolean[connections.size()];
		for (int i = 0; i < routable.length; i++) {
			routable[i] = true;
			if (origin != null && connections.get(i) == origin)
				routable[i] = false;
			if (connections.get(i) == null)
				routable[i] = false;
		}

		while (true) {
			int dest = routeTo(r, routable);
			SimpleNode destNode;

			if (dest == -1) {
				destNode = null;
			} else {
				destNode = connections.get(dest);
			}

			r.sink(this, isSink(r, routable, origin, destNode));

			if (dest == -1)
				return RESULT_RNF;
			assert destNode != null;

			if (destNode.isLoop(r)) {
				routable[dest] = false;
				r.decrementHTL();
				if (r.getHTL() == 0) {
					//r.sink(this, isSink(r, routable, origin, destNode));
					return RESULT_SUCCESS;
				}
				continue;
			}

			int routeResult = destNode.route(r, this);
			if (routeResult == RESULT_SUCCESS) {
				return RESULT_SUCCESS;
			} else if (routeResult == RESULT_RNF || routeResult == RESULT_INSTANT_REJECT) {
				routable[dest] = false;
			} else {
				assert false;
			}

			if (r.getHTL() == 0) {
				//ran out of HTL while routing here or downstream (ie downstream RNF)
				return RESULT_SUCCESS;
			}
		}
	}

	//Not in routeTo in order to reduce GC thrash
	double[] ignoreLoc = new double[1];

	/**
	 * Get the node that the request will be routed to.
	 *
	 * @param r The request to be routed
	 * @param routable Which of the connections we have should be considered routable
	 * @return The index of the node that is to be routed to next
	 */
	private int routeTo(Request r, boolean[] routable) {
		if (routable == null || routable.length != connections.size())
			throw new IllegalArgumentException();

		if (r.routePolicy >= 0 && r.routePolicy <= 6) {
			SimpleNode bestNode = null;
			int bestIdx = -1;
			double bestDist = 1.0;
			double bestLocalDist = 1.0;
			double myDist = distanceToLoc(r);

			if (r.routePolicy >= 3) {
				//loop fix
				if (ignoreLoc.length < connections.size() + 1) {
					ignoreLoc = new double[connections.size() + 1];
				}
				for (int i = 0; i < connections.size(); i++) {
					if (routable[i]) {
						ignoreLoc[i] = -1.0;
					} else {
						ignoreLoc[i] = connections.get(i).location;
					}
				}
				ignoreLoc[ignoreLoc.length - 1] = location;
			}

			for (int i = 0; i < connections.size(); i++) {
				if (!routable[i]) continue;
				SimpleNode target = connections.get(i);
				if (target == null) continue;
				double dist = 1.0;
				double localDist = target.distanceToLoc(r);
				double foafDist = 1.0;
				if (r.routePolicy == 0) {
				} else if (r.routePolicy == 1 || r.routePolicy == 2) {
					foafDist = target.minFOAFDist(r, null);
				} else if (r.routePolicy >= 3) {
					foafDist = target.minFOAFDist(r, ignoreLoc);
				} else {
					assert false;
				}

				if (r.routePolicy == 0) {
					dist = localDist;
				} else if (r.routePolicy == 1 || r.routePolicy == 2 || r.routePolicy == 3) {
					dist = foafDist;
				} else if (r.routePolicy == 4 || r.routePolicy == 5 || r.routePolicy == 6) {
					assert localDist >= foafDist;
					assert myDist != foafDist;
					if (myDist < foafDist) {
						dist = foafDist;
					} else {
						double p = 1.0;
						if (r.routePolicy == 4) {
							p = 3.0 / 6.0;
						} else if (r.routePolicy == 5) {
							p = 4.0 / 6.0;
						} else if (r.routePolicy == 6) {
							p = 5.0 / 6.0;
						} else {
							assert false;
						}
						dist = Math.pow(localDist, 1.0 - p) * Math.pow(foafDist, p);
					}
				} else {
					assert false;
				}
				if (dist < bestDist) {
					bestDist = dist;
					bestNode = target;
					bestIdx = i;
					bestLocalDist = localDist;
				} else if ((r.routePolicy >= 2) && dist == bestDist) {
					if (localDist < bestLocalDist) {
						bestLocalDist = localDist;
						bestIdx = i;
						bestNode = target;
					}
				} else {
				}
			}

			if (bestNode == null) return -1;
			assert bestDist <= 0.5 && bestDist >= 0.0;
			assert bestIdx >= 0;
			return bestIdx;
		} else {
			return -1;
		}
	}

	/**
	 * Check whether a connection exists to a node.
	 *
	 * @param other Node to check for connectivity
	 * @return Whether a connection exists
	 */
	public boolean isConnected(SimpleNode other) {
		return connections.contains(other);
	}

	/**
	 * Form a connection to a node.  Connections are undirected, so the
	 * complementary connection is also formed.
	 *
	 * @param other Node to connect to
	 */
	public void connect(SimpleNode other) {
		if (other == this)
			throw new IllegalArgumentException();
		if (connections.contains(other) || other.connections.contains(this))
			throw new IllegalArgumentException();

		connections.add(other);
		other.connections.add(this);
	}

	/**
	 * Compute the degree (number of connections) of this node.
	 *
	 * @return Degree of this node
	 */
	public int degree() {
		return connections.size();
	}

	public ArrayList<SimpleNode> getConnections() {
		return connections;
	}

	/**
	 * Count closed triplets centered on this node.
	 * See http://en.wikipedia.org/wiki/Clustering_coefficient
	 * Local clustering coefficient =
	 * closedTriplets() / (degree() * (degree() - 1) / 2)
	 *
	 * @return Count of closed triplets
	 */
	public int closedTriplets() {
		int cTrips = 0;
		int degree = connections.size();
		if (degree < 2) return 0;
		for (int i = 0; i < degree; i++) {
			for (int j = i + 1; j < degree; j++) {
				if (connections.get(i).isConnected(connections.get(j))) {
					cTrips++;
				}
			}
		}
		assert cTrips <= (degree * (degree - 1)) / 2;
		return cTrips;
	}

	/**
	 * Compute the local clustering coefficient.
	 * See http://en.wikipedia.org/wiki/Clustering_coefficient
	 *
	 * @return Local clustering coefficient
	 */
	public double localClusterCoeff() {
		int degree = connections.size();
		if (degree < 2) return 0.0;
		return ((double)closedTriplets())/((double)((degree * (degree - 1)) / 2));
	}

	/**
	 * Get the final node of a random walk starting here.
	 *
	 * @param hops Number of hops to walk
	 * @param uniform Whether to perform a uniform (normal) walk, or attempt to correct for high-degree bias
	 * @param rand Randomness source to use
	 * @return The final node of the walk
	 */
	public SimpleNode randomWalk(int hops, boolean uniform, Random rand) {
		if (hops < 0) throw new IllegalArgumentException("Must have positive hops.");
		if (hops == 0) return this;

		SimpleNode next;
		if (uniform) {
			//uniform walk: choose a neighbor with uniform probability
			next = connections.get(rand.nextInt(connections.size()));
		} else {
			//non-uniform walk: attempt to correct for high-degree node bias
			//Metroplis-Hastings approach
			next = this;
			while (hops > 0 && next == this) {
				next = connections.get(rand.nextInt(degree()));
				double beta = ((double)(degree())) / ((double)(next.degree()));
				//accept with probability beta
				if (rand.nextDouble() > beta) {
					//not accepted; decrement hops and try again
					next = this;
					hops--;
				}
			}
			if (hops == 0) {
				assert next == this;
				return this;
			}
		}
		if (hops == 1) {
			return next;
		} else {
			return next.randomWalk(hops - 1, uniform, rand);
		}
	}

	/**
	 * Attempt a M-H type location swap with the target.
	 *
	 * @param target The node to possibly swap locations with
	 * @return Whether the swap was actually performed
	 */
	public boolean attemptSwap(SimpleNode target) {
		double log = logSumLinkLengths(location) + target.logSumLinkLengths(target.location);
		log = log - logSumLinkLengths(target.location) + target.logSumLinkLengths(location);
		double acceptProb = Math.exp(log);
		return false;
	}

	/**
	 * Return the sum of the natural logarithms of the lengths this node's
	 * connections would be, if it had location l.
	 *
	 * @param l Hypothetical location for this node
	 * @return Log sum of link lengths
	 */
	public double logSumLinkLengths(double l) {
		double s = 0.0;
		for (int i = 0; i < connections.size(); i++) {
			s += Math.log(connections.get(i).distanceToLoc(l));
		}
		return s;
	}
}
