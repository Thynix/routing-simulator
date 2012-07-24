package org.freenetproject.routing_simulator.graph.node;

import org.freenetproject.routing_simulator.PathFolding;
import org.freenetproject.routing_simulator.util.lru.LRUQueue;
import org.freenetproject.routing_simulator.graph.Location;
import org.freenetproject.routing_simulator.Request;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.ListIterator;
import java.util.Random;

/**
 * A simple node model.  Has a location and a set of connections.
 */
public class SimpleNode implements Serializable {
	private double location;
	private ArrayList<SimpleNode> connections;
	public Request lastRequest;
	public boolean lowUptime;
	public int desiredDegree;

	private double pInstantReject;
	private Random rand;

	/**Index of this node in the graph; purely for convenience, not used in any decision making.*/
	public int index;

	private LRUQueue<SimpleNode> lruQueue;

	//Not in routeTo in order to reduce GC thrash
	double[] ignoreLoc = new double[1];

	private void writeObject(ObjectOutputStream out) throws IOException {
		out.writeDouble(location);
		out.writeInt(index);
		out.writeBoolean(lowUptime);
		out.writeDouble(pInstantReject);
		out.writeInt(desiredDegree);
	}

	private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
		location = in.readDouble();
		//Connections must be initialized later from the network view where other nodes are visible.
		connections = new ArrayList<SimpleNode>();
		index = in.readInt();
		lowUptime = in.readBoolean();
		pInstantReject = in.readDouble();
		lruQueue = new LRUQueue<SimpleNode>();
		ignoreLoc = new double[1];
		desiredDegree = in.readInt();
	}

	public void setRand(Random rand) {
		this.rand = rand;
	}

	/**
	 * Public constructor.
	 *
	 * @param location The routing location of this node.
	 * @param lowUptime Whether this is a low-uptime node
	 */
	public SimpleNode(double location, boolean lowUptime, double pInstantReject, Random rand, int candidate) {
		if (location < 0.0 || location >= 1.0)
			throw new IllegalArgumentException("Location must be in [0,1).");

		this.location = location;
		connections = new ArrayList<SimpleNode>();
		lastRequest = null;
		this.lowUptime = lowUptime;
		this.pInstantReject = pInstantReject;
		this.rand = rand;
		index = -1;
		lruQueue = new LRUQueue<SimpleNode>();
		if (candidate < 1) desiredDegree = 1;
		else desiredDegree = candidate;
	}

	/**
	 * @return True if the node is at (or above) its desired degree.
	 */
	public boolean atDegree() {
		return desiredDegree <= degree();
	}

	/**
	 * Get the routing distance to a location.
	 *
	 * @param l Location to compute distance of
	 * @return The circular routing distance
	 * @see org.freenetproject.routing_simulator.graph.Location distance
	 */
	public double distanceToLoc(double l) {
		return Location.distance(location, l);
	}

	/**
	 * Get the routing distance to a request location.
	 *
	 * @param r Request to compute distance of
	 * @return The circular routing distance
	 * @see Location distance
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

	/*private enum ConnectionState {
		CONNECTED,
		DISCONNECTED
	}*/

	/**
	 * Check that invariants hold between this node and the given. For debugging.
	 * @param other Other node to check this node's relationship with.
	 * @param connected Whether the two should be connected or not.
	 */
	/*private void checkInvariants(SimpleNode other, ConnectionState connected) {
		switch (connected) {
		case CONNECTED:
			assert lruQueue.contains(other) && other.lruQueue.contains(this);
			assert isConnected(other) && other.isConnected(this);
			break;
		case DISCONNECTED:
			assert !lruQueue.contains(other) && !other.lruQueue.contains(this);
			assert !isConnected(other) && !other.isConnected(this);
			break;
		}
	}*/

	/**
	 * Called to offer a connection between the specified peer and this one during path folding.
	 * The lowest peer in the LRU queue is dropped to make room. Refuses to fold to self or to a node which is
	 * already connected.
	 * @param peer peer to consider a connection to.
	 * @param acceptanceRate probability of accepting the fold.
	 * @return whether the offered connection was accepted.
	 */
	public boolean offerPathFold(final SimpleNode peer, double acceptanceRate) {
		// Do not path fold to self.
		if (peer == this) return false;
		// Do not path fold to a node which is already connected. TODO: Connections undirected - both should be true if either is.
		if (this.connections.contains(peer) || peer.connections.contains(this)) return false;
		if ((atDegree() || !peer.atDegree()) && rand.nextDouble() < (1.0 - acceptanceRate)) return false;

		// Disconnect from least, but restore invariant that both are in each other's LRU queues.
		if (atDegree()) {
		final SimpleNode least = lruQueue.pop();
		if (least == null) return false;
		lruQueue.pushLeast(least);

		//checkInvariants(least, ConnectionState.CONNECTED);
		disconnect(least);
		}
		//checkInvariants(least, ConnectionState.DISCONNECTED);

		// Peers added via path folding are added to the end.
		//checkInvariants(peer, ConnectionState.DISCONNECTED);
		connect(peer);
		lruQueue.remove(peer);
		lruQueue.pushLeast(peer);
		peer.lruQueue.remove(this);
		peer.lruQueue.pushLeast(this);
		//checkInvariants(peer, ConnectionState.CONNECTED);

		return true;
	}

	/**
	 * Folds with FREENET policy.
	 * Starting with the last node on the chain, nodes successively offers previous
	 * nodes a connection. If the connection is accepted the process starts again from the node previous to the
	 * accepting node.
	 * @param nodeChain Nodes which make up the path the request has followed. First element is the origin of the
	 *                  request; last is the endpoint.
	 */
	private static void successFreenet(final ArrayList<SimpleNode> nodeChain) {
		//TODO: Include HTL in the chain to test not path folding at high HTL.
		//TODO: Immediate peers should also be bumped up in the LRU upon success.
		// Iterate starting at the end.
		final ListIterator<SimpleNode> iterator = nodeChain.listIterator(nodeChain.size() - 1);

		SimpleNode foldingFrom;
		if (iterator.hasPrevious()) foldingFrom = iterator.previous();
		else return;

		//Start from the second-to-last node.
		while (iterator.hasPrevious()) {
			//If the path fold is accepted, the one before the accepting one starts another fold.
			//Use 7% acceptance as rough result from my node. TODO: Model Freenet's behavior for more accurate bootstrapping simulation.
			if (iterator.previous().offerPathFold(foldingFrom, 0.07) && iterator.hasPrevious()) {
				foldingFrom = iterator.previous();
			}
		}
	}

	/**
	 * Folds with SANDBURG policy.
	 * Offers a connection to the endpoint to every other node.
	 * @param nodeChain  Nodes which make up the path the request has followed. First element is the origin of the
	 *                   request; last is the endpoint.
	 */
	private static void successSandburg(final ArrayList<SimpleNode> nodeChain) {
		final ListIterator<SimpleNode> iterator = nodeChain.listIterator(nodeChain.size() - 1);

		final SimpleNode endpoint;
		if (iterator.hasPrevious()) endpoint = iterator.previous();
		else return;

		//Accept path fold with 1/path length probability.
		final double acceptProbability = 1 / nodeChain.size();

		while (iterator.hasPrevious()) {
			iterator.previous().offerPathFold(endpoint, acceptProbability);
		}
	}

	private static void success(final ArrayList<SimpleNode> nodeChain, PathFolding policy) {
		switch (policy) {
		case NONE: return;
		case FREENET: successFreenet(nodeChain);
		case SANDBERG: successSandburg(nodeChain);
		}
	}

	/**
	 * Route a request.  Routing policy is chosen by
	 * <code>r.routePolicy</code>.
	 *
	 * @param r The request to route
	 * @param origin the previous node in the chain.
	 * @param folding path folding policy to use on success
	 * @return Code indicating result
	 */
	public int route(final Request r, final SimpleNode origin, final PathFolding folding) {
		return route(r, origin, folding, new ArrayList<SimpleNode>());
	}

	/**
	 * Route a request.  Routing policy is chosen by
	 * <code>r.routePolicy</code>.
	 *
	 * @param r The request to route
	 * @param origin the previous node in the chain.
	 * @return Code indicating result
	 */
	public int route(final Request r, final SimpleNode origin, final PathFolding folding, final ArrayList<SimpleNode> previousNodes) {
		assert r.getHTL() <= Request.MAX_HTL;
		assert r.getHTL() > 0;
		assert !isLoop(r);
		//the originating node didn't instant reject the request
		if (origin != null && rand.nextDouble() < pInstantReject) {
			return RESULT_INSTANT_REJECT;
		}
		previousNodes.add(this);
		lastRequest = r;
		r.handleArrival(this);
		r.decrementHTL();
		if (r.getHTL() == 0) {
			//ran out of HTL; successfully routed.
			//r.sink(this, isSink(r, routable, origin, destNode));
			success(previousNodes, folding);
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

			int routeResult = destNode.route(r, this, folding, previousNodes);
			if (routeResult == RESULT_SUCCESS) {
				success(previousNodes, folding);
				return RESULT_SUCCESS;
			} else if (routeResult == RESULT_RNF || routeResult == RESULT_INSTANT_REJECT) {
				routable[dest] = false;
			} else {
				assert false;
			}

			if (r.getHTL() == 0) {
				//ran out of HTL while routing here or downstream (ie downstream RNF)
				success(previousNodes, folding);
				return RESULT_SUCCESS;
			}
		}
	}

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
	 * Check whether a connection from this node to the other exists.
	 *
	 * @param other Node to check for connectivity
	 * @return Whether a connection exists
	 */
	public boolean isConnected(SimpleNode other) {
		return connections.contains(other);
	}

	/**
	 * Form a bidirectional connection to the given node.
	 *
	 * @param other Node to connect with.
	 */
	public void connect(SimpleNode other) {
		this.connectOutgoing(other);
		other.connectOutgoing(this);
	}

	/**
	 * Form a one-way connection from this node to the given node.
	 * @param other Node to connect to.
	 */
	public void connectOutgoing(SimpleNode other) {
		if (other == this)
			throw new IllegalArgumentException("Cannot connect to self.");
		assert other != null;
		if (isConnected(other))
			throw new IllegalArgumentException("Cannot connect: already connected.");

		connections.add(other);
		lruQueue.push(other);
	}

	/**
	 * Disconnect from a node which is already connected to this one. As connections are
	 * undirected, the other node also disconnects.
	 * @param other node to disconnect from.
	 */
	public void disconnect(SimpleNode other) {
		if (other == this)
			throw new IllegalArgumentException("Cannot disconnect from self.");
		if (!isConnected(other) || !other.isConnected(this))
			throw new IllegalArgumentException("Cannot disconnect: not connected.");

		// Should be connected.
		//checkInvariants(other, ConnectionState.CONNECTED);

		connections.remove(other);
		lruQueue.remove(other);
		other.connections.remove(this);
		other.lruQueue.remove(this);

		// Now disconnected.
		//checkInvariants(other, ConnectionState.DISCONNECTED);
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
	 * Get a trace of a random walk starting here.
	 *
	 * @param hops Number of hops to walk
	 * @param uniform Whether to perform a uniform (normal) walk, or attempt to correct for high-degree bias
	 * @param rand Randomness source to use
	 * @return list of nodes along the way. The list's first element is the this node, its last is the endpoint.
	 */
	public ArrayList<SimpleNode> randomWalkList(int hops, boolean uniform, Random rand) {
		if (hops < 0) throw new IllegalArgumentException("Must have positive hops.");
		ArrayList<SimpleNode> list = new ArrayList<SimpleNode>();
		//System.out.println(hops + " HTL: At " + this.index);
		list.add(this);
		if (hops == 0) return list;

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
					/* Walk stayed on this node; had it ended here it would be the endpoint.
					 * Don't add extra if this is the last in the chain case though; will be added
					 * in the call for zero hops.
					 */
					if (hops > 0) list.add(this);
				}
			}
			//HTL ran out while looking for somewhere to route to: current node is endpoint.
			//It was already added above.
			if (hops == 0) {
				assert next == this;
				//Hops is currently zero - can't be negative so make it zero for the recursive call.
				hops = 1;
			}
		}
		list.addAll(next.randomWalkList(hops - 1, uniform, rand));
		return list;
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
		return randomWalkList(hops, uniform, rand).get(hops);
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
