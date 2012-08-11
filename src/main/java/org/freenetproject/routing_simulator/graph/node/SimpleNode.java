package org.freenetproject.routing_simulator.graph.node;

import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.FoldingPolicy;
import org.freenetproject.routing_simulator.RouteResult;
import org.freenetproject.routing_simulator.RoutingPolicy;
import org.freenetproject.routing_simulator.graph.Location;
import org.freenetproject.routing_simulator.util.DistanceEntry;
import org.freenetproject.routing_simulator.util.lru.LRUQueue;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.ListIterator;

/**
 * A simple node model.  Has a location and a set of connections.
 */
public class SimpleNode {
	private final double location;
	private final ArrayList<SimpleNode> connections;
	private final int desiredDegree;

	private final RandomGenerator rand;

	/**Index of this node in the graph; purely for convenience, not used in any decision making.*/
	public final int index;

	final LRUQueue<SimpleNode> lruQueue;

	public void write(DataOutputStream out) throws IOException {
		out.writeDouble(location);
		out.writeInt(desiredDegree);
	}

	/**
	 * The last request ID to be routed by this node.
	 */
	private long lastRouted = -1;
	/**
	 * The current request ID.
	 */
	static long requestID = 0;

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SimpleNode)) return false;

		final SimpleNode other = (SimpleNode)o;

		if (this.degree() != other.degree()) return false;

		for (int i = 0; i < this.degree(); i++) {
			if (this.getConnections().get(i).index != other.getConnections().get(i).index) return false;
		}

		return this.index == other.index && this.location == other.location &&
		       this.desiredDegree == other.desiredDegree;
	}

	@Override
	public int hashCode() {
		return index;
	}

	public SimpleNode(DataInputStream in, int index, RandomGenerator rand) throws IOException {
		location = in.readDouble();
		desiredDegree = in.readInt();

		//Connections must be initialized later from the network view where other nodes are visible.
		connections = new ArrayList<SimpleNode>();
		this.index = index;
		lruQueue = new LRUQueue<SimpleNode>();
		this.rand = rand;
	}

	/**
	 * Public constructor.
	 *
	 * @param location The routing location of this node.
	 * @param rand Used for random numbers in decision making.
	 * @param desiredDegree Desired degree of the node.
	 */
	public SimpleNode(double location, RandomGenerator rand, int desiredDegree, int index) {
		if (location < 0.0 || location >= 1.0)
			throw new IllegalArgumentException("Location must be in [0,1).");

		this.location = location;
		connections = new ArrayList<SimpleNode>();
		this.rand = rand;
		this.index = index;
		lruQueue = new LRUQueue<SimpleNode>();
		if (desiredDegree < 1) this.desiredDegree = 1;
		else this.desiredDegree = desiredDegree;
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
	 * @param other node to consider.
	 * @return routing distance between this node and the other.
	 */
	public double distanceTo(SimpleNode other) {
		return distanceToLoc(other.getLocation());
	}

	/**
	 * Get the location of this node.
	 *
	 * @return Node location
	 */
	public double getLocation() {
		return location;
	}

	/**
	 * @return a peer which can be disconnected when path folding.
	 */
	SimpleNode disconnectCandidate() {
		final SimpleNode least = lruQueue.pop();
		lruQueue.pushLeast(least);
		return least;
	}

	/**
	 * Called to offer a connection between the specified peer and this one during path folding.
	 * The lowest peer in the LRU queue is dropped to make room.
	 *
	 * @param peer peer to consider a connection to.
	 * @param acceptanceRate probability of accepting the fold.
	 * @return whether the offered connection was accepted.
	 */
	private boolean offerPathFold(final SimpleNode peer, double acceptanceRate) {
		// If already at degree, don't fold if not accepted.
		if (atDegree() && rand.nextDouble() > acceptanceRate) return false;

		final SimpleNode least = peer.disconnectCandidate();

		// Initial degree sum over all nodes involved. This should remain invariant.
		final int initialDegree = degree() + least.degree() + peer.degree();

		peer.disconnect(least);

		// Peers added via path folding are added to the end.
		connect(peer);
		lruQueue.remove(peer);
		lruQueue.pushLeast(peer);
		peer.lruQueue.remove(this);
		peer.lruQueue.pushLeast(this);

		// Path folding should not change the total connection count.
		assert initialDegree == degree() + least.degree() + peer.degree();

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
	private static ArrayList<SimpleNode> successFreenet(final ArrayList<SimpleNode> nodeChain) {
		//TODO: Include HTL in the chain to test not path folding at high HTL.
		// Iterate uses previous(); first node is size() - 2: second-to-last as intended.
		ListIterator<SimpleNode> iterator = nodeChain.listIterator(nodeChain.size() - 1);

		// TODO: More general variable name - from?
		// Get final element.
		SimpleNode foldingFrom = nodeChain.get(nodeChain.size() - 1);

		// Promote successful peers in the LRU queue.
		while (iterator.hasPrevious()) {
			SimpleNode to = iterator.previous();
			assert to.isConnected(foldingFrom);
			to.lruQueue.push(foldingFrom);
			foldingFrom = to;
		}

		final ArrayList<SimpleNode> disconnected = new ArrayList<SimpleNode>();
		iterator = nodeChain.listIterator(nodeChain.size() - 1);
		// Fold starting from the second-to-last node.
		while (iterator.hasPrevious()) {
			final SimpleNode foldingTo = iterator.previous();
			// Do not path fold to self.
			if (foldingFrom == foldingTo) continue;
			// Do not path fold to a node which is already connected.
			if (foldingFrom.isConnected(foldingTo)) continue;
			//If the path fold is accepted, the one before the accepting one starts another fold.
			//Use 7% acceptance as rough result from my node.
			// If the fold is accepted and there was a candidate for disconnection, add it.
			final SimpleNode candidate = foldingFrom.disconnectCandidate();
			if (foldingTo.offerPathFold(foldingFrom, 0.07)) {
				if (candidate != null && candidate.degree() == 0) disconnected.add(candidate);
				if (iterator.hasPrevious()) {
					foldingFrom = iterator.previous();
				}
			}
		}
		return disconnected;
	}

	/**
	 * Rewires this node's shortcut to the endpoint. Expects every node to have only two outgoing connections:
	 * its lattice connection and its shortcut. Assumes the lattice connection is the first; shortcut the second.
	 * @param endpoint endpoint to fold to.
	 * @param foldingPolicy If SANDBERG, first two are lattice and not touched, and shortcuts are undirected.
	 *                      If SANDBERG_DIRECTED, first one is lattice and not touched, and shortcuts are directed.
	 *                      Other values are not accepted.
	 * @return the node disconnected by the fold. null if no node was disconnected or the fold did not occur.
	 */
	private SimpleNode offerShortcutFold(final SimpleNode endpoint, final double acceptanceRate, final FoldingPolicy foldingPolicy) {
		if (foldingPolicy != FoldingPolicy.SANDBERG && foldingPolicy != FoldingPolicy.SANDBERG_DIRECTED
		     && foldingPolicy != FoldingPolicy.SANDBERG_NO_LATTICE) {
			throw new IllegalArgumentException("Attempted to use shortcut folding with policy " + foldingPolicy);
		}
		// Do not path fold to self.
		if (endpoint == this) return null;
		// Do not path fold to a node which is already connected.
		if (this.connections.contains(endpoint)) return null;
		if (atDegree() && rand.nextDouble() < acceptanceRate) return null;
		// Always accept if not at degree.


		final int latticeLinks;
		if (foldingPolicy == FoldingPolicy.SANDBERG) latticeLinks = 2;
		else if (foldingPolicy == FoldingPolicy.SANDBERG_DIRECTED) latticeLinks = 1;
		else /*if(foldingPolicy == FoldingPolicy.SANDBERG_NO_LATTICE)*/ latticeLinks = 0;

		// A node should always have its lattice links.
		assert degree() >= latticeLinks;
		assert endpoint.degree() >= latticeLinks;

		// No shortcuts remain - cannot add this fold because there is no connection to drop.
		if (endpoint.degree() == latticeLinks) {
			return null;
		}

		final SimpleNode disconnected;

		final int initialDegree;

		if (foldingPolicy == FoldingPolicy.SANDBERG || foldingPolicy == FoldingPolicy.SANDBERG_NO_LATTICE) {
			// Endpoint drops connection.
			int disconnectedShortcut = rand.nextInt(endpoint.degree() - latticeLinks) + latticeLinks;
			disconnected = endpoint.connections.get(disconnectedShortcut);

			// Can't leave node that's being disconnected without its lattice links, if any.
			assert disconnected.degree() > latticeLinks;

			initialDegree = degree() + disconnected.degree() + endpoint.degree();

			endpoint.disconnect(disconnected);
			endpoint.connect(this);
		} else /*if (foldingPolicy == FoldingPolicy.SANDBERG_DIRECTED)*/ {
			int disconnectedShortcut = rand.nextInt(degree() - latticeLinks) + latticeLinks;
			disconnected = connections.get(disconnectedShortcut);

			initialDegree = degree() + disconnected.degree() + endpoint.degree();

			disconnectOutgoing(disconnected);
			connectOutgoing(endpoint);
		}

		// Total connection count should remain invariant.
		assert initialDegree == degree() + disconnected.degree() + endpoint.degree();

		return disconnected;
	}

	/**
	 * Offers a connection to the endpoint to every other node. Does not affect lattice links, the number of
	 * which is determined based on the folding policy. Number of lattice links is documented on each folding
	 * policy enum.
	 * @param nodeChain  Nodes which make up the path the request has followed. First element is the origin of the
	 *                   request; last is the endpoint.
	 * @param foldingPolicy Folds with either SANDBERG (Undirected and first two connections lattice) or
	 *                      SANDBERG_DIRECTED policy. (Directed and first connection lattice.)
	 * @see org.freenetproject.routing_simulator.FoldingPolicy
	 */
	private static ArrayList<SimpleNode> successSandberg(final ArrayList<SimpleNode> nodeChain, FoldingPolicy foldingPolicy) {
		final ListIterator<SimpleNode> iterator = nodeChain.listIterator(nodeChain.size() - 1);

		final SimpleNode endpoint;
		if (iterator.hasNext()) endpoint = iterator.next();
		else return new ArrayList<SimpleNode>();

		// Should have final element as endpoint.
		assert endpoint.equals(nodeChain.get(nodeChain.size() - 1));

		// Back up iterator to avoid attempting to fold endpoint to endpoint.
		iterator.previous();

		final ArrayList<SimpleNode> disconnectedNodes = new ArrayList<SimpleNode>();
		while (iterator.hasPrevious()) {
			SimpleNode disconnected = iterator.previous().offerShortcutFold(endpoint, 0.07, foldingPolicy);
			if (disconnected != null && disconnected.degree() == 0) disconnectedNodes.add(disconnected);
		}

		/*
		 * If a node loses its sole connection due to a fold, if it is present earlier on in the chain it could
		 * fold to the endpoint too and gain a connection again. This behavior does not have a real-world
		 * analog.
		 */
		Iterator<SimpleNode> disconnectIterator = disconnectedNodes.iterator();
		while (disconnectIterator.hasNext()) {
			if (disconnectIterator.next().degree() != 0) disconnectIterator.remove();
		}

		for (SimpleNode node : disconnectedNodes) assert node.degree() == 0;

		return disconnectedNodes;
	}

	/**
	 *
	 * @param nodeChain
	 * @param policy
	 * @return List of nodes which lost all their peers through folding.
	 */
	private static ArrayList<SimpleNode> success(final ArrayList<SimpleNode> nodeChain, FoldingPolicy policy) {
		// Can't fold if no nodes involved, (local node was closest right off) or if one node nowhere to fold to.
		if (nodeChain.size() < 2) return new ArrayList<SimpleNode>();
		switch (policy) {
		case NONE: return new ArrayList<SimpleNode>();
		case FREENET: return successFreenet(nodeChain);
		case SANDBERG_NO_LATTICE:
		case SANDBERG:
		case SANDBERG_DIRECTED: return successSandberg(nodeChain, policy);
		default: throw new IllegalStateException("Missing folding implementation for policy " + policy.name());
		}
	}

	/**
	 * Routes greedily from this node to the node at the target location. Stops when hopsToLive hits zero, or the
	 * local node is the closest to the target out of itself and all its neighbors. When that happens, it performs
	 * path folds with the specified policy.
	 *
	 * @param target Node to route to.
	 * @param hopsToLive Maximum number of additional hops.
	 * @param foldingPolicy Path folding policy to use on success.
	 * @return Routing was successful: the target location was reached.
	 *
	 */
	public RouteResult route(final SimpleNode target, final int hopsToLive, final RoutingPolicy routingPolicy, final FoldingPolicy foldingPolicy) {
		/*
		 * NOTE: This is static and package-local: not thread-safe! The simulator is as of this writing strictly
		 * single-threaded. The request ID could be an argument otherwise.
		 */
		requestID++;
		// TODO: Duplicate argument value determination between these methods: chain and target.
		switch (routingPolicy) {
			case GREEDY: return greedyRoute(target.getLocation(), hopsToLive, false, new greedy(), foldingPolicy, new ArrayList<SimpleNode>());
			case LOOP_DETECTION: return greedyRoute(target.getLocation(), hopsToLive, false, new loopDetection(), foldingPolicy, new ArrayList<SimpleNode>());
			case BACKTRACKING: return greedyRoute(target.getLocation(), hopsToLive, true, new loopDetection(), foldingPolicy, new ArrayList<SimpleNode>());
			default: throw new IllegalStateException("Routing for policy " + routingPolicy.name() + " not implemented.");
		}
	}

	private interface PeerSelector {
		public SimpleNode selectPeer(final double target, final SimpleNode from, final ArrayList<SimpleNode> chain);
	}

	private static class loopDetection implements PeerSelector {
		@Override
		public SimpleNode selectPeer(final double target, final SimpleNode from, final ArrayList<SimpleNode> chain) {
			SimpleNode next = from;
			ArrayList<DistanceEntry> distances = new ArrayList<DistanceEntry>(from.degree());
			for (SimpleNode peer : from.getConnections()) {
				distances.add(new DistanceEntry(peer.distanceToLoc(target), peer));
			}
			Collections.sort(distances);

			if (distances.size() > 0) {
				next = distances.get(0).node;

				while (next.lastRouted == requestID) {
					distances.remove(0);
					if (distances.size() > 0) next = distances.get(0).node;
					else {
						// No peers remaining that would not result in a loop.
						break;
					}
				}
			}

			return next;
		}
	}

	private static class greedy implements PeerSelector {
		@Override
		public SimpleNode selectPeer(double target, SimpleNode from, ArrayList<SimpleNode> chain) {
			SimpleNode next = from;
			double closest = from.distanceToLoc(target);

			// Check peer distances.
			for (SimpleNode peer : from.getConnections()) {
				double distance = peer.distanceToLoc(target);
				if (distance < closest) {
					next = peer;
					closest = distance;
				}
			}

			return next;
		}
	}

	private RouteResult greedyRoute(final double target, int hopsToLive, final boolean backtracking, final PeerSelector peerSelector, final FoldingPolicy foldingPolicy, final ArrayList<SimpleNode> chain) {
		if (hopsToLive <= 0) throw new IllegalStateException("hopsToLive must be positive. It is " + hopsToLive);

		/*
		 * Check whether the request reached its destination, which was selected from among node
		 * locations.
		 */
		if (this.getLocation() == target) {
			chain.add(this);
			return new RouteResult(true, success(chain, foldingPolicy), chain.size());
		}

		// Find node next node to route to.
		final SimpleNode next = peerSelector.selectPeer(target, this, chain);

		// Nowhere is closer or available, and this node is not the target one.
		if (next == this) return new RouteResult(false, chain.size());

		//TODO: Probabilistic decrement
		hopsToLive--;

		chain.add(this);
		if (backtracking) lastRouted = requestID;
		if (hopsToLive == 0) {
			return new RouteResult(true, success(chain, foldingPolicy), chain.size());
		} else {
			final RouteResult result = next.greedyRoute(target, hopsToLive, backtracking, peerSelector, foldingPolicy, chain);
			// If the routing did not succeed and did not use all remaining hops, backtrack if enabled.
			final int additionalHops = result.pathLength - chain.size();
			if (backtracking && !result.success && additionalHops < hopsToLive) {
				System.out.println("BACKTRACKING");
				return this.greedyRoute(target, hopsToLive - additionalHops, backtracking, peerSelector, foldingPolicy, chain);
			} else {
				return result;
			}
		}
	}

	/**
	 * Disconnect from a random peer, and connect to the given one. Must have at least one connection which can be
	 * dropped.
	 *
	 * @param peer node to connect to.
	 *
	 * @return Peer disconnected from.
	 */
	public SimpleNode swapConnections(SimpleNode peer) {
		// Must have at least one connection to drop one.
		assert degree() > 0;
		final SimpleNode disconnected = connections.get(rand.nextInt(connections.size()));

		final int initalDegree = degree() + disconnected.degree() + peer.degree();

		disconnect(disconnected);
		connect(peer);

		// Connection count should remain invariant.
		assert initalDegree == degree() + disconnected.degree() + peer.degree();

		return disconnected;
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
	 * Disconnect in both directions from a node which is already connected to this one.
	 * @param other node to disconnect from.
	 */
	public void disconnect(SimpleNode other) {
		this.disconnectOutgoing(other);
		other.disconnectOutgoing(this);
	}

	/**
	 * Removes the connection going from this node to the given one.
	 * @param other node to disconnect from.
	 */
	public void disconnectOutgoing(SimpleNode other) {
		if (other == this)
			throw new IllegalArgumentException("Cannot disconnect from self.");
		if (!isConnected(other))
			throw new IllegalArgumentException("Cannot disconnect: not connected.");

		connections.remove(other);
		lruQueue.remove(other);
	}

	/**
	 * @return Number of connections outgoing from this node.
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
	public ArrayList<SimpleNode> randomWalkList(int hops, boolean uniform, RandomGenerator rand) {
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
	public SimpleNode randomWalk(int hops, boolean uniform, RandomGenerator rand) {
		return randomWalkList(hops, uniform, rand).get(hops);
	}
}
