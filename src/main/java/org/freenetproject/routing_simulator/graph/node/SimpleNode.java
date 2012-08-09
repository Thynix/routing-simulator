package org.freenetproject.routing_simulator.graph.node;

import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.FoldingPolicy;
import org.freenetproject.routing_simulator.RoutingPolicy;
import org.freenetproject.routing_simulator.graph.Location;
import org.freenetproject.routing_simulator.util.DistanceEntry;
import org.freenetproject.routing_simulator.util.lru.LRUQueue;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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

	private final LRUQueue<SimpleNode> lruQueue;

	public void write(DataOutputStream out) throws IOException {
		out.writeDouble(location);
		out.writeInt(desiredDegree);
	}

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
	 * Called to offer a connection between the specified peer and this one during path folding.
	 * The lowest peer in the LRU queue is dropped to make room. Refuses to fold to self or to a node which is
	 * already connected.
	 * @param peer peer to consider a connection to.
	 * @param acceptanceRate probability of accepting the fold.
	 * @return whether the offered connection was accepted.
	 */
	private boolean offerPathFold(final SimpleNode peer, double acceptanceRate) {
		// Do not path fold to self.
		if (peer == this) return false;
		// Do not path fold to a node which is already connected. TODO: Connections undirected - both should be true if either is.
		if (this.isConnected(peer) || peer.isConnected(this)) return false;
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
	 * Rewires this node's shortcut to the endpoint. Expects every node to have only two outgoing connections:
	 * its lattice connection and its shortcut. Assumes the lattice connection is the first; shortcut the second.
	 * @param endpoint endpoint to fold to.
	 * @param foldingPolicy If SANDBERG, first two are lattice and not touched, and shortcuts are undirected.
	 *                      If SANDBERG_DIRECTED, first one is lattice and not touched, and shortcuts are directed.
	 *                      Other values are not accepted.
	 * @return whether the fold was executed.
	 */
	private boolean offerShortcutFold(final SimpleNode endpoint, final double acceptanceRate, final FoldingPolicy foldingPolicy) {
		if (foldingPolicy != FoldingPolicy.SANDBERG && foldingPolicy != FoldingPolicy.SANDBERG_DIRECTED) {
			throw new IllegalArgumentException("Attempted to use shortcut folding with policy " + foldingPolicy);
		}
		//TODO: connections size is number of shortcuts + 1 (1 for the lattice link)
		//assert connections.size() == 2;
		// Do not path fold to self.
		if (endpoint == this) return false;
		// Do not path fold to a node which is already connected.
		if (this.connections.contains(endpoint)) return false;
		if (atDegree() && rand.nextDouble() < (1.0 - acceptanceRate)) return false;
		// More likely to accept if not at degree.
		else if (!atDegree() && rand.nextDouble() < (1.2 - acceptanceRate)) return false;


		final int latticeLinks;
		if (foldingPolicy == FoldingPolicy.SANDBERG) latticeLinks = 2;
		else /*if (foldingPolicy == FoldingPolicy.SANDBERG_DIRECTED)*/ latticeLinks = 1;

		// A node should always have its lattice links.
		assert degree() >= latticeLinks;
		assert endpoint.degree() >= latticeLinks;

		// No shortcuts remain - accept the incoming.
		if (degree() == latticeLinks) {
			// The time at which a node should be able to lose shortcuts
			assert foldingPolicy == FoldingPolicy.SANDBERG;
			connect(endpoint);
			return true;
		}

		int disconnectedShortcut = rand.nextInt(degree() - latticeLinks) + latticeLinks;
		final SimpleNode disconnected = connections.get(disconnectedShortcut);

		// Can't leave node that's being disconnected without its lattice links.
		assert disconnected.degree() > latticeLinks;

		if (foldingPolicy == FoldingPolicy.SANDBERG) {
			assert !this.isConnected(endpoint) && !endpoint.isConnected(this);
			assert this.isConnected(disconnected) && disconnected.isConnected(this);
			disconnect(disconnected);
			assert !this.isConnected(disconnected) && !disconnected.isConnected(this);
			connect(endpoint);
			assert this.isConnected(endpoint) && endpoint.isConnected(this);
		} else /*if (foldingPolicy == FoldingPolicy.SANDBERG_DIRECTED)*/ {
			disconnectOutgoing(disconnected);
			connectOutgoing(endpoint);
		}

		return true;
	}

	/**
	 * Offers a connection to the endpoint to every other node.
	 * @param nodeChain  Nodes which make up the path the request has followed. First element is the origin of the
	 *                   request; last is the endpoint.
	 * @param foldingPolicy Folds with either SANDBERG (Undirected and first two connections lattice) or
	 *                      SANDBERG_DIRECTED policy. (Directed and first connection lattice.)
	 */
	private static void successSandberg(final ArrayList<SimpleNode> nodeChain, FoldingPolicy foldingPolicy) {
		final ListIterator<SimpleNode> iterator = nodeChain.listIterator(nodeChain.size() - 1);

		final SimpleNode endpoint;
		if (iterator.hasNext()) endpoint = iterator.next();
		else return;

		// Should have final element as endpoint.
		assert endpoint.equals(nodeChain.get(nodeChain.size() - 1));

		// Back up iterator to avoid attempting to fold endpoint to endpoint.
		iterator.previous();

		while (iterator.hasPrevious()) {
			iterator.previous().offerShortcutFold(endpoint, 0.07, foldingPolicy);
		}
	}

	private static void success(final ArrayList<SimpleNode> nodeChain, FoldingPolicy policy) {
		// Can't fold if no nodes involved, (local node was closest right off) or if one node nowhere to fold to.
		if (nodeChain.size() < 2) return;
		switch (policy) {
		case NONE: return;
		case FREENET: successFreenet(nodeChain); break;
		case SANDBERG: successSandberg(nodeChain, FoldingPolicy.SANDBERG); break;
		case SANDBERG_DIRECTED: successSandberg(nodeChain, FoldingPolicy.SANDBERG_DIRECTED); break;
		}
	}

	/**
	 * Routes greedily from this node to the node at the target location. Stops when hopsToLive hits zero, or the
	 * local node is the closest to the target out of itself and all its neighbors. When that happens, it performs
	 * path folds with the specified policy.
	 *
	 * @param target Location to route to.
	 * @param hopsToLive Maximum number of additional hops.
	 * @param foldingPolicy Path folding policy to use on success.
	 *
	 */
	public void route(final double target, final int hopsToLive, final RoutingPolicy routingPolicy, final FoldingPolicy foldingPolicy) {
		switch (routingPolicy) {
			case GREEDY: greedyRoute(target, hopsToLive, false, foldingPolicy, new ArrayList<SimpleNode>()); break;
			//TODO: might be cleaner to use different routing method internally? Some degree of duplication but might read more nicely.
			case LOOP_DETECTION: greedyRoute(target, hopsToLive, true, foldingPolicy, new ArrayList<SimpleNode>()); break;
		}
	}

	private void greedyRoute(final double target, int hopsToLive, final boolean loopDetection, final FoldingPolicy foldingPolicy, final ArrayList<SimpleNode> chain) {
		if (hopsToLive <= 0) throw new IllegalStateException("hopsToLive must be positive. It is " + hopsToLive);
		// Find node closest to target. Start out assuming this node is the closest.
		SimpleNode next = this;
		if (loopDetection) {
			ArrayList<DistanceEntry> distances = new ArrayList<DistanceEntry>(this.degree());
			for (SimpleNode peer : this.getConnections()) {
				distances.add(new DistanceEntry(peer.distanceToLoc(target), peer));
			}
			Collections.sort(distances);

			if (distances.size() > 0) {
				next = distances.get(0).node;

				while (chain.contains(next)) {
					distances.remove(0);
					if (distances.size() > 0) next = distances.get(0).node;
					else {
						// No peers remaining that would not result in a loop.
						next = this;
						break;
					}
				}
			}
		} else {
			double closest = distanceToLoc(target);

			// Check peer distances.
			for (SimpleNode peer : connections) {
				double distance = peer.distanceToLoc(target);
				if (distance < closest) {
					next = peer;
					closest = distance;
				}
			}
		}

		/*
		 * Local node is the closest, (without loop detection) or has no peers which have not already been
		 * visited. (loop detection) Dead end.
		 */
		if (next == this) {
			/*
			 * Check whether the request reached its destination, which was selected from among node
			 * locations.
			 */
			if (loopDetection && this.getLocation() != target) return;
			chain.add(this);
			success(chain, foldingPolicy);
			return;
		}

		//TODO: Probabilistic decrement
		hopsToLive--;

		if (hopsToLive == 0) {
			if (loopDetection && this.getLocation() != target) return;
			chain.add(this);
			success(chain, foldingPolicy);
		} else {
			chain.add(this);
			next.greedyRoute(target, hopsToLive, loopDetection, foldingPolicy, chain);
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
