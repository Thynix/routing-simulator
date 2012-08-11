package org.freenetproject.routing_simulator;

public enum RoutingPolicy {
	/**
	 * Simple greedy routing to a node location: at each hop, route to the peer closest to the destination. The
	 * routing ends when the hopsToLive value hits zero or the local node is closer than all its peers. In a network
	 * with lattice links,the latter is sufficient to reach the global optimum.
	 */
	GREEDY,
	/**
	 * Greedy routing with loop detection. A peer will reject a request it has already received, and the routing
	 * node will route to the next-best. The routing ends when the hopsToLive value hits zero, or the local node
	 * has no peers left to route to.
	 */
	LOOP_DETECTION,
	/**
	 * Loop detection with backtracking: when there are no non-loop nodes remaining the request is passed up the
	 * chain for another routing attempt.
	 */
	BACKTRACKING
}
