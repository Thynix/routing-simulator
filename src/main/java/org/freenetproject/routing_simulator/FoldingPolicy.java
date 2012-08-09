package org.freenetproject.routing_simulator;

public enum FoldingPolicy {
	/**
	 * Do not perform path folding.
	 */
	NONE,
	/**
	 * Path fold with 7% acceptance - each node along the chain, multiple times.
	 */
	FREENET,
	/**
	 * Path fold to endpoint only. Undirected network with lattice edges.
	 */
	SANDBERG,
	/**
	 * Path fold to endpoint only. Directed network with lattice edges.
	 */
	SANDBERG_DIRECTED,
	/**
	 * Path fold to endpoint only. Undirected network with no lattice edges.
	 */
	SANDBERG_NO_LATTICE
}
