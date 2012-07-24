package org.freenetproject.routing_simulator;

public enum PathFolding {
	/**
	 * Do not perform path folding.
	 */
	NONE,
	/**
	 * Path fold with 7% acceptance - each node along the chain, multiple times.
	 */
	FREENET,
	/**
	 * Path fold only between origin and endpoint with probability inverse of path length.
	 */
	SANDBERG
}
