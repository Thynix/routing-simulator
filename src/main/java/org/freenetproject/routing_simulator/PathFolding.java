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
	 * Path fold to endpoint only.
	 */
	SANDBERG,
	/**
	 * Path fold to endpoint only with directed shortcut edges.
	 */
	SANDBERG_DIRECTED
}
