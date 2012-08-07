package org.freenetproject.routing_simulator;

/**
 * Indicates which type of graph generator to use.
 */
public enum GraphGenerator {
	/**
	 * The graph will be loaded from a file.
	 */
	LOAD,
	/**
	 * The network will be generated with directed edges according to section 2.2.1 of Oskar Sandberg's
	 * "Searching in a Small World." Shortcut edges are added conforming to the link length source. In an extension
	 * of the Sandberg model more than one shortcut edge can be added per node. Lattice links are required and
	 * negatively oriented; from X to X - 1 mod N for all X = 0 .. N - 1.
	 */
	SANDBERG,
	/**
	 * The graph will be generated with undirected edges from every other node to the super node. Lattice links
	 * optional.
	 */
	SUPER_NODE,
	/**
	 * The graph will be generated with undirected edges and conform to the degree and link length sources.
	 * Lattice links optional.
	 */
	STANDARD
}
