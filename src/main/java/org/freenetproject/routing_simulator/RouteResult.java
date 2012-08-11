package org.freenetproject.routing_simulator;

import org.freenetproject.routing_simulator.graph.node.SimpleNode;

import java.util.ArrayList;

/**
 * Gives details on the results of routing a request.
 */
public class RouteResult {

	/**
	 * List of nodes caused to have zero degree by path folding.
	 */
	public final ArrayList<SimpleNode> disconnected;
	/**
	 * True if and only if the routing arrived at its exact target.
	 */
	public final boolean success;
	/**
	 * The length of the path taken, valid if and only if routing was successful.
	 */
	public final int pathLength;

	public RouteResult(boolean success, int pathLength) {
		this(success, new ArrayList<SimpleNode>(), pathLength);
	}

	public RouteResult(boolean success, ArrayList<SimpleNode> disconnected, int pathLength) {
		this.success = success;
		this.disconnected = disconnected;
		this.pathLength = pathLength;
	}
}
