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

	public RouteResult(boolean success) {
		this.success = success;
		this.disconnected = new ArrayList<SimpleNode>();
	}

	public RouteResult(boolean success, ArrayList<SimpleNode> disconnected) {
		this.success = success;
		this.disconnected = disconnected;
	}
}
