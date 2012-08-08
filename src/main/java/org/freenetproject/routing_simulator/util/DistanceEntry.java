package org.freenetproject.routing_simulator.util;

import org.freenetproject.routing_simulator.graph.node.SimpleNode;

/**
 * Convenience class for ranking nodes by distance.
 */
public class DistanceEntry implements Comparable<DistanceEntry> {
	public final double distance;
	public final SimpleNode node;

	public DistanceEntry(double distance, SimpleNode node) {
		this.distance = distance;
		this.node = node;
	}

	@Override
	public int compareTo(DistanceEntry other) {
		return Double.compare(this.distance, other.distance);
	}
}