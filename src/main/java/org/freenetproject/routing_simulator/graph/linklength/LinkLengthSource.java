package org.freenetproject.routing_simulator.graph.linklength;

import org.freenetproject.routing_simulator.graph.node.SimpleNode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

/**
 * Used to generate graphs conforming to a link length distribution. Not thread-safe.
 */
public abstract class LinkLengthSource {

	static class DistanceEntry implements Comparable<DistanceEntry> {
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

	final Random random;

	/**
	 * Stores a list of the distances to each other node for every node submitted with a peer query.
	 */
	private final HashMap<SimpleNode, ArrayList<DistanceEntry>> linkLengths;

	/**
	 * The nodes which make up the network being connected.
	 */
	final ArrayList<SimpleNode> nodes;

	/**
	 * @param random To make decisions.
	 * @param nodes Nodes which make up the network being wired.
	 */
	LinkLengthSource(Random random, ArrayList<SimpleNode> nodes) {
		this.random = random;
		this.nodes = nodes;
		/*
		 * TODO: If memory is a concern and a subclass is not using closestTo() it could set a flag to
		 * initialize linkLengths to a capacity of 0.
		 */
		this.linkLengths = new HashMap<SimpleNode, ArrayList<DistanceEntry>>(nodes.size());
	}

	/**
	 * Determines the node which provides closest to the desired link length out of the given nodes.
	 *
	 * @param from node the link is coming from.
	 * @param length desired link length.
	 * @return the node from the network providing the link closest to the specified length.
	 */
	SimpleNode closestTo(final SimpleNode from, final double length) {
		assert nodes.size() > 1;
		// Check if the link lengths have already been computed, and if not compute them.
		if (!linkLengths.containsKey(from)) {
			final ArrayList<DistanceEntry> distances = new ArrayList<DistanceEntry>(nodes.size() - 1);
			for (SimpleNode peer : nodes) {
				if (peer == from) continue;
				distances.add(new DistanceEntry(from.distanceTo(peer), peer));
			}
			Collections.sort(distances);
		}

		assert linkLengths.containsKey(from);
		final ArrayList<DistanceEntry> distances = linkLengths.get(from);
		int index = Collections.binarySearch(distances, new DistanceEntry(length, null));
		// Choose closest match. If not found index = -(insertion point) - 1, so insertion point = -1 - index.
		if (index < 0) index = -1 - index;
		if (index == distances.size()) index = distances.size() - 1;

		return distances.get(index).node;
	}

	/**
	 * Find a suitable peer which fits this distribution. Assumes that node locations do not change during
	 * generation.
	 * @param from Node to form a link from.
	 * @return a node for which the link length matches the link length distribution scheme.
	 */
	public abstract SimpleNode getPeer(SimpleNode from);
}
