package org.freenetproject.routing_simulator.graph.node;

import org.freenetproject.routing_simulator.graph.node.SimpleNode;

import java.util.Random;

/**
 * This is a SimpleNode which also determines a desired degree given a WeightedDistribution.
 */
public class WeightedDegreeNode extends SimpleNode {
	public final int desiredDegree;

	/**
	 * Generates a node with a desired degree which is guaranteed to be at least 1.
	 * @param location
	 * @param lowUptime
	 * @param pInstantReject
	 * @param rand
	 * @param candidate
	 */
	public WeightedDegreeNode(double location, boolean lowUptime, double pInstantReject, Random rand, int candidate) {
		super(location, lowUptime, pInstantReject, rand);
		if (candidate == 0) desiredDegree = 1;
		else desiredDegree = candidate;
	}

	/**
	 * @return True if the node is at (or above) its desired degree.
	 */
	public boolean atDegree() {
		return desiredDegree <= degree();
	}
}
