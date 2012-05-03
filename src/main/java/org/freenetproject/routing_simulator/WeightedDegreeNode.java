package org.freenetproject.routing_simulator;

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
	 * @param distribution
	 */
	WeightedDegreeNode(double location, boolean lowUptime, double pInstantReject, Random rand, WeightedDistribution distribution) {
		super(location, lowUptime, pInstantReject, rand);
		int candidate = distribution.randomValue();
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
