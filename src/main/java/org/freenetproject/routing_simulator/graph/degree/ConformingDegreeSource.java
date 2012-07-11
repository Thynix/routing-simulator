package org.freenetproject.routing_simulator.graph.degree;

import org.freenetproject.routing_simulator.util.WeightedDistribution;

import java.util.Random;

/**
 * Provides link lengths conforming with the distribution described by a list of link lengths.
 */
public class ConformingDegreeSource implements DegreeSource {
	private final WeightedDistribution distribution;

	public ConformingDegreeSource(String filename, Random random) {
		this.distribution  = new WeightedDistribution(filename, random);
	}

	@Override
	public int getDegree() {
		return distribution.randomValue();
	}
}
