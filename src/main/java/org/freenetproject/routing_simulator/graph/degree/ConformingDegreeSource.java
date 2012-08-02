package org.freenetproject.routing_simulator.graph.degree;

import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.util.WeightedDistribution;

import java.io.DataInputStream;

/**
 * Provides link lengths conforming with the distribution described by a list of link lengths.
 */
public class ConformingDegreeSource implements DegreeSource {
	private final WeightedDistribution distribution;

	public ConformingDegreeSource(DataInputStream input, RandomGenerator random) {
		this.distribution  = new WeightedDistribution(input, random);
	}

	@Override
	public int getDegree() {
		return distribution.randomValue();
	}
}
