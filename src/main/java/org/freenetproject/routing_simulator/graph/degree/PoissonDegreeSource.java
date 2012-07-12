package org.freenetproject.routing_simulator.graph.degree;

import org.apache.commons.math3.distribution.PoissonDistribution;

/**
 * Provides degrees conforming to a Poisson distribution with the given mean.
 */
public class PoissonDegreeSource implements DegreeSource {

	private final PoissonDistribution distribution;

	public PoissonDegreeSource(int mean) {
		distribution = new PoissonDistribution(mean);
	}

	@Override
	public int getDegree() {
		return distribution.sample();
	}
}
