package org.freenetproject.routing_simulator.graph.degree;

public interface DegreeSource {
	/**
	 * @return degree conforming to the distribution.
	 */
	public int getDegree();
}