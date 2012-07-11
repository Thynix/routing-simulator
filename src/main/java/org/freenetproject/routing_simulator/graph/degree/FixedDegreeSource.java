package org.freenetproject.routing_simulator.graph.degree;

/**
* Degree source which provides a single constant number.
*/
public class FixedDegreeSource implements DegreeSource {
	private final int degree;

	public FixedDegreeSource(int degree) {
		this.degree = degree;
	}

	@Override
	public int getDegree() {
		return degree;
	}
}
