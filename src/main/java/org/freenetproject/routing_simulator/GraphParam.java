package org.freenetproject.routing_simulator;

/**
 * Class to specify graph parameters for Kleinberg etc graphs.
 */
public class GraphParam {
	public final int n;
	public final int p;
	public final int q;
	public final double pLowUptime;
	public final double pInstantReject;
	public final boolean fastGeneration;

	GraphParam(int n, int p, int q, double pLowUptime, double pInstantReject, boolean fastGeneration) {
		if (n <= 2*p + 2*q + 1)
			throw new IllegalArgumentException("Not enough nodes.");
		if (q < 0 || p < 0)
			throw new IllegalArgumentException("Must have non-negative outgoing edges.");
		if (pLowUptime < 0.0 || pLowUptime > 1.0)
			throw new IllegalArgumentException("Probabilities are in [0,1].");
		if (pInstantReject < 0.0 || pInstantReject > 1.0)
			throw new IllegalArgumentException("Probabilities are in [0,1].");

		this.n = n;
		this.p = p;
		this.q = q;
		this.pLowUptime = pLowUptime;
		this.pInstantReject = pInstantReject;
		this.fastGeneration = fastGeneration;
	}
}
