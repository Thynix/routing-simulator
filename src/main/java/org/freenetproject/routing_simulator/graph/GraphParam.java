package org.freenetproject.routing_simulator.graph;

/**
 * Class to specify graph parameters for Kleinberg etc graphs.
 */
public class GraphParam {
	public final int n;
	public final double pLowUptime;
	public final double pInstantReject;
	public final boolean fastGeneration;

	public GraphParam(int n, double pLowUptime, double pInstantReject, boolean fastGeneration) {
		if (pLowUptime < 0.0 || pLowUptime > 1.0)
			throw new IllegalArgumentException("Probabilities are in [0,1]; got " + pLowUptime + " .");
		if (pInstantReject < 0.0 || pInstantReject > 1.0)
			throw new IllegalArgumentException("Probabilities are in [0,1]; got " + pInstantReject + ".");

		this.n = n;
		this.pLowUptime = pLowUptime;
		this.pInstantReject = pInstantReject;
		this.fastGeneration = fastGeneration;
	}
}
