package org.freenetproject.routing_simulator.graph.linklength;

import java.util.Random;

/**
 * Used to generate graphs conforming to a link length distribution.
 */
public interface LinkLengthSource {
	/**
	 * @return a desired link length for a connection determined by the link length distribution scheme.
	 * This will be attempted to be matched as closely as location distribution allows.
	 */
	//TODO: Don't take random - implementations can do that in the constructor if needed.
	public double getLinkLength(Random random);
}
