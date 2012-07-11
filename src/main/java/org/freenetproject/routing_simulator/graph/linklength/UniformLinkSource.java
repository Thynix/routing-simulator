package org.freenetproject.routing_simulator.graph.linklength;

import java.util.Random;

/**
 * Generates link lengths with uniform / flat probability. Terrible distribution.
 */
public class UniformLinkSource implements LinkLengthSource {
	@Override
	public double getLinkLength(Random random) {
		return random.nextDouble() * 0.5;
	}
}