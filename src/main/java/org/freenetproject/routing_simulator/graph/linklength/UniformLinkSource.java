package org.freenetproject.routing_simulator.graph.linklength;

import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;

import java.util.ArrayList;

/**
 * Generates link lengths with uniform / flat probability. Terrible distribution.
 */
public class UniformLinkSource extends LinkLengthSource {

	/**
	 * @see LinkLengthSource#LinkLengthSource(org.apache.commons.math3.random.RandomGenerator, java.util.ArrayList
	 */
	public UniformLinkSource(RandomGenerator random, ArrayList<SimpleNode> nodes) {
		super(random, nodes);
	}

	@Override
	public SimpleNode getPeer(SimpleNode from) {
		return closestTo(from, random.nextDouble() * 0.5);
	}
}