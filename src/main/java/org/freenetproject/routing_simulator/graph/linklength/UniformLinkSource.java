package org.freenetproject.routing_simulator.graph.linklength;

import org.freenetproject.routing_simulator.graph.node.SimpleNode;

import java.util.ArrayList;
import java.util.Random;

/**
 * Generates link lengths with uniform / flat probability. Terrible distribution.
 */
public class UniformLinkSource extends LinkLengthSource {

	/**
	 * @see LinkLengthSource#LinkLengthSource(java.util.Random, java.util.ArrayList)
	 */
	public UniformLinkSource(Random random, ArrayList<SimpleNode> nodes) {
		super(random, nodes);
	}

	@Override
	public SimpleNode getPeer(SimpleNode from) {
		return closestTo(from, random.nextDouble() * 0.5);
	}
}