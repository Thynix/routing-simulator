package org.freenetproject.routing_simulator.graph.node;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.FoldingPolicy;
import org.freenetproject.routing_simulator.RoutingPolicy;
import org.freenetproject.routing_simulator.graph.Graph;
import org.freenetproject.routing_simulator.graph.GraphTest;
import org.testng.annotations.Test;

/**
 * Tests node equality,
 */
public class NodeTest {

	@Test
	public void testEquality() {
		final Graph graph = GraphTest.generateKleinberg();
		final SimpleNode node = graph.getNode(0);
		final SimpleNode second = graph.getNode(1);

		// Reflexive.
		assert node.equals(node);

		// Symmetric within the class.
		assert !node.equals(second);
		assert !second.equals(node);

		// Symmetric with other classes.
		final Object testObject = new Object();
		assert !node.equals(testObject);
		assert !testObject.equals(node);

		// Not equal to null.
		assert !node.equals(null);
	}

	@Test
	public void testDisconnectCandidate() {
		RandomGenerator random = new MersenneTwister(0);

		// Manually make some nodes - easier to follow than graph generation.
		SimpleNode A = new SimpleNode(0.1, random, 0, 0);
		SimpleNode B = new SimpleNode(0.2, random, 0, 1);
		SimpleNode C = new SimpleNode(0.3, random, 0, 2);

		// Connect in a triangle.
		A.connect(C);
		A.connect(B);
		B.connect(C);

		// In A's LRU queue promote C; not B.
		A.lruQueue.push(C);

		/*
		 * The LRU queue should not be modified by disconnectCandidate().
		 */
		SimpleNode before[] = new SimpleNode[A.lruQueue.size()];
		A.lruQueue.toArrayOrdered(before);

		// B was not promoted, and so should be the least recently used.
		assert A.disconnectCandidate().equals(B);

		SimpleNode after[] = new SimpleNode[A.lruQueue.size()];
		A.lruQueue.toArrayOrdered(after);

		assert before.length == after.length;
		for (int i = 0; i < before.length; i++) assert before[i].equals(after[i]);
	}
}
