package org.freenetproject.routing_simulator.graph;

import org.freenetproject.routing_simulator.graph.node.SimpleNode;
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
}
