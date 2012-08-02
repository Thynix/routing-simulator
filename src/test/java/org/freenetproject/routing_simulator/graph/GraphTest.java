package org.freenetproject.routing_simulator.graph;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.graph.degree.FixedDegreeSource;
import org.freenetproject.routing_simulator.graph.linklength.KleinbergLinkSource;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;
import org.testng.annotations.Test;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Test graph saving and loading.
 */
public class GraphTest {

	private final File temporary ;

	public GraphTest() {
		temporary = new File("temp");
	}

	/**
	 * @return A consistent, fresh randomness source.
	 */
	private RandomGenerator getRandom() {
		return new MersenneTwister(0);
	}

	private Graph generateKleinberg() {
		final RandomGenerator random = getRandom();
		final ArrayList<SimpleNode> nodes = Graph.generateNodes(100, random, false, new FixedDegreeSource(5));
		return Graph.connectGraph(nodes, random, new KleinbergLinkSource(random, nodes));
	}

	/**
	 * Deletes the destination file, if it exists, and writes the graph to it.
	 * @param graph Graph to write.
	 * @param destination File to write the graph to.
	 * @throws IOException
	 */
	private void writeToFile(final Graph graph, final File destination) throws IOException {
		assert !destination.exists() || destination.delete();
		final DataOutputStream outputStream = new DataOutputStream(new FileOutputStream(destination));
		graph.write(outputStream);
	}

	private Graph readFromFile(final File source) throws IOException {
		assert source.exists();
		final DataInputStream inputStream = new DataInputStream(new FileInputStream(source));
		return Graph.read(inputStream, getRandom());
	}

	/**
	 * Tests that the two graphs are exactly the same. The graphs must have:
	 * <ul>
	 *         <li>the same number of nodes.</li>
	 *         <li>the same number of edges.</li>
	 *         <li>nodes with the same locations at the same indexes.</li>
	 *         <li>edges which connect nodes at the same locations.</li>
	 * </ul>
	 *
	 * @param graph1 First graph to consider.
	 * @param graph2 Second graph to consider.
	 *
	 * @return True if and only if the graphs are exactly the same.
	 */
	private boolean equal(final Graph graph1, final Graph graph2) {
		if (graph1.size() != graph2.size()) return false;

		if (graph1.nEdges() != graph2.nEdges()) return false;

		// Graphs are known to be the same size; check that the same locations are at the same indexes.
		for (int i = 0; i < graph1.size(); i++) {
			if (graph1.getNode(i).getLocation() != graph2.getNode(i).getLocation()) return false;
		}

		/*
		 * Indexes and locations being the same too, now index is analogous to location.
		 * Check that the same nodes are connected.
		 */
		for (int i = 0; i < graph1.size(); i++) {
			final ArrayList<SimpleNode> connections1 = graph1.getNode(i).getConnections();
			final ArrayList<SimpleNode> connections2 = graph2.getNode(i).getConnections();

			if (!connections1.equals(connections2)) return false;
		}

		return true;
	}

	/**
	 * A graph is the same whether it has been generated this run or read from a written graph.
	 *
	 * @throws IOException Error writing to or reading from temporary file.
	 */
	@Test
	public void saveLoad() throws IOException {
		// TODO: Write to memory instead of a temporary file.
		final Graph written = generateKleinberg();
		writeToFile(written, temporary);

		final Graph read = readFromFile(temporary);

		assert equal(written, read);
	}

	/**
	 * A graph read from a written graph can be written again, and produces the same graph.
	 *
	 * @throws IOException Error writing to or reading from temporary file.
	 */
	@Test
	public void saveLoaded() throws IOException {
		final Graph original = generateKleinberg();
		writeToFile(original, temporary);

		final Graph firstRead = readFromFile(temporary);
		writeToFile(firstRead, temporary);

		final Graph secondRead = readFromFile(temporary);

		assert equal(firstRead, secondRead);
	}

	/**
	 * The equality check returns true for graphs which are the same, and false for those which are not.
	 */
	@Test
	public void testEquality() {
		final Graph one = generateKleinberg();
		final Graph two = generateKleinberg();

		assert equal(two, one);

		two.getNode(0).disconnect(two.getNode(0).getConnections().get(0));

		assert !equal(two, one);
	}
}
