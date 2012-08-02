package org.freenetproject.routing_simulator.graph.linklength;

import org.apache.commons.math3.random.RandomGenerator;
import org.freenetproject.routing_simulator.graph.node.SimpleNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Generates links conforming to the 1/d Kleinberg distribution.
 */
public class KleinbergLinkSource extends LinkLengthSource {
	private final HashMap<SimpleNode, double[]> probabilityCDF;

	public KleinbergLinkSource(RandomGenerator random, ArrayList<SimpleNode> nodes) {
		super(random, nodes);
		probabilityCDF = new HashMap<SimpleNode, double[]>(nodes.size());
	}

	@Override
	public SimpleNode getPeer(final SimpleNode from) {
		/* Build probability array if it does not already exist.
		 *
		 * Find normalizing constant for this node - sum distance probabilities so that they
		 * they are in non-decreasing order and may be searched through to find the closest link.
		 * Note that this means here the probability is proportional to 1/distance.
		 * sumProb is a non-normalized CDF of probabilities by node index.
		 */
		if (!probabilityCDF.containsKey(from)) {
			final double[] sumProb = new double[nodes.size()];
			double norm = 0.0;
			for (int j = 0; j < nodes.size(); j++) {
				if (j != from.index) {
					norm += 1.0 / from.distanceTo(nodes.get(j));
				}
				sumProb[j] = norm;
				//CDF must be non-decreasing
				if (j > 0) assert sumProb[j] >= sumProb[j-1];
			}
			probabilityCDF.put(from, sumProb);
		}
		assert probabilityCDF.containsKey(from);

		/*
		 * sumProb is a CDF, so to weight by it pick a "Y value" and find closest index.
		 * norm is now the highest (and last) value in the CDF, so this is picking
		 * a distance probability sum and finding the closest node for that distance.
		 * Because there are more nodes which match values in highly represented domains
		 * (steeper in the CDF) a random value is more likely to be in those areas.
		 */
		final double[] sumProb = probabilityCDF.get(from);
		final double norm = sumProb[sumProb.length - 1];

		int idx;
		double x;
		do {
			x = random.nextDouble() * norm;
			assert x <= norm;
			idx = Arrays.binarySearch(sumProb, x);

			/*
			 * If such value is not actually present, as it might not be due to being
			 * floating point, use the index where it would be inserted:
			 * idx = -insertion point - 1
			 * insertion point = -1 - idx
			 * The insertion point would be the length of the array and thus out of bounds
			 * if all elements were less than it, but this will not happen as norm is the
			 * greatest element and nextDouble() is [0, 1). This does not mean it will not
			 * choose the greatest element as insertion point is the index of the first
			 * greater element.
			 */
			if (idx < 0) idx = -1 - idx;
		} while (idx == from.index);

		/*
		 * idx is index of the first greater element, but use the lesser if it is closer, unless that
		 * would mean the origin node connecting to itself.
		 */
		if (idx > 0 && Math.abs(x - sumProb[idx - 1]) < Math.abs(x - sumProb[idx]) && idx - 1 != from.index) idx--;

		// A node cannot connect to itself - do not offer it as a possibility.
		assert nodes.get(idx) != from;

		return nodes.get(idx);
	}
}
