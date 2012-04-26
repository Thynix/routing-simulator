import java.util.Random;

/**
 * This is a SimpleNode which also determines a desired degree given a WeightedDistribution.
 */
public class WeightedDegreeNode extends SimpleNode {
	public final int desiredDegree;

	WeightedDegreeNode(double location, boolean lowUptime, double pInstantReject, Random rand, WeightedDistribution distribution) {
		super(location, lowUptime, pInstantReject, rand);
		desiredDegree = distribution.randomValue();
	}

	public boolean atDegree() {
		return desiredDegree == degree();
	}
}
