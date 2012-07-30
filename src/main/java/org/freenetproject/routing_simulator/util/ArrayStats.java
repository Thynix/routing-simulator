package org.freenetproject.routing_simulator.util;

import static java.lang.Math.sqrt;

/**
 * A class to compute stats about an array.  For speed, it caches a variety of
 * calculated properties of the array.  Therefore, if the array changes after
 * creation of the ArrayStats object, the invalidate() method must be called.
 */
public class ArrayStats {
	private final int[] iArray;
	private final double[] dArray;
	private final boolean isInt;

	private boolean valid;
	private double[] centralMoments;
	private double[] cumulants;
	private double mean;

	/**
	 * Construct a stats object for an int array.
	 *
	 * @param a The array to compute stats for
	 */
	public ArrayStats(int[] a) {
		this(a, null);
	}

	/**
	 * Construct a stats object for a double array.
	 *
	 * @param a The array to compute stats for
	 */
	public ArrayStats(double[] a) {
		this(null, a);
	}

	/**Private constructor handles initialization.*/
	private ArrayStats(int[] iArray, double[] dArray) {
		if (iArray == null && dArray == null)
			throw new NullPointerException();
		if (iArray != null && dArray != null)
			throw new IllegalArgumentException();
		valid = false;
		centralMoments = null;
		cumulants = null;
		mean = 0.0;
		this.isInt = iArray != null;
		this.iArray = iArray;
		this.dArray = dArray;
	}

	/**Calculate all cached data if required*/
	private void makeValid() {
		if (valid) return;
		computeMean();
		computeMoments();
		computeCumulants();
		assert !valid;
		valid = true;
	}

	/**Compute the mean*/
	private void computeMean() {
		if (valid) throw new IllegalStateException();
		mean = 0.0;
		if (isInt) {
			for (int anIArray : iArray) {
				mean += anIArray;
			}
			mean /= iArray.length;
		} else {
			for (double aDArray : dArray) {
				mean += aDArray;
			}
			mean /= dArray.length;
		}
	}

	/**Compute the moments about the mean*/
	private void computeMoments() {
		if (valid) throw new IllegalStateException();
		centralMoments = new double[6];
		centralMoments[0] = 1.0;
		centralMoments[1] = 0.0;
		int n = isInt ? iArray.length : dArray.length;
		for (int i = 0; i < n; i++) {
			double x = isInt ? iArray[i] : dArray[i];
			x -= mean;
			double m = x;
			for (int j = 2; j < centralMoments.length; j++) {
				m *= x;
				centralMoments[j] += m;
			}
		}
		for (int i = 2; i < centralMoments.length; i++) {
			centralMoments[i] /= n;
		}
	}

	/**Compute the cumulants*/
	private void computeCumulants() {
		if (valid) throw new IllegalStateException();
		cumulants = new double[centralMoments.length];

		assert cumulants.length == 6;	//TODO: generalized computation
		cumulants[0] = 0.0;
		cumulants[1] = mean;
		cumulants[2] = centralMoments[2];
		cumulants[3] = centralMoments[3];
		cumulants[4] = centralMoments[4] -
			3 * cumulants[2] * cumulants[2];
		cumulants[5] = centralMoments[5] -
			10 * cumulants[2] * cumulants[3];
	}

	/**
	 * Find the mean of the array.
	 *
	 * @return Arithmetic mean of the array
	 */
	public double mean() {
		makeValid();
		return mean;
	}

	/**
	 * Find the standard deviation of the array.
	 *
	 * @return Standard deviation of the array
	 */
	public double stdDev() {
		makeValid();
		return sqrt(cumulants[2]);
	}

	/**
	 * Find a standardized moment.
	 *
	 * @param n The moment to find
	 * @return The nth standardized moment of the array
	 */
	private double standardMoment(int n) {
		if (n < 0) throw new IllegalArgumentException();
		if (n == 0) return 1.0;
		if (n == 1) return 0.0;
		if (n == 2) return 1.0;
		makeValid();
		double sigma = stdDev();
		double sigmaN = sigma;
		for (int i = 1; i < n; i++) sigmaN *= sigma;
		return centralMoments[n] / sigmaN;
	}

	/**
	 * Find the skewness of the array.
	 *
	 * @return The skewness of the array
	 */
	public double skewness() {
		makeValid();
		return standardMoment(3);
	}
}
