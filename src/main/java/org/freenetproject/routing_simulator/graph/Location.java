package org.freenetproject.routing_simulator.graph;

/**
 * Simple locations.
 */

public class Location {
	/**
	 * Return the minimum circular distance between two locations.
	 * Locations must be in [0,1).
	 *
	 * @param a First location
	 * @param b Second location
	 * @return Distance between a and b
	 */
	public static double distance(double a, double b) {
		if (a < 0.0 || a >= 1.0)
			throw new IllegalArgumentException("Invalid distance: " + a);
		if (b < 0.0 || b >= 1.0)
			throw new IllegalArgumentException("Invalid distance: " + b);
		double d = Math.min(Math.abs(a - b), 1.0 - Math.abs(b - a));
		assert d >= 0.0 && d <= 0.5: "Location error: "+a+", "+b+": "+d;
		return d;
	}
}
