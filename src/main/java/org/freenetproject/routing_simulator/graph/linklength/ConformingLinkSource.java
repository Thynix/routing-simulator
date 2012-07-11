package org.freenetproject.routing_simulator.graph.linklength;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class ConformingLinkSource implements LinkLengthSource {
	private final ArrayList<Double> lengths;
	public ConformingLinkSource(String filename) {
		lengths = new ArrayList<Double>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
			//TODO: Read all, put into ArrayList, get a link length selects from that.
			String line;
			//TODO: This seems like a C++ way of doing things. What's the Java way?
			while ( (line = reader.readLine()) != null) {
				//File format has link length as first value, separated by a space.
				lengths.add(Double.valueOf(line.split(" ")[0]));
			}
		} catch (FileNotFoundException e) {
			System.out.println(e);
			System.out.println("Unable to open file \"" + filename + "\".");
			System.exit(1);
		} catch (IOException e) {
			System.out.println(e);
			System.exit(2);
		}
	}

	@Override
	public double getLinkLength(Random random) {
		return lengths.get(random.nextInt(lengths.size()));
	}
}
