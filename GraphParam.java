/**
 * Class to specify graph parameters for Kleinberg etc graphs.
 */
public class GraphParam {
	int n;
	int p;
	int q;
	double pLowUptime;
	double pInstantReject;
	boolean evenSpacing;
	boolean fastGeneration;

	GraphParam(int n, int p, int q, double pLowUptime, double pInstantReject, boolean evenSpacing, boolean fastGeneration) {
		this.n = n;
		this.p = p;
		this.q = q;
		this.pLowUptime = pLowUptime;
		this.pInstantReject = pInstantReject;
		this.evenSpacing = evenSpacing;
		this.fastGeneration = fastGeneration;
	}
}
