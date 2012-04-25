/**
 * Class to specify graph parameters for Kleinberg etc graphs.
 */
public class GraphParam {
	public final int n;
	public final int p;
	public final int q;
	public final double pLowUptime;
	public final double pInstantReject;
	public final boolean evenSpacing;
	public final boolean fastGeneration;

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
