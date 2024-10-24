package util;

import java.util.HashMap; 
import java.util.Set; 

public class Interval implements Comparable<Interval>{
	public final int start;
	public final int stop;
	public final int strand;
	private String id= null;
	static private HashMap<String,Integer> bc2inx = new HashMap<String, Integer>();
	private HashMap<Integer,Integer> cov;
	private boolean unstranded=false;
	
	public Interval(int start,int stop, int strand) {
		this.start = start;
		this.stop = stop;
		this.strand = strand;
		this.cov = new HashMap<>();
		unstranded = strand == 0;
	}
	
	public Interval(int start,int stop, int strand,String id) {
		this(start,stop,strand);
		this.id = id;
	}
	
	public int getCov(String bc) {
		if(!bc2inx.containsKey(bc))
			return(0);
		return cov.get(bc2inx.get(bc));
	}
	
	public void addCov(String bc) {
		if(bc == null)
			bc = "";
		Integer inx = bc2inx.get(bc);
		if(inx == null) {
			inx = bc2inx.size();
			bc2inx.put(bc, inx);
		}
		Integer val = cov.get(inx);
		if(val == null)
			val = 0;
		cov.put(inx, val+1);
	}
	
	public String getId(){
		return id;
	}
	
	public Set<String> getBarcodes(){
		return bc2inx.keySet();
	}
	
	public int getTotalCov() {
		int r = 0;
		for(Integer k : cov.keySet())
			r += cov.get(k);
		return r;
	}

	public boolean equals(Object arg0) {
		if(arg0 == null || this.getClass() != arg0.getClass())
			return false;
		if(this == arg0)
			return true;
		if(arg0 instanceof Interval){
			Interval i = (Interval) arg0;
			return (unstranded || i.unstranded || i.strand == strand) && i.start == start && i.stop == stop && ((id == null && i.id == null) || id.equals(i.id)); 
		}
		return false;
	}
	
	/**
	 * if true, strand doesn't affect equals and compareTo methods;
	 * @param us
	 
	public void setUnstranded(boolean us){
		unstranded = us;
	}
	*/

	public int compareTo(Interval o) {
		if(!unstranded && !o.unstranded && o.strand != strand)
			return strand - o.strand;
		if(o.start != start)
			return start - o.start;
		return stop-o.stop;
	}

	public boolean overlap(Interval i){
		return i.start <= stop && i.stop>= start;
	}

	public int hashCode() {
		return start+stop;
	}
	
	public String toString() {
		return getClass().getName()+(id != null?" (id="+id+")":"")+": "+strand+":"+start+"-"+stop;
	}
	
	public int length(){
		return stop-start+1;
	}
	
}
