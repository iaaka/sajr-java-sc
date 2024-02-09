package util;

import java.util.HashMap; 
import java.util.Set; 

public class Interval implements Comparable<Interval>{
	public final int start;
	public final int stop;
	public final int strand;
	private String id= null;
	private HashMap<String,Integer> cov;
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
		if(!cov.containsKey(bc))
			return(0);
		return cov.get(bc);
	}
	
	public void addCov(String bc) {
		if(bc == null)
			bc = "";
		if(!cov.containsKey(bc))
			cov.put(bc, 0);
		cov.put(bc, cov.get(bc)+1);
	}
	
	public String getId(){
		return id;
	}
	
	public Set<String> getBarcodes(){
		return cov.keySet();
	}
	
	public int getTotalCov() {
		int r = 0;
		for(String k : cov.keySet())
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
