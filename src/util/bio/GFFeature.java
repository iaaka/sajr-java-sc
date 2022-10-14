package util.bio;

import java.util.ArrayList;
import java.util.HashMap;

import util.Util;

public class GFFeature{
	public final String seqname;
	public final String source;
	public final String feature;
	public final int start;
	public final int stop;
	public final double score;
	public final int strand;
	public final int phase;
	private final HashMap<String, String> attrs;
	
	public GFFeature(String seqname,String feature,int start,int stop,int strand,String source,double score,int phase){
		this.seqname = seqname;
		this.feature = feature;
		this.start = start;
		this.stop = stop;
		this.strand = strand;
		this.source = source;
		this.score = score;
		this.phase = phase;
		attrs = new HashMap<>();
	}
	
	public GFFeature(String seqname,String feature,int start,int stop,int strand){
		this(seqname,feature,start,stop,strand,"",0,0);
	}
	
	/**
	 * remember, it doesn't use attributes
	 */
	public int hashCode() {
		return (seqname+"|"+source+"|"+feature+"|"+start+"|"+stop+"|"+score+"|"+strand+"|"+phase).hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		GFFeature gff;
		if(obj instanceof GFFeature)
			 gff = (GFFeature) obj;
		else 
			return false;
		return seqname.equals(gff.seqname) && source.equals(gff.source) && feature.equals(gff.feature) && start == gff.start && stop == gff.stop && score == gff.score && strand == gff.strand && phase == gff.phase;
	}
	
	public GFFeature(String l) throws GFFException {
		String[] t = l.split("\t",-1);
		if(t.length != 9)
			throw new GFFException("Wrong gff file: there should be 9 columns: ",l);
		seqname = t[0];
		source = t[1];
		feature = t[2];
		start = Integer.parseInt(t[3]);
		stop = Integer.parseInt(t[4]);
		if(start> stop)
			throw new GFFException("Wrong gff file: starts cannot be greater than end",l);
		if(!t[6].equals(".") && !t[6].equals("-") && !t[6].equals("+"))
			throw new GFFException("Wrong gff file: strand should be in [.+-]",l);
		score = t[5].equals(".")?0:Double.parseDouble(t[5]);
		if(t[6].equals("."))
			strand = 0;
		else
			strand = t[6].equals("-")?-1:1;
		attrs = new HashMap<String, String>();
		if(t[7].equals("."))
			phase = -1;
		else
			phase = Integer.parseInt(t[7]);
		t[8] = Util.removeLeadingSpaces(t[8]);
		t = parseAttrs(t[8]);
		if(t.length%2 != 0)
			throw new GFFException("Wrong gff file: cannot parse attributes. There should be even number of elements separated by [\\s;=]+",l);
		for(int i=0;i<t.length;) {
			if(t[i].equals("")) {
				i++;
				continue;
			}
			attrs.put(t[i], t[i+1]);
			i+=2;
		}
	}
	
	public void addAttr(String n,String v){
		attrs.put(n, v);
	}
	
	/**
	 * Splits string by any sequence of "[\\s;=]+". all sequences within double quotes are considered as single token
	 * 
	 * @param a
	 * @return
	 */
	private String[] parseAttrs(String a) {
		String sep = "[\\s;=]+";
		char q = '"';
		ArrayList<String> r = new ArrayList<>();
		ArrayList<Integer> quotes = new ArrayList<>();		
		int t = -1;
		while(true) {
			t = a.indexOf(q, t+1);
			if(t == -1)
				break;
			quotes.add(t);
		}
		int prev = 0;
		for(int i=0;i<quotes.size()-1;i+=2) {
			if(quotes.get(i) > prev) {
				String f = a.substring(prev,quotes.get(i));
				String[] tmp =  f.split(sep);
				for(String m : tmp)
					if(m.length()>0)
						r.add(m);
			}
			r.add(a.substring(quotes.get(i)+1,quotes.get(i+1)));
			prev = quotes.get(i+1)+1;
		}
		if(a.length()>prev) {
			String f = a.substring(prev);
			String[] tmp =  f.split(sep);
			for(String m : tmp)
				if(m.length()>0)
					r.add(m);
		}
		//System.out.println("=====================");
		//for(String m : r)
		//	System.out.println("\t\t'"+m+"'");
		return r.toArray(new String[0]);
	}
	
	public String getAttr(String name) {
		return attrs.get(name);
	}
	
	public int getAttrInt(String name) {
		return Integer.parseInt(getAttr(name));
	}
	
	public String toString() {
		String r = seqname+"\t"+source+"\t"+feature+"\t"+start+"\t"+stop+"\t"+score+"\t"+(strand==0?".":(strand==1?"+":"-"))+"\t"+(phase==-1?".":phase)+"\t";
		for(String n : attrs.keySet())
			r = r + n+"="+attrs.get(n)+"; ";
		if(attrs.size() != 0)
			r = r.substring(0,r.length()-2);
		return r;
	}	
}
