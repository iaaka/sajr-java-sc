package util;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import util.bio.Gene;


public class SSIHashMap {
	private HashMap<String,Integer> k2 = new HashMap<>();
	private ArrayList<String> inx2k2 = new ArrayList<String>();
	private HashMap<String,HashMap<Integer,Integer>> h = new HashMap<>();
	private final static String EMPTYSTRING = "-";
	
	
	public HashSet<String> getKeys1(){
		HashSet<String> res = new HashSet<>();
		res.addAll(h.keySet());
		return  res;
	}
	
	public HashSet<String> getKeys2(){
		HashSet<String> res = new HashSet<>();
		res.addAll(k2.keySet());
		return res;
	}
	
	private Integer get(String key1,Integer inx2) {
		return h.get(key1).get(inx2);
	}
	
	public Integer get(String key1,String key2) {
		return this.get(key1,k2.get(key2));
	}
	
	public void set(String key1,String key2,int v) {
		if(key1.equals(""))
			key1 = EMPTYSTRING;
		if(key2.equals(""))
			key2 = EMPTYSTRING;
		if(v == 0)
			return;
		Integer inx2 = k2.get(key2);
		if(inx2 == null) {
			inx2 = k2.size();
			k2.put(key2, inx2);
			inx2k2.add(key2);
		}
		if(!h.containsKey(key1)) {
			h.put(key1, new HashMap<>());
		}
		h.get(key1).put(inx2,v);
	}
	
	public void addAll(SSIHashMap v) {
		for(String k1 : v.h.keySet()) {
			for(Integer inx2 : v.h.get(k1).keySet()) {
				this.set(k1,v.inx2k2.get(inx2), v.get(k1,inx2));
			}
		}
	}
	
	public int size() {
		int size = 0;
		for(HashMap<Integer, Integer> v : h.values())
			size += v.size();
		return size;
	}
	
	public void writeMM(String out,String[] ka1,String[] ka2) throws FileNotFoundException {
		HashMap<String,Integer> ka1h = new HashMap<>();
		HashMap<String,Integer> ka2h = new HashMap<>();
		// key1
		PrintStream o = new PrintStream(out+"_key1.csv");
		int i = 1;
		for(String s : ka1) {
			o.println(s);
			ka1h.put(s, i);
			i++;
		}
		o.close();
		// key2
		o = new PrintStream(out+"_key2.csv");
		i = 1;
		for(String s : ka2) {
			o.println(s);
			ka2h.put(s, i);
			i++;
		}
		o.close();
		//values
		o = new PrintStream(out+".mtx");
		String type = "integer";
		o.println("%%MatrixMarket matrix coordinate "+type+" general");
		o.println(ka1.length+" "+ka2.length+" "+this.size());
		for(String k1 : h.keySet()) {
			if(ka1h.get(k1) == null) {
				o.close();
				throw new FileNotFoundException("unknown key1: '"+k1+"'");
			}
			for(Integer inx2 : h.get(k1).keySet()) {
				String k2 = inx2k2.get(inx2);
				o.println(ka1h.get(k1)+" "+ka2h.get(k2)+" "+this.get(k1,k2));
			}
		}
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		SSIHashMap h = new SSIHashMap();
		h.set("c1", "b1", 5);
		h.set("c1", "b2", 2);
		h.set("c2", "b2", 3);
		h.set("c2", "b3", 7);
		h.writeMM("test", new String[] {"c1","c2","c3","c4"}, new String[] {"b1","b2","b3","b4"});
	}
}
