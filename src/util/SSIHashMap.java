package util;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import util.bio.Gene;


public class SSIHashMap {
	private HashSet<String> k1 = new HashSet<>();
	private HashSet<String> k2 = new HashSet<>();
	private HashMap<String,Integer> h = new HashMap<>();
	private final static String DELIMETER = ":";
	
	

	@SuppressWarnings("unchecked")
	public HashSet<String> getKeys1(){
		return (HashSet<String>) k1.clone();
	}
	
	@SuppressWarnings("unchecked")
	public HashSet<String> getKeys2(){
		return (HashSet<String>) k2.clone();
	}
	
	public void set(String key1,String key2,int v) {
		if(v == 0)
			return;
		k1.add(key1);
		k2.add(key2);
		h.put(key1+DELIMETER+key2, v);
	}
	
	public void addAll(SSIHashMap v) {
		for(String k : v.h.keySet())
			h.put(k, v.h.get(k));
		k1.addAll(v.k1);
		k2.addAll(v.k2);
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
		String type = h.values().iterator().next() instanceof Integer?"integer":"real";
		o.println("%%MatrixMarket matrix coordinate "+type+" general");
		o.println(ka1.length+" "+ka2.length+" "+h.size());
		for(String s : h.keySet()) {
			String[] ss = s.split(DELIMETER);
			o.println(ka1h.get(ss[0])+" "+ka2h.get(ss[1])+" "+h.get(s));
		}
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		System.out.println("afafgew");
		String[] ss = ("aad"+DELIMETER+"ad").split(DELIMETER);
		for(String s : ss)
			System.out.println(s);
//		SSIHashMap h = new SSIHashMap();
//		h.set("c1", "b1", 5);
//		h.set("c1", "b2", 2);
//		h.set("c2", "b2", 3);
//		h.set("c2", "b3", 7);
//		h.writeMM("test", new String[] {"c1","c2","c3","c4"}, new String[] {"b1","b2","b3","b4"});
	}
}
