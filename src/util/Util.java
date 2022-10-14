package util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;



public class Util {
	private Util(){};
	
	public static String join(int[] l,String sep){
		StringBuffer b = new StringBuffer();
		for(int e : l){
			b.append(e).append(sep);
		}
		return b.substring(0,Math.max(0,b.length()-sep.length())).toString();
	}
	
	public static int[] sum(int[] a,int[] b){
		int[] r = new int[a.length];
		for(int i=0;i<a.length;i++)
			r[i] = a[i]+b[i];
		return r; 
	}
	
	public static <E> String join(E[] l,String sep){
		StringBuffer b = new StringBuffer();
		for(E e : l){
			b.append(e).append(sep);
		}
		return b.substring(0,Math.max(0, b.length()-sep.length())).toString();
	}
	

	
	public static <E> String join(Iterable<E> l,String sep){
		StringBuffer b = new StringBuffer();
		for(E e : l){
			b.append(e.toString()).append(sep);
		}
		return b.substring(0,Math.max(0, b.length()-sep.length())).toString();
	}
	
	public static int max(int[] v) {
		int r = v[0];
		for(int i=1;i<v.length;i++)
			if(v[i] > r)
				r = v[i];
		return r;
	}
	
	public static String removeLeadingSpaces(String s) {
		int i = 0;
		for(;i<s.length();i++)
			if(s.charAt(i)!=' ')
				break;
		return s.substring(i);
	}
	
	public static int[] parseInt(String[] d){
		int[] r = new int[d.length];
		for(int i=0;i<d.length;i++)
			r[i] = Integer.parseInt(d[i]);
		return r;
	}
	
	public static int[] parseIntSequence(String s){
		ArrayList<Integer> r = new ArrayList<Integer>();
		String[] ss = s.split(",");
		for(String i : ss){
			if(i.contains(":")){
				String[] range = i.split(":");
				if(range.length != 2)
					throw new RuntimeException("Wrong format of Integer sequence: '"+s+"'");
				for(int j = Integer.parseInt(range[0]);j<= Integer.parseInt(range[1]);j++)
					r.add(j);
			}else
				r.add(Integer.parseInt(i));
		}
		return Util.toArray(r);
	}
	
	/**]
	 * 
	 * @param a
	 * @param b
	 * @return hashset that contains elements from a that do not exist in b
	 */
	public static <E> HashSet<E> diff(Set<E> a,Set<E> b){
		HashSet<E> r = new HashSet<>();
		for(E e : a) {
			if(!b.contains(e))
				r.add(e);
		}
		return r;
	}
	
	public static HashMap<String,String> readFasta(String f) throws IOException{
		HashMap<String, String> r = new HashMap<>();
		BufferedReader i = new BufferedReader(new FileReader(f));
		StringBuffer seq = null;
		String name = null;
		for(String l=i.readLine();l != null;l = i.readLine()){
			if(l.length() == 0)
				continue;
			if(l.charAt(0) == '>'){
				if(name != null){
					r.put(name, seq.toString());
				}
				name = l.substring(1);
				seq = new StringBuffer(seq==null?10000:seq.length());
			}else
				seq.append(l);
		}
		r.put(name, seq.toString());
		i.close();
		return r;
	}
	
	public static int[] getMapIntervals(SAMRecord r) {
		ArrayList<int[]> res = new ArrayList<int[]>();
		int start = r.getAlignmentStart();
		for(CigarElement c : r.getCigar().getCigarElements()) {
			if(c.getOperator().consumesReferenceBases()) {
				if(!c.getOperator().equals(CigarOperator.N))
					res.add(new int[] {start,start+c.getLength()-1});
				start += c.getLength();
			}
		}
		
		//join neighbor blocks
		res = union(res);
		int[] fr = new int[res.size()*2];
		for(int i=0;i<res.size();i++) {
			fr[i*2] = res.get(i)[0];
			fr[i*2+1] = res.get(i)[1];
		}
		return fr;
	}
	
	public static <E> HashSet<E> intersect(Collection<E> a,Iterable<E> b){
		HashSet<E> r = new HashSet<>();
		for(E e : b)
			if(a.contains(e))
				r.add(e);
		return r;
	}
	
	public static int[] toArray(ArrayList<Integer> d){
		int[] r = new int[d.size()];
		for(int i =0;i<r.length;i++)
			r[i] = d.get(i);
		return r;
	}
	
	/**
	 * calculates union of supplied segments
	 * 
	 * @param s intervals sorted by starts and stops, each element should have length not less than 2. Only first two values will be used. 
	 * They expected to be start and stop (inclusive), start sould be <= stop. If start >= stp method behavior is unpredictable  
	 * @return
	 */
	public static ArrayList<int[]> union(ArrayList<int[]> s){
		ArrayList<int[]> r = new ArrayList<>();
		if(s.size() == 0)
			return r;
		int start = s.get(0)[0];
		int stop = s.get(0)[1];
		for(int i=1;i<s.size();i++) {
			if(s.get(i)[0]<=stop+1) {
				stop = Math.max(stop,s.get(i)[1]);
			}else {
				r.add(new int[] {start,stop});
				start = s.get(i)[0];
				stop = s.get(i)[1];
			}
		}
		r.add(new int[] {start,stop});
		return r;
	}
	
	public static String intSequenceToShortString(int[] d){
		StringBuilder r = new StringBuilder();
		int start = d[0];
		for(int i=1;i<=d.length;i++){
			if(i<d.length && d[i-1]+1 == d[i])
				continue;
			else if(d[i-1] == start){
				r.append(d[i-1]);
			}else{
				r.append(start).append(":").append(d[i-1]);
			}
			r.append(",");
			if(i<d.length)
				start = d[i];
		}
		return r.substring(0, r.length()-1).toString();
	}
	
	public static void main(String[] args) {
//		System.out.println(Util.join(parseIntSequence("1,2,10"),","));
//		System.out.println(Util.join(parseIntSequence("1:10"),","));
//		System.out.println(Util.join(parseIntSequence("3,2:4,1"),","));
		System.out.println(intSequenceToShortString(new int[]{2,10,1,3,11,18}));
	}

}