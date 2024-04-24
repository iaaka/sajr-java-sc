package util.bio;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import util.Graph;
import util.Interval;
import util.Log;
import util.Settings;
import util.StopIntervalComparator;
import util.SSIHashMap;
import util.Util;
import util.bio.Seg.segPos;
import util.bio.Seg.segType;

public class Gene extends Interval {
	public final String chr_id;
	private ArrayList<Intron> introns;
	private ArrayList<Seg> segs;
	int max_intron_length = 0;
	
	public Gene(int start, int stop, int strand,String chr_id) {
		this(start, stop, strand, chr_id,null);
		
	}

	public Gene(int start, int stop, int strand,String chr_id,String id) {
		super(start, stop, strand,id);
		this.chr_id = chr_id;
		introns = new ArrayList<>();
		segs = new ArrayList<>();
	}
	public void addSeg(Seg s){
		if(s.strand != strand)
			Log.closeWithError("Cannot add segment: segment from another strand.\n"+s+"\n"+this, new RuntimeException());
		if(s.start < start || s.stop > stop)
			Log.closeWithError("Cannot add segment: segment is not within gene:\n"+s+"\n"+this, new RuntimeException());
		segs.add(s);
	}
	
	public void addIntron(Intron i){
		if(i.strand!=0 && i.strand != strand)
			Log.closeWithError("Cannot add intron: intron from another strand.\n"+i+"\n"+this, new RuntimeException());
		if(i.start < start || i.stop > stop)
			Log.closeWithError("Cannot add intron: intron is not within gene:\n"+i+"\n"+this, new RuntimeException());
		introns.add(i);
		max_intron_length = Math.max(max_intron_length, i.length());
	}
	
	public HashSet<Intron> getIntronHash(){
		return new HashSet<Intron>(introns);
	}	
	
	@SuppressWarnings({ "unchecked", "incomplete-switch" })
	public SSIHashMap[] getSegCoverage(){
		SSIHashMap ir = new SSIHashMap();
		SSIHashMap er = new SSIHashMap();
		
		HashMap<Seg,LinkedList<Intron>> seg2intron = linkSeg2Introns();
		
		HashSet<String> barcodes = new HashSet<String>();
		for(Intron j : introns)
			barcodes.addAll(j.getBarcodes());
		for(Seg j : segs)
			barcodes.addAll(j.getBarcodes());
		
		for(String b : barcodes) {
			int[] ecov = new int[segs.size()];
			int tot_first = 0;
			int tot_last = 0;
			for(int i=0;i<segs.size();i++) {
				LinkedList<Intron> ints = seg2intron.get(segs.get(i));
				for(Intron j : ints)
					ecov[i] += j.getCov(b);
				switch(segs.get(i).segpos) {
				case FIRST:
					tot_first += ecov[i];
					break;
				case LAST:
					tot_last += ecov[i];
					break;
				
				}
			}

			for(int i =0;i<segs.size();i++) {
				Seg s = segs.get(i);
				if(s.getId() == null)
					continue;
				switch(s.segpos) {
				case FIRST:
					ir.set(s.getId(), b, ecov[i]);
					er.set(s.getId(), b, (tot_first-ecov[i]));
					break;
				case LAST:
					ir.set(s.getId(), b, ecov[i]);
					er.set(s.getId(), b, (tot_last-ecov[i]));
					break;
				case INTERNAL:
					ir.set(s.getId(), b, s.getCov(b));
					er.set(s.getId(), b, ecov[i]);
					break;
				case ONLY:
					ir.set(s.getId(), b, s.getCov(b));
					er.set(s.getId(), b, ecov[i]);
					break;
				}
			}
		}
		
		return new SSIHashMap[] {ir,er};
	}
	
	
	/**
	 * adds intron and sort introns. used to add new intron while read counting
	 * shouldn't be used frequently
	 * 
	 * probably there are no need to sort, but just for consistency
	 * i'll change it if it will be too slow
	 * @param i
	 */
	public void addIntronAndSort(Intron i) {
		addIntron(i);
		Collections.sort(introns);
	}
	
	/**
	 * 
	 * @param e
	 * @param i
	 * @param l
	 * @return NaN if e and i are zero
	 */
	private Double calcIR(double e, double i, double l) {
		if(e == 0 && i == 0)
			return Double.NaN;
		int rl = Settings.S().getInt(Settings.EFFECTIVE_READ_LENGTH);
		e = e/(rl-1);
		if(Settings.S().getBoolean(Settings.COUNT_ONLY_BORDER_READS) && rl<=l) {
			i = i/(2*rl-2);
		}else {
			i = i/(l+rl-1);
		}
		return i/(i+e);
	}
	
	/**
	 * introns and exons should be sorted
	 * @return for each internal segment gives introns that span it (if any). For first and last segments gives introns that link it to the rest of gene.
	 */
	@SuppressWarnings("incomplete-switch")
	private HashMap<Seg,LinkedList<Intron>> linkSeg2Introns(){
		HashMap<Seg,LinkedList<Intron>> r = new HashMap<>();
		@SuppressWarnings("unchecked")
		ArrayList<Intron> stop_introns = (ArrayList<Intron>)introns.clone();
		StopIntervalComparator comp = new StopIntervalComparator();
		Collections.sort(stop_introns,comp);
		//assumes that gene leftmost segment is first if strand is 1 and last otherwise
		//need to treat unstranded data.
		int rstrand = segs.get(0).segpos == Seg.segPos.FIRST?1:-1;
		for(Seg s : segs) {
			LinkedList<Intron> ints = new LinkedList<>();
			r.put(s, ints);
			//set position by genome rather than by transcript
			Seg.segPos tmp = Seg.segPos.INTERNAL;
			switch(s.segpos){
			case FIRST:
				tmp = rstrand==1?Seg.segPos.FIRST:Seg.segPos.LAST;
				break;
			case LAST:
				tmp = rstrand==1?Seg.segPos.LAST:Seg.segPos.FIRST;
			}
			
			switch(tmp) {
			case FIRST:
				int inx = -Collections.binarySearch(introns, new Interval(s.stop+1,s.stop-1,strand))-1;
				for(;inx<introns.size() && introns.get(inx).start == s.stop+1;inx++)
					ints.add(introns.get(inx));
				break;
			case INTERNAL:
				inx = -Collections.binarySearch(introns, new Interval(s.start,s.start-2,strand))-1;
				for(int i=inx;i<introns.size();i++) {
					if(introns.get(i).start != s.start)
						break;
					ints.add(introns.get(i));
				}
				inx--;
				for(;inx>=0;inx--){
					Intron i= introns.get(inx);
					if(i.start + max_intron_length < s.start)
						break;
					if(s.start <= i.stop)
						ints.add(i);
				}
				break;
			case LAST:
				inx = -Collections.binarySearch(stop_introns, new Interval(s.start+1,s.start-1,strand),comp)-2;
				for(;inx>=0 && introns.get(inx).stop == s.start-1;inx--)
					ints.add(introns.get(inx));
				break;
			
			}
		}
		return r;
	}
	
	public void setSegTypes(){
		Collections.sort(segs);
		Collections.sort(introns);
		int lcnt=0,fcnt=0;
		for(int s =0;s<segs.size();s++){
			if(segs.get(s).segpos == Seg.segPos.FIRST){
				fcnt++;
				continue;
			}
			if(segs.get(s).segpos == Seg.segPos.LAST){
				lcnt++;
				continue;
			}
			if(segs.get(s).segpos == Seg.segPos.ONLY){
				Seg t = segs.get(s);
				segs.set(s, new Seg(t.start,t.stop,t.strand,Seg.segType.EXN,t.segpos,t.getId()));
				continue;
			}
			int inx = Collections.binarySearch(introns, segs.get(s));
			Seg t = segs.get(s);
			if(inx >= 0){		
				segs.set(s, new Seg(t.start,t.stop,t.strand,Seg.segType.INT,t.segpos,t.getId()));
				continue;
			}
			inx = -inx-1;
			for(int i=0;;i++){
				if(inx+i<introns.size()){
					if(introns.get(inx+i).overlap(t)){
						segs.set(s, new Seg(t.start,t.stop,t.strand,Seg.segType.ALT,t.segpos,t.getId()));
						break;
					}
				}
				if(inx-i-1>=0){
					if(introns.get(inx-i-1).overlap(t)){
						segs.set(s, new Seg(t.start,t.stop,t.strand,Seg.segType.ALT,t.segpos,t.getId()));
						break;
					}
				}else if(inx+i>=introns.size()){
					segs.set(s, new Seg(t.start,t.stop,t.strand,Seg.segType.EXN,t.segpos,t.getId()));
					break;
				}
			}
		}
		//set first and last
		for(int s =0;s<segs.size();s++){
			Seg t = segs.get(s);
			if(t.segpos == Seg.segPos.FIRST){
				segs.set(s, new Seg(t.start,t.stop,t.strand,fcnt>1?Seg.segType.ALT:Seg.segType.EXN,t.segpos,t.getId()));
				continue;
			}
			if(t.segpos == Seg.segPos.LAST){
				segs.set(s, new Seg(t.start,t.stop,t.strand,lcnt>1?Seg.segType.ALT:Seg.segType.EXN,t.segpos,t.getId()));
				continue;
			}
		}
		
	}
	
	
	public int getIntronCount(){
		return introns.size();
	}
	
	public int getSegCount(){
		return segs.size();
	}
	public Seg getSeg(int i){
		return segs.get(i);
	}
	
	public Intron getIntron(int i){
		return introns.get(i);
	}
	
	/**
	 * prepare gene for read counting:
	 * sorts segments and introns
	 * adds pseudo introns
	 * checks for ONLY segments
	 * check whether all splice sites (from introns) are borders of exons
	 */
	public void prepare(){
		Collections.sort(segs);
		Collections.sort(introns);
		//fill not-retained introns with pseudo-segments
		boolean has_only = segs.get(0).segpos == Seg.segPos.ONLY;
		ArrayList<Seg> t = new ArrayList<>();
		for(int i=1;i<segs.size();i++) {
			has_only |= segs.get(i).segpos == Seg.segPos.ONLY;
			if(segs.get(i-1).stop+1 != segs.get(i).start) {
				if(segs.get(i).start < segs.get(i-1).stop+1)
					Log.closeWithError("Segments "+segs.get(i-1).getId()+" and "+segs.get(i).getId()+" from gene " +getId()+" overlap. Segments are regions between nearest splice sites, they cannot overlap by defenition!", new RuntimeException());
				t.add(new Seg(segs.get(i-1).stop+1,segs.get(i).start-1,strand,Seg.segType.INT,Seg.segPos.INTERNAL));
			}
		}
		if(has_only && segs.size() != 1)
			Log.warn("Gene "+getId()+" is multiexon, but have ONLY segment.");
		if(t.size() != 0) {
			segs.addAll(t);
			Collections.sort(segs);
		}
		//check whether all splice sites (from introns) are borders of exons
		@SuppressWarnings("unchecked")
		ArrayList<Seg> segs_end = (ArrayList<Seg>) segs.clone();
		StopIntervalComparator comp = new StopIntervalComparator();
		Collections.sort(segs_end,comp);
		for(Intron i : introns) {
			int start_inx = Collections.binarySearch(segs_end,new Interval(i.start-1,i.start-1,i.strand),comp);
			int stop_inx = Collections.binarySearch(segs,new Interval(i.stop+1,i.stop+1,i.strand));
			if(start_inx<0) start_inx = -start_inx-2;
			if(stop_inx<0) stop_inx = -stop_inx-1;
				
			if(start_inx == -1 || stop_inx == segs.size() || segs_end.get(start_inx).stop != i.start-1 ||
			   segs.get(stop_inx).start != i.stop+1) {
				Log.closeWithError("Intron '"+i+"' is not from gene "+getId()+": there are no adjusted segments.", new RuntimeException());
			}
		}
	}
		
	public HashSet<Seg> getSegForRead(int[] r){
		HashSet<Seg> res = new HashSet<>();
		if(Settings.S().getBoolean(Settings.COUNT_ONLY_BORDER_READS) & r.length==2) {
			return res;
		}
		for(int i=0;i<r.length;i+=2) {
			int inx = Collections.binarySearch(segs, new Interval(r[i],r[i+1],strand)); 
			if(inx < 0)
				inx = -inx-1;
			if(inx>0 && segs.get(inx-1).stop>= r[i])
				res.add(segs.get(inx-1));
			for(;inx<segs.size();inx++){
				if(segs.get(inx).start > r[i+1])
					break;
				res.add(segs.get(inx));
			}
		}
		// I will require that at least one segment border should be splicing site used by read
		// something unexpected can happen for site bordered not by splicing sites....
		if(Settings.S().getBoolean(Settings.COUNT_ONLY_BORDER_READS)){
			// prepare junction borders
			HashSet<Integer> starts = new HashSet<>();
			HashSet<Integer> stops = new HashSet<>();
			for(int i=1;i<r.length-1;i+=2) {
				stops.add(r[i]);
				starts.add(r[i+1]);	
			}
			HashSet<Seg> res_ = new HashSet<>();
			for(Seg s : res) {
				if(starts.contains(s.start) || stops.contains(s.stop))
					res_.add(s);
			}
			res = res_;
		}
		return res;
	}
	
	/**
	 * transform exon representation (exons should be added to gene through addSeg method) to seg representation
	 * skips alternative TSS and polyA if they are in exon.
	 * Changes objects, but coordinates cannot be changes because they are final. 
	 * So, gene with corrected coordinates is returned. 
	 * 
	 * @return new gene with corrected coordinates (it is likely that original coordinates where wrong).
	 * Intron and seg list are the same. It is save to use only new copy of the gene.
	 */
	public Gene exon2seg() {
		if(segs.size()==0)
			Log.closeWithError("Wrong annotation file format: Gene "+getId()+" has no exons.",null);
		//make splice sites
		HashSet<Intron> ints = new HashSet<>();
		HashSet<Integer> sites_ = new HashSet<>();
		for(Intron i : introns) {
			if(i.length() == 0) //skipp zero length introns
				continue;
			sites_.add(i.start-1);
			sites_.add(i.stop);
			ints.add(i);
		}
		//remove duplicate introns if any
		introns = new ArrayList<>(ints);
		ArrayList<Integer> sites = new ArrayList<>(sites_);
		Collections.sort(sites);
		//make exonic regions
		Collections.sort(segs);
		ArrayList<int[]> regs = new ArrayList<>(segs.size());
		for(Seg s : segs) {
			regs.add(new int[] {s.start,s.stop});
		}
		regs = Util.union(regs);
		segs = new ArrayList<>();
		int sinx = 0;
		int rinx = 0;
		int start = regs.get(0)[0];
		Seg.segPos startpos = Seg.segPos.FIRST;
		//make segments
		for(;;) {
			//site within exonic region
			if(sinx != sites.size() && sites.get(sinx)<regs.get(rinx)[1]) {
				if(strand == -1 && startpos == Seg.segPos.FIRST)
					startpos = Seg.segPos.LAST;
				segs.add(new Seg(start,sites.get(sinx),strand,segType.NA,startpos));
				start = sites.get(sinx)+1;
				startpos = Seg.segPos.INTERNAL;
				sinx++;
				continue;
			}
			if(sinx != sites.size() && sites.get(sinx) == regs.get(rinx)[1]) {//constitutive intron
				if(strand == -1 && startpos == Seg.segPos.FIRST)
					startpos = Seg.segPos.LAST;
				segs.add(new Seg(start,sites.get(sinx),strand,segType.NA,startpos));
				sinx++;
			}else {//last exon
				if(startpos == Seg.segPos.INTERNAL)
					if(strand == -1)
						startpos = Seg.segPos.FIRST;
					else
						startpos = Seg.segPos.LAST;
				else
					startpos = Seg.segPos.ONLY;
				segs.add(new Seg(start,regs.get(rinx)[1],strand,segType.NA,startpos));
			}
			rinx++;
			if(sinx == sites.size())
				break;
			if(rinx == regs.size())
				break;
			start = regs.get(rinx)[0];
			if(start == sites.get(sinx)+1) {
				startpos = Seg.segPos.INTERNAL;
				sinx++;
			}else
				startpos = Seg.segPos.FIRST;
		}
		for(;rinx<regs.size();rinx++) {
			segs.add(new Seg(regs.get(rinx)[0],regs.get(rinx)[1],strand,segType.NA,segPos.ONLY));
		}
		setSegTypes();
		
		Gene g = new Gene(regs.get(0)[0], regs.get(regs.size()-1)[1], strand, chr_id,getId());
		g.introns = introns;
		g.segs = segs;
		return g;
	}
	
	public void printTranscripts(PrintStream out,HashMap<String,Double> seg2psi,double irPSI, double asPSI){
		GFFeature f = new GFFeature(chr_id, "gene",start ,stop , strand,"SAJR",0,0);
		f.addAttr("gene_id", getId());
		out.println(f);
		Graph<Integer, Interval> g = makeSpliceGraph(seg2psi,irPSI,asPSI);
		ArrayList<ArrayList<Interval>> transcs = g.getPaths(start-1, stop);
		int tr_id = 1;
		for(ArrayList<Interval> t : transcs){
			if(t.get(0) instanceof Intron || t.get(t.size()-1) instanceof Intron)
				Log.closeWithError("gene '"+getId()+"' starts or ends with intron. It isn't possible!", new RuntimeException());
			//save transcript
			f = new GFFeature(chr_id, "transcript",t.get(0).start ,t.get(t.size()-1).stop , strand,"SAJR",0,0);
			f.addAttr("gene_id", getId());
			f.addAttr("transcript_id", getId()+".t"+tr_id);
			out.println(f);
			int exn_no = 1;
			int start = t.get(0).start;
			int stop = t.get(0).start-1;
			for(int i = 0;i<=t.size();i++){
				if(i != t.size() && t.get(i) instanceof Intron) continue;
				if(i == t.size() || stop != t.get(i).start-1){
					f = new GFFeature(chr_id, "exon",start ,stop , strand,"SAJR",0,0);
					f.addAttr("gene_id", getId());
					f.addAttr("transcript_id", getId()+".t"+tr_id);
					f.addAttr("exon_number", ""+exn_no);
					out.println(f);
					
					if(i < t.size()){
						start = t.get(i).start;
						stop = t.get(i).stop;
						exn_no++;
					}
				}else
					stop = t.get(i).stop;
			}
			tr_id++;
		}
	}
	
	/**
	 * 
	 * @param seg2psi use null to include all segments
	 * @param irPSI
	 * @param asPSI
	 * @return
	 */
	public Graph<Integer, Interval> makeSpliceGraph(HashMap<String,Double> seg2psi,double irPSI, double asPSI){
		Graph<Integer, Interval> r = new Graph<>();
		for(Intron i : introns){
			r.add(i.start-1);
			r.add(i.stop);
			r.addDirectedEdge(i.start-1,i.stop, i);
		}
		// beginning of the gene (by chr coordinates) is gene start, end is gene.stop
		r.add(start-1);
		r.add(stop);
		for(Seg s : segs){
			if(s.getId() == null) continue;			
			if(seg2psi != null){
				Double psi = seg2psi.get(s.getId());
				if(psi == null)
					Log.closeWithError("Segment '"+s+"' absents in read count file.", new RuntimeException());
				if(!psi.isNaN() && s.segtype.equals(Seg.segType.ALT) && psi < asPSI) continue;
				if(!psi.isNaN() && s.segtype.equals(Seg.segType.INT) && psi < irPSI) continue;
			}
			switch(s.segpos){
			case ONLY:
				r.addDirectedEdge(start-1, stop, s);
				break;
			case FIRST:
				if(s.strand==1){
					r.add(s.stop);
					r.addDirectedEdge(start-1, s.stop, s);
				}else{
					r.add(s.start-1);
					r.addDirectedEdge(s.start-1, stop, s);
				}break;				
			case INTERNAL:
				r.add(s.start-1);
				r.add(s.stop);
				r.addDirectedEdge(s.start-1, s.stop, s);
				break;
			case LAST:
				if(s.strand==1){
					r.add(s.start-1);
					r.addDirectedEdge(s.start-1, stop, s);
				}else{		
					r.add(s.stop);
					r.addDirectedEdge(start-1, s.stop, s);
				}break;
			}
		}
		return r;
	}
}
