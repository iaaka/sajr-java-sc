package util.bio;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import util.Log;
import util.Settings;
import util.Util;

public class ChrAnnotation {
	int max_gene_length=0;
	ArrayList<Gene> genes;
	HashMap<Intron,List<Gene>> intron2genes;
	HashMap<Integer,HashSet<Gene>> leftSS2genes;
	HashMap<Integer,HashSet<Gene>> rightSS2genes;
	HashMap<Intron,Intron> introns;
	String chr_id;

	public ChrAnnotation(String chr_id) {
		genes = new ArrayList<>();
		this.chr_id = chr_id;
	}
	
	@SuppressWarnings("unchecked")
	public  ArrayList<Gene> getGenes(){
		return (ArrayList<Gene>)genes.clone();
	}
	
	public  ArrayList<Intron> getIntrons(){
		return new ArrayList<>(introns.keySet());
	}

	/**
	 * if two genes share the same introns, these introns should be the same object.
	 * @param g
	 */
	public void addGene(Gene g){
		genes.add(g);
		max_gene_length = Math.max(max_gene_length, g.length());
	}
	
	/**
	 * the method should be called when all genes from the chromosome were loaded
	 * sorts genes, introns and segments within genes and makes introns to gene hash. 
	 */
	public void loaded(){
		Collections.sort(genes);
		intron2genes = new HashMap<>();
		introns = new HashMap<>();
		if(Settings.S().getBoolean(Settings.LOOK_FOR_GENE_FOR_UNKNOWN_JUNCTIONS)) {
			leftSS2genes = new HashMap<>();
			rightSS2genes = new HashMap<>();
		}
		for(Gene g : genes){
			g.prepare();
			for(int i=0;i<g.getIntronCount();i++){
				Intron in = g.getIntron(i);
				introns.put(in, in);
				List<Gene> gs = intron2genes.get(in);
				if(gs == null){
					gs = new LinkedList<>();
					intron2genes.put(in, gs);
				}
				gs.add(g);
				if(Settings.S().getBoolean(Settings.LOOK_FOR_GENE_FOR_UNKNOWN_JUNCTIONS)) {
					HashSet<Gene> sgs = leftSS2genes.get(in.start);
					if(sgs == null) {
						sgs = new HashSet<>();
						leftSS2genes.put(in.start, sgs);
					}
					sgs.add(g);
					sgs = rightSS2genes.get(in.stop);
					if(sgs == null) {
						sgs = new HashSet<>();
						rightSS2genes.put(in.stop, sgs);
					}
					sgs.add(g);
				}
			}
		}
	}
	
	/**
	 * 
	 * @param r
	 * @param strand
	 * @return set of introns that are with this read
	 */
	private HashSet<Intron> getIntronsForRead(int[] r,int strand){
		HashSet<Intron> res = new HashSet<>();
		for(int i=2;i<r.length;i+=2) {
			Intron newInt = new Intron(r[i-1]+1,r[i]-1,strand);
			Intron in = introns.get(newInt);
			if(in != null)
				res.add(in);
			else {
				Log.addStat(Log.NEW_JUNCTIONS_FOUND, 1);
				res.add(newInt);
				introns.put(newInt, newInt);
				LinkedList<Gene> int2genes = new LinkedList<>();
				intron2genes.put(newInt, int2genes);
				if(Settings.S().getBoolean(Settings.LOOK_FOR_GENE_FOR_UNKNOWN_JUNCTIONS)){
					HashSet<Gene> lg = leftSS2genes.get(r[i-1]+1);
					HashSet<Gene> rg = rightSS2genes.get(r[i]-1);
					HashSet<Gene> gs = new HashSet<>();
					if(lg != null && rg != null) {
						gs = Util.intersect(lg, rg);
					}else { // search for gene by overlap
						gs = getGenesByOverlap(newInt.start,newInt.stop,strand);	
					}
					for(Gene g: gs) {
						if((strand==0 || g.strand == strand) && g.start <= newInt.start && g.stop>=newInt.stop){
							g.addIntronAndSort(newInt);
							int2genes.add(g);
						}
					}
				}
			}
		}
		return res;
	}
	
	/**
	 * 
	 * @param ints
	 * @return set of genes that contains ALL/'at least one' (depends on ONLY_JUNCTIONS_FROM_SAME_GENE) introns
	 */
	private HashSet<Gene> getGenesByIntrons(HashSet<Intron> ints){
		if(ints.size() == 0)
			return new HashSet<>();
		HashSet<Gene> g = new HashSet<Gene>();
		int j = 0;
		for(Intron i : ints) {
			List<Gene> gns = intron2genes.get(i);
			if(gns != null){
				if(j>0 && Settings.S().getBoolean(Settings.ONLY_JUNCTIONS_FROM_SAME_GENE))
					g = Util.intersect(g, gns);
				else
					g.addAll(gns);
			}
			j++;
		}
		return g;
	}
	
	/**
	 * 
	 * @param ints
	 * @return set of genes that contains at least one intron
	 */
	protected HashSet<Gene> getAllGenesByIntrons(HashSet<Intron> ints){
		HashSet<Gene> g = new HashSet<Gene>();
		for(Intron i : ints){
			if(intron2genes.containsKey(i)){
				g.addAll(intron2genes.get(i));
			}
		}
		return g;
	}
	
	protected HashSet<Gene> getGenesByOverlap(int start, int stop, int strand){
		if(strand == 0){
			HashSet<Gene> r = getGenesByOverlapStranded(start,stop,1);
			r.addAll(getGenesByOverlapStranded(start,stop,-1));
			return r;
		}else{
			return getGenesByOverlapStranded(start,stop,strand);
		}
	}
	/**
	 * @param strand shouldn't be 0
	 * @return returns genes that overlap with given region
	 */
	private HashSet<Gene> getGenesByOverlapStranded(int start, int stop, int strand){
		if(strand != 1 && strand != -1)
			throw new RuntimeException("strand should be either 1 or -1");
		HashSet<Gene> gs = new HashSet<>();
		int inx = Collections.binarySearch(genes, new Gene(start,stop,strand,chr_id));
		if(inx < 0)
			inx = -inx-1;
		//forward
		for(int i=inx;i<genes.size();i++){
			Gene g = genes.get(i);
			if(g.start > stop || g.strand != strand)
				break;
			gs.add(g);
		}
		//backward
		for(int i=inx-1;i>=0;i--){
			Gene g = genes.get(i);
			if(g.start + max_gene_length < start || g.strand != strand)
				break;
			if(start <= g.stop)
				gs.add(g);
		}
		return gs;
	}
	
	
	private void addCov2Genes(HashMap<Gene,HashSet<Seg>> gene2segs,HashSet<Intron> cintrons,boolean paired,String bc)	{
		if(gene2segs.size() > 0)
			Log.addStat(Log.GENE_READS, paired?2:1);
		//filter genes were reads overlap only intron (if other exists)
		HashMap<Gene,HashSet<Seg>> gene2segs_ = new HashMap<>();
		for(Gene g : gene2segs.keySet()) {
			for(Seg s : gene2segs.get(g))
				if(s.segtype != Seg.segType.INT) {
					gene2segs_.put(g, gene2segs.get(g));
					break;
				}
		}
		if(gene2segs_.size()>0) 
			gene2segs = gene2segs_;
		//count read
		boolean count_junc = false;
		for(Gene g : gene2segs.keySet()) {
			boolean has_intron = false;
			boolean has_exn = false;
			boolean has_internal_exn = false;
			for(Seg s : gene2segs.get(g)) {
				has_intron = has_intron || s.segtype == Seg.segType.INT;
				has_exn = has_exn || s.segtype == Seg.segType.EXN;
				has_internal_exn = has_internal_exn || (s.segtype == Seg.segType.EXN && (s.segpos == Seg.segPos.INTERNAL || s.segpos == Seg.segPos.ONLY));
			}
			if(!Settings.S().getBoolean(Settings.COUNT_ONLY_BORDER_READS) || gene2segs.get(g).size()!=1){
				for(Seg s : gene2segs.get(g)) {
					if(s.segtype == Seg.segType.INT)
						s.addCov(bc);
					else if(!has_intron || Settings.S().getBoolean(Settings.COUNT_INTRON_READS)) {
						s.addCov(bc);
						count_junc = true;
					}
				}
			}
			if((!has_intron || Settings.S().getBoolean(Settings.COUNT_INTRON_READS)) && has_exn &&
				(has_internal_exn || !Settings.S().getBoolean(Settings.COUNT_ONLY_INTERNAL))) {
					g.addCov(bc);
			}				
		}
		 
		for(Intron i : cintrons)
			if(count_junc || intron2genes.get(i).size()==0) // if it is unknown junction then we can count any read for it
				i.addCov(bc);
			
		if(count_junc)
			Log.addStat(Log.EXON_READS, paired?2:1);
	}
	
	private boolean allIntronsHaveGenes(HashSet<Intron> ints){
		for(Intron in : ints)
			if(intron2genes.get(in).size() == 0)
				return false;
		return true;
	}
	
	public void addRead(int[] r,int strand,String bc) {
		Log.addStat(Log.USED_READS, 1);
		HashSet<Intron> cintrons = new HashSet<>();
		HashSet<Gene> cgenes = null;
		if(r.length > 2) {
			Log.addStat(Log.JUNCTIONS_CNT, 1);
			cintrons = getIntronsForRead(r, strand);
			if(!allIntronsHaveGenes(cintrons)) {
				Log.addStat(Log.UNKNOWN_JUNCTION, 1);
				if(!Settings.S().getBoolean(Settings.USE_READS_WITH_UNKNOWN_JUNCTIONS))
					return;
			}
			cgenes = getGenesByIntrons(cintrons);
			if(cgenes.size()==0 && cintrons.size()!=0)  
				Log.addStat(Log.UNKNOWN_JUNCTION_COMB, 1);
		}else
			cgenes = getGenesByOverlap(r[0], r[1], strand);
		
		HashMap<Gene,HashSet<Seg>> gene2segs = new HashMap<>();
		for(Gene g : cgenes) {
			gene2segs.put(g,g.getSegForRead(r));
		}
		addCov2Genes(gene2segs,cintrons,false,bc);
	}
	
	
	public void addReads(int[] r1,int[] r2, int strand,String bc) {
		Log.addStat(Log.USED_READS, 2);
		HashSet<Intron> cintrons = new HashSet<>();
		HashSet<Gene> cgenes = null;
		if(r1.length > 2 || r2.length > 2) {
			Log.addStat(Log.JUNCTIONS_CNT, 1);
			cintrons = getIntronsForRead(r1, strand);
			HashSet<Intron> tmp = getIntronsForRead(r2, strand);
			boolean unknown_junc = !allIntronsHaveGenes(cintrons) || !allIntronsHaveGenes(tmp);
			if(unknown_junc)
				Log.addStat(Log.UNKNOWN_JUNCTION, 2);
			if(unknown_junc && !Settings.S().getBoolean(Settings.USE_READS_WITH_UNKNOWN_JUNCTIONS)) 
				return;
			cintrons.addAll(tmp);
			cgenes = getGenesByIntrons(cintrons);
			if(cgenes.size()==0 && cintrons.size() != 0)
				Log.addStat(Log.UNKNOWN_JUNCTION_COMB, 2);
		}else {
			cgenes = getGenesByOverlap(r1[0], r1[1], strand);
			cgenes.addAll(getGenesByOverlap(r2[0], r2[1], strand));
		}
		HashMap<Gene,HashSet<Seg>> gene2segs = new HashMap<>();
		for(Gene g : cgenes) {
			HashSet<Seg> segs = g.getSegForRead(r1);
			segs.addAll(g.getSegForRead(r2));
			gene2segs.put(g,segs);
		}
		addCov2Genes(gene2segs,cintrons,true,bc);
	}
	
	public String getID() {
		return chr_id;
	}
}
