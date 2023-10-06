package rc;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;

import run.Run;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import util.Log;
import util.Settings;
import util.SSIHashMap;
import util.Util;
import util.bio.ChrAnnotation;
import util.bio.GFFException;
import util.bio.GFFeature;
import util.bio.GFFParser;
import util.bio.Gene;
import util.bio.Intron;
import util.bio.Seg;

class SingleReadReader {
	HashMap<String, ChrAnnotation> chrs;
	
	public SingleReadReader(HashMap<String, ChrAnnotation> chrs) {
		this.chrs = chrs;
	}

	public void read(SAMRecord r) {
		ChrAnnotation c = chrs.get(r.getReferenceName());
		if(c==null){
			c = new ChrAnnotation(r.getReferenceName());
			c.loaded();
			chrs.put(r.getReferenceName(),c);
		}
		int[] inters = Util.getMapIntervals(r);
		int strand = Settings.S().getInt(Settings.STRANDED)*(r.getReadNegativeStrandFlag()?-1:1);
		if(r.getReadPairedFlag() && !r.getFirstOfPairFlag())
			strand = -strand;
		if(r.getCharacterAttribute("XS") != null) //prabably it whould be better to check that it isn't contradict to getReadNegativeStrandFlag...
			strand = r.getCharacterAttribute("XS") == '+'?1:-1;
		String bc = r.getStringAttribute("CB");
		c.addRead(inters, strand,bc);
	}
}

class PairedReadReader {
	HashMap<String, ChrAnnotation> chrs;
	String read_name = null;
	//list of records by: read name;first_mate;chr;start;mate_start
	HashMap<String,LinkedList<SAMRecord>> records = new HashMap<>(); 
	
	public PairedReadReader(HashMap<String, ChrAnnotation> chrs) {
		this.chrs = chrs;
	}

	public void read(SAMRecord r) {
		String n = r.getReadName();
		//truncate read mate information if exists
		if((n.charAt(n.length()-1) == '1' || n.charAt(n.length()-1) == '2') &&
		   (n.charAt(n.length()-2) == ';' || n.charAt(n.length()-2) == '\\' || n.charAt(n.length()-2) == '/'))
			n = n.substring(0,n.length()-2);
		String mate_key = n+";"+!r.getFirstOfPairFlag()+";"+r.getMateReferenceName()+";"+r.getMateAlignmentStart()+";"+r.getAlignmentStart();
		LinkedList<SAMRecord> mates = records.get(mate_key);
		if(mates == null) {
			String its_key = n+";"+r.getFirstOfPairFlag()+";"+r.getReferenceName()+";"+r.getAlignmentStart()+";"+r.getMateAlignmentStart();
			mates = records.get(its_key);
			if(mates == null) {
				mates = new LinkedList<>();
				records.put(its_key,mates);
			}
			mates.add(r);
		}else {
			SAMRecord f = r.getFirstOfPairFlag()?r:mates.removeFirst();
			SAMRecord s = r.getFirstOfPairFlag()?mates.removeFirst():r;
			if(mates.size() == 0)
				records.remove(mate_key);
			ChrAnnotation c = chrs.get(f.getReferenceName());
			if(c==null){
				c = new ChrAnnotation(f.getReferenceName());
				c.loaded();
				chrs.put(f.getReferenceName(),c);
			}
			int[] r1 = Util.getMapIntervals(f);
			int[] r2 = Util.getMapIntervals(s);
			int strand = Settings.S().getInt(Settings.STRANDED)*(f.getReadNegativeStrandFlag()?-1:1);
			if(f.getCharacterAttribute("XS") != null) //Probably it would be better to check that it isn't contradict to getReadNegativeStrandFlag...
				strand = f.getCharacterAttribute("XS") == '+'?1:-1;
			if(s.getCharacterAttribute("XS") != null) //Probably it would be better to check that it isn't contradict to getReadNegativeStrandFlag...
				strand = s.getCharacterAttribute("XS") == '+'?1:-1;
			String bc = r.getStringAttribute("CB");
			c.addReads(r1, r2, strand,bc);

		}
	}
	
	public void finish() {
		if(records.size()!=0) {
			int i = 0;
			for(LinkedList<SAMRecord> t : records.values())
				i += t.size();
			SAMRecord t = records.values().iterator().next().getFirst();
			Log.throwUncrucialExc("There are "+i+" read locations that doesn't have expected mate records. " +
					"For example read "+t.getReadName()+", location "+t.getReferenceName()+":"+t.getAlignmentStart()+" "+
					"should have mate mapped to "+t.getMateReferenceName()+":"+t.getMateAlignmentStart()+" "+
					"but it cannot be found. It could happen if one mate was filtered out. " +
					"For eaxample if NH attributes is one (or absent) in one of mates, " +
					"while in other it is more than 1 and use_mult is set to false.");
		}
	}
}

public class ReadCounter {
	HashMap<String, ChrAnnotation> chrs;
	ArrayList<Gene> genes;
	// to collapse UMIs
	int pos = -1;
	HashSet<String> reads = new HashSet<String>();

	public ReadCounter() {
		try {
			loadGff(Settings.S().getString(Settings.ANN_IN));
		} catch (Exception e) {
			Log.closeWithError("Cannot read annotation file: "+Settings.S().getString(Settings.ANN_IN),e);
		}
	}
	
	/**
	 * if data is unstranded, loaded annotation will be also unstranded
	 * @param f
	 * @throws IOException
	 * @throws GFFException 
	 */
	private void loadGff(String fname) throws IOException, GFFException{
		String[] black_list = Settings.S().getString(Settings.GENE_BLACK_LIST).split("@");
		Arrays.sort(black_list);
		chrs = new HashMap<>();
		genes = new ArrayList<>();
		int[] stat = new int[3];
		Gene g = null;
		ChrAnnotation chr = null;
		HashMap<Intron,Intron> ints = null;
		//read
		GFFParser gffp = new GFFParser(fname);
		ArrayList<GFFeature> gff = new ArrayList<GFFeature>();
		for(GFFeature f=gffp.next();f != null;f=gffp.next())
			gff.add(f);
		gffp.close();
		//sort
		Collections.sort(gff,new Comparator<GFFeature>() {
			public int compare(GFFeature o1, GFFeature o2) {
				if(!o1.seqname.equals(o2.seqname))
					return o1.seqname.compareTo(o2.seqname);
				String gid1 = o1.getAttr("gene_id");
				String gid2 = o2.getAttr("gene_id");
				if(!gid1.equals(gid2))
					return gid1.compareTo(gid2);
				//genes first
				if(o1.feature.equals("gene") && !o2.feature.equals("gene"))
					return -1;
				if(!o1.feature.equals("gene") && o2.feature.equals("gene"))
					return 1;
				return 0;
				 
			}
		});
		//parse
		for(GFFeature f : gff){
			if(f.getAttr("gene_id") == null)
				Log.closeWithError("Wrong annotation file format. Features should have gene_id attribute.", null);
			if(Arrays.binarySearch(black_list, f.getAttr("gene_id")) >= 0)
				continue;
			int strand = f.strand;
			switch(f.feature) {
			case "gene":
				stat[0]++;
				if(chr == null || !chr.getID().equals(f.seqname)) {
					if(chr != null)
						chr.loaded();
					chr = new ChrAnnotation(f.seqname);
					if(chrs.containsKey(f.seqname))
						Log.closeWithError("Annotation file isn't sorted by chr: '"+f.seqname+"' meet at least twice!",new RuntimeException());
					chrs.put(f.seqname, chr);
					ints = new HashMap<Intron, Intron>();
				}
				g = new Gene(f.start, f.stop,strand, f.seqname,f.getAttr("gene_id"));
				chr.addGene(g);
				genes.add(g);
				break;
			case "segment":
				try {
					g.addSeg(new Seg(f.start, f.stop,strand, Seg.segType.valueOf(f.getAttr("type")), Seg.segPos.valueOf(f.getAttr("position")),f.getAttr("segment_id")));
				}catch(IllegalArgumentException e) {
					Log.closeWithError("Unknown segment type = '"+f.getAttr("type")+
							"' or position = '"+f.getAttr("position")+"'. " +
							"Type should be in: "+Util.join(Seg.segType.values(), ", ")+
							", position should be in: "+Util.join(Seg.segPos.values(), ", "), e);
				}
				stat[1]++;
				break;
			case "intron":
				Intron i = new Intron(f.start, f.stop,strand);
				if(!ints.containsKey(i))
					ints.put(i, i);
				g.addIntron(ints.get(i));
				stat[2]++;
				break;
			default:
				Log.closeWithError("Annotation contains unknown feature: '"+f.feature+"'. Only gene,segment and intron are allowed.", new RuntimeException());
			}
			
		}
		if(chr!= null)
			chr.loaded();
		Log.println("Annotation loaded: #chr="+chrs.size()+"; #genes="+stat[0]+"; #segs="+stat[1]+"; #introns="+stat[2]);
	}
	
	private void countReads() throws IOException  {
		if(!(new File(Settings.S().getString(Settings.IN))).exists())
			Log.closeWithError("Input file '"+Settings.S().getString(Settings.IN)+"' doesn't exists",new RuntimeException());
		SamReader in = SamReaderFactory.makeDefault().open(SamInputResource.of(new BufferedInputStream(new FileInputStream(Settings.S().getString(Settings.IN)),10000000)));
		if(!in.getFileHeader().getAttribute("SO").equals("coordinate"))
			Log.closeWithError("Input file '"+Settings.S().getString(Settings.IN)+"' is not sorted",new RuntimeException());
		SingleReadReader sreader = new SingleReadReader(chrs);
		PairedReadReader preader = new PairedReadReader(chrs);
		int i = 0;
		Iterator<SAMRecord> samIterator = in.iterator();
		for(;samIterator.hasNext();) {
			try{
				SAMRecord r = samIterator.next();
				Log.addStat(Log.TOTAL_READS, 1);
				if(i % 10000000 == 0) {
					Log.println(i+" lines parsed");
				}
				i++;
				if(!accept(r)) {
					continue;
				}
				if(Settings.S().getInt(Settings.PAIRED)==0 || !r.getReadPairedFlag() || !r.getProperPairFlag())
					sreader.read(r);
				else
					preader.read(r);
			}catch(SAMFormatException e){
				Log.throwUncrucialExc("Something wrong with read SAM format: "+e.getMessage()+".\n"
						+ "The read was skipped.");
			}
		}
		Log.println(i+" lines parsed");
		preader.finish();
		in.close();
	}
	
	public static void countAndPrint() throws IOException {
		Log.addStat(Log.UNMAPPED, 0);
		Log.addStat(Log.MULTI_READS, 0);
		Log.addStat(Log.EXON_READS, 0);
		Log.addStat(Log.GENE_READS, 0);
		Log.addStat(Log.JUNCTIONS_CNT, 0);
		Log.addStat(Log.PAIRED, 0);
		Log.addStat(Log.SINGLETONS, 0);
		Log.addStat(Log.TOTAL_READS, 0);
		Log.addStat(Log.UNKNOWN_JUNCTION, 0);
		Log.addStat(Log.UNKNOWN_JUNCTION_COMB, 0);
		Log.addStat(Log.USED_READS, 0);
		Log.addStat(Log.PCR_DUPLICATES, 0);
		if(Settings.S().getBoolean(Settings.LOOK_FOR_GENE_FOR_UNKNOWN_JUNCTIONS))
			Log.addStat(Log.NEW_JUNCTIONS_FOUND, 0);
		Log.println("Count reads: "+Settings.S().getString(Settings.IN)+" -> "+Settings.S().getString(Settings.OUT_BASE));
		ReadCounter r = new ReadCounter();
		r.countReads();
		try {
//			r.printGeneCov(r.genes);
			r.printSegCov(r.genes);
//			r.printIntronCov();
			Log.printStat();
			Log.cleanStat();
		}catch(IOException e) {
			Log.closeWithError("Cannot write output: "+e.getMessage(),e);
		}
		
	}
	
	public boolean accept(SAMRecord r) {
		if(Settings.S().getInt(Settings.PAIRED)!=0 && r.getReadPairedFlag() && r.getProperPairFlag() && !r.getReferenceIndex().equals(r.getMateReferenceIndex())) { 
			Log.warn("Mates of read '"+r.getReadName()+"' are from different chromosomes, while bam FLAG says that they are properly paired. They will be treated as singletons.");
			r.setProperPairFlag(false);
		}
		
		if(r.getReadUnmappedFlag()) {
			Log.addStat(Log.UNMAPPED, 1);
			return false;
		}
		Integer nh = (Integer)r.getAttribute("NH");
		if(nh == null) {
			Log.throwUncrucialExc("Read "+r.getReadName()+" doesn't have NH attribute!");
			nh = 1;
		}
		if(nh > 1)
			Log.addStat(Log.MULTI_READS, 1);
		if(r.getReadPairedFlag() && r.getProperPairFlag())
			Log.addStat(Log.PAIRED, 1);
		else
			Log.addStat(Log.SINGLETONS, 1);	
		
//	if(r.getDuplicateReadFlag()) {
//		Log.addStat(Log.PCR_DUPLICATES, 1);
//		return false;
//	}
	// starsolo doesn't set duplicate read flag, so I'll use chr-start-cigar-cb-ub to collapse reads to umi
	if(r.getStringAttribute("CB") != null && r.getStringAttribute("UB") != null) {
		int start = r.getAlignmentStart();
		if(start != pos) {
			pos = start;
			reads.clear();
		}
		if(!reads.add(r.getContig()+"|"+start+"|"+r.getCigarString()+"|"+r.getStringAttribute("CB")+"|"+r.getStringAttribute("UB"))) {
			Log.addStat(Log.PCR_DUPLICATES, 1);
			return false;
		}
	}
	//
		
		return (Settings.S().getBoolean(Settings.USE_MULT) || nh == 1) 
				&& (Settings.S().getInt(Settings.PAIRED)==0 || Settings.S().getBoolean(Settings.USE_SINGLETONS) || (r.getReadPairedFlag() && r.getProperPairFlag()));
	}
	
	
	private void printSegCov(ArrayList<Gene> genes) throws FileNotFoundException {
		SSIHashMap ir = new SSIHashMap();
		SSIHashMap er = new SSIHashMap();
		
		for(Gene g : genes) {
			SSIHashMap[] t = g.getSegCoverage();
			ir.addAll(t[0]);
			er.addAll(t[1]);
		}
		HashSet<String> t = ir.getKeys1();
		t.addAll(er.getKeys1());
		String[] segids = t.toArray(new String[0]);
		
		t= ir.getKeys2();
		t.addAll(er.getKeys2());
		String[] cells = t.toArray(new String[0]);
		
		ir.writeMM(Settings.S().getString(Settings.OUT_BASE)+".i", segids, cells);
		er.writeMM(Settings.S().getString(Settings.OUT_BASE)+".e", segids, cells);
	}
}
