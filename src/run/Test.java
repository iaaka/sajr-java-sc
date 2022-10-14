package run;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import util.Settings;

public class Test {
	String version = "$LastChangedDate: 2014-09-08 12:04:56 +0400 (Mon, 08 Sep 2014) $\n$LastChangedRevision: 5128 $".replace("$", "");
	
	public static void main(String[] args) throws FileNotFoundException  {
//		Run.main(new String[]{"count_reads",
//				"-ann_in=/Users/pm19/tmp/sajr.ss/genes.sajr",
//				"-batch_in=/Users/pm19/tmp/sajr.ss/WS_D_SKNsp11418114.bam",
//				"-batch_out=/Users/pm19/tmp/sajr.ss/WS_D_SKNsp11418114"});
		//Run.main(new String[]{"annotate"});
		//Run.main(new String[]{"sajrcomp"});
		
		SamReader in = SamReaderFactory.makeDefault().open(SamInputResource.of(new BufferedInputStream(new FileInputStream("/Users/pm19/Downloads/sajr/1sorted.bam"),10000000)));
		System.out.println(!in.getFileHeader().getAttribute("SO").equals("coordinate"));
		System.out.println(in.iterator().next().getCharacterAttribute("as")==null);
		
	}
}
