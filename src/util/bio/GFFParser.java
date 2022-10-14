package util.bio;

import java.io.BufferedReader;
import java.io.*;
import java.io.FileReader;
import java.io.IOException;
import java.util.zip.GZIPInputStream;

public class GFFParser{
	private BufferedReader in;
	private String line;
	
	/**
	 * parses GFF file, omits all lines started with '#' 
	 * @param fname
	 * @throws IOException
	 */
	public GFFParser(String fname) throws IOException {
		if(fname.endsWith(".gz")){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fname))),1000000);
		}else
			in = new BufferedReader(new FileReader(fname),1000000);
		readTillNextLine();
	}
	
	private void readTillNextLine() throws IOException {
		for(line = in.readLine();line !=null && line.startsWith("#");line = in.readLine()) {}
	}
	
	public void close() throws IOException {
		in.close();
	}

	/**
	 * 
	 * @return null if it's the end of the file
	 * @throws GFFException
	 */
	public GFFeature next() throws GFFException {
		if(line == null)
			return null;
		GFFeature t = new GFFeature(line);
		try {
			readTillNextLine();
		} catch (IOException e) {
			throw new GFFException(e.getMessage(),line);
		}
		return t;
	}
}
