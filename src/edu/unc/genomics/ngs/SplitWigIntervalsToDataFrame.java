package edu.unc.genomics.ngs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.*; // for Path, Paths, newBufferedWriter
//import java.nio.file.Path;
//import java.nio.file.Paths;

import org.apache.log4j.Logger;

import com.beust.jcommander.Parameter;

import edu.unc.genomics.CommandLineTool;
import edu.unc.genomics.CommandLineToolException;
import edu.unc.genomics.Interval;
import edu.unc.genomics.Contig;
import edu.unc.genomics.ReadablePathValidator;
import edu.unc.genomics.io.IntervalFileReader;
import edu.unc.genomics.io.WigFileReader;
import edu.unc.genomics.io.WigFileWriter;
import edu.unc.genomics.io.WigFileException;
import edu.ucsc.genome.TrackHeader;

/**
 * For each interval in Loci file (Bed), output min, mean, max, N of data overlapping from Input file (Wig) to a BedGraph that includes the interval name/id.
 * @author davidcking
 * modified from timothypalpant
 */
public class SplitWigIntervalsToBedGraphPlus extends CommandLineTool {

	private static final Logger log = Logger.getLogger(SplitWigIntervalsToBedGraphPlus.class);

	@Parameter(names = {"-i", "--input"}, description = "Input file (Wig)", 
             required = true, validateWith = ReadablePathValidator.class)
	public Path inputFile;
	@Parameter(names = {"-l", "--loci"}, description = "Loci file (Bed)", 
             required = true, validateWith = ReadablePathValidator.class)
	public Path lociFile;
	@Parameter(names = {"-o", "--output"}, description = "Output bedgraph-plus file name")
	public Path outputFile;
	
	@Override
	public void run() throws IOException {
        log.debug("Initializing output file");
        Charset charset = Charset.forName("US-ASCII");
        BufferedWriter writer = Files.newBufferedWriter(outputFile, charset);
		log.debug("Initializing input file");
		int count = 0, skipped = 0;
		try (
               WigFileReader wig = WigFileReader.autodetect(inputFile);
               IntervalFileReader<? extends Interval> intervals = IntervalFileReader.autodetect(lociFile) // no ending semicolon in try-with-resources block
            )
        {
			log.debug("Iterating over all intervals and writing Wig min,mean,max,N for each");
            int interval_i = 0;
			for (Interval interval : intervals) {
                log.info("encountered interval: " + interval);
                try {
                    Contig query = wig.query(interval);
                    float[] values = query.getValues();
                    final String chr = interval.getChr();
                    String id = interval.getId();
                    if (id == null) {
                        id = "interval_" + interval_i;
                    }
                    float i_max = Float.MIN_VALUE;
                    float i_min = Float.MAX_VALUE;
                    float i_sum = 0;
                    int i_N = 0;
                    for (int startPos = interval.low(); startPos < interval.high(); startPos++) {
                        int i = startPos - interval.low();
                        i_max = Math.max(i_max, values[i]);
                        i_min = Math.min(i_min, values[i]);
                        i_sum += values[i];
                        i_N++;
                    }
                    float i_mean = i_sum / i_N;
                    int startPos = interval.low();
                    String line = String.format("%1$s\t%2$d\t%3$d\t%4$f\t%5$f\t%6$f\t%7$d\t%8$s\n", chr, startPos, startPos+1, i_min, i_mean, i_max, i_N, id);
                                                                                                  /* 1      2          3        4        5      6     7    8 */
                    writer.write(line);
                    
                    interval_i++; 
                } catch (WigFileException e) {
                    log.info("Skipping interval "+interval+" which has no data");
                    skipped++;
                }
            count++;
            }
		}
        writer.flush();
        writer.close();		
		log.info(count + " intervals processed");
		log.info(skipped + " interval skipped");
	}
	
	public static void main(String[] args) {
		new SplitWigIntervalsToBedGraphPlus().instanceMain(args);
	}

}
