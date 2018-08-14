package edu.unc.genomics.ngs;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

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
 * For each interval in Loci file (Bed), output overlapping from Input file (Wig) to a BedGraph that includes the interval name/id.
 * @author davidcking
 *
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
		log.debug("Initializing input file");
		int count = 0, skipped = 0;
		try (
               WigFileReader wig = WigFileReader.autodetect(inputFile);
               IntervalFileReader<? extends Interval> intervals = IntervalFileReader.autodetect(lociFile) // no ending semicolon in try-with-resources block
            )
        {
			log.debug("Iterating over all intervals and writing Wig for each");
            //TrackHeader header = TrackHeader.newWiggle();
			for (Interval interval : intervals) {
                System.out.printf("encountered interval: %s\n", interval);
                try {
                    Contig query = wig.query(interval);
                    /*header.setName(interval.getId());
                    try (WigFileWriter writer = new WigFileWriter(outputFile, header)) {
                      writer.write(query);
                    }*/
                    float[] values = query.getValues();
                    System.out.printf("interval %1$s has %2$d datapoints%n", interval, values.length);
                    final String chr = interval.getChr();
                    final String id = interval.getId();
                    for (int startPos = interval.low(); startPos < interval.high(); startPos++) {
                        int i = startPos - interval.low();
                        System.out.printf("%1$s\t%2$d\t%3$d\t%4$f\t%5$s\n", chr, startPos, startPos+1, values[i], id); 
                        
                    }
            
                } catch (WigFileException e) {
                    log.info("Skipping interval "+interval+" which has no data");
                    skipped++;
                }
            count++;
            }
		}
		
		log.info(count + " intervals processed");
		log.info(skipped + " interval skipped");
	}
	
	public static void main(String[] args) {
		new SplitWigIntervalsToBedGraphPlus().instanceMain(args);
	}

}
