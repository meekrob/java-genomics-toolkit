package edu.unc.genomics.wigmath;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;
import org.broad.igv.bbfile.WigItem;

import com.beust.jcommander.Parameter;

import edu.unc.genomics.CommandLineTool;
import edu.unc.genomics.CommandLineToolException;
import edu.unc.genomics.io.WigFile;
import edu.unc.genomics.io.WigFileException;
import edu.unc.utils.FloatCorrelation;

public class WigCorrelate extends CommandLineTool {

	private static final Logger log = Logger.getLogger(WigCorrelate.class);

	@Parameter(description = "Input files", required = true)
	public List<String> inputFiles = new ArrayList<String>();
	@Parameter(names = {"-w", "--window"}, description = "Window size (bp)")
	public Integer windowSize = 100;
	@Parameter(names = {"-t", "--type"}, description = "Correlation metric to use (pearson/spearman)")
	public String type = "pearson";
	@Parameter(names = {"-o", "--output"}, description = "Output file")
	public Path outputFile;
	
	private Correlation corr;
	private List<WigFile> wigs = new ArrayList<>();
	private List<String> chromosomes;
	int[] chrStarts, chrStops, chrLengths, nBins;
	int totalNumBins = 0;
	private float[][] correlationMatrix;
	
	@Override
	public void run() throws IOException {
		if (inputFiles.size() < 2) {
			throw new CommandLineToolException("Cannot correlate < 2 input files.");
		}
		
		corr = Correlation.fromName(type);
		if (corr == null) {
			log.error("Unknown correlation metric: "+type);
			throw new CommandLineToolException("Unknown correlation metric: "+type+". Options are pearson, spearman");
		} else {
			log.debug("Computing "+type+" correlation");
		}
		
		log.debug("Initializing input files");
		for (String inputFile : inputFiles) {
			try {
				wigs.add(WigFile.autodetect(Paths.get(inputFile)));
			} catch (IOException | WigFileException e) {
				log.error("Error initializing input Wig file: " + inputFile);
				e.printStackTrace();
				throw new CommandLineToolException(e.getMessage());
			}
		}
		log.debug("Initialized " + wigs.size() + " input files");
		correlationMatrix = new float[wigs.size()][wigs.size()];
		// Put ones down the diagonal
		for (int i = 0; i < correlationMatrix.length; i++) {
			correlationMatrix[i][i] = 1;
		}
		
		// Get the maximum extent for each chromosome
		chromosomes = new ArrayList<>(WigMathTool.getCommonChromosomes(wigs));
		chrStarts = new int[chromosomes.size()];
		Arrays.fill(chrStarts, Integer.MAX_VALUE);
		chrStops = new int[chromosomes.size()];
		Arrays.fill(chrStops, Integer.MIN_VALUE);
		for (int i = 0; i < chromosomes.size(); i++) {
			String chr = chromosomes.get(i);
			for (WigFile w : wigs) {
				if (w.getChrStart(chr) < chrStarts[i]) {
					chrStarts[i] = w.getChrStart(chr);
				}
				if (w.getChrStop(chr) > chrStops[i]) {
					chrStops[i] = w.getChrStop(chr);
				}
			}
		}
		// Calculate the number of bins for each chromosome
		chrLengths = new int[chromosomes.size()];
		nBins = new int[chromosomes.size()];
		for (int i = 0; i < chromosomes.size(); i++) {
			chrLengths[i] = chrStops[i] - chrStarts[i] + 1;
			nBins[i] = (int) Math.ceil(((double)chrLengths[i])/windowSize);
			totalNumBins += nBins[i];
		}
		log.debug("Total number of bins for all chromosomes = "+totalNumBins);
		
		// Compute the pairwise correlations between all files
		// only keeping two files in memory at any one time
		for (int i = 0; i < wigs.size()-1; i++) {
			log.debug("Loading data from file "+i);
			float[] binsI = getDataVector(wigs.get(i));
			
			for (int j = i+1; j < wigs.size(); j++) {
				log.debug("Loading data from file "+j);
				float[] binsJ = getDataVector(wigs.get(j));
				
				log.debug("Correlating ("+i+","+j+")");
				switch (corr) {
				case PEARSON:
					correlationMatrix[i][j] = FloatCorrelation.pearson(binsI, binsJ);
					break;
				case SPEARMAN:
					correlationMatrix[i][j] = FloatCorrelation.spearman(binsI, binsJ);
					break;
				}
				
				// Copy the correlation to the symmetric point in the matrix
				correlationMatrix[j][i] = correlationMatrix[i][j];
			}
		}
		
		// Write the correlation matrix to output
		StringBuilder output = new StringBuilder(type);
		// Header row
		for (int i = 0; i < wigs.size(); i++) {
			output.append("\t"+wigs.get(i).getPath().getFileName());
		}
		// Data rows
		for (int i = 0; i < wigs.size(); i++) {
			output.append("\n"+wigs.get(i).getPath().getFileName());
			for (int j = 0; j < wigs.size(); j++) {
				output.append("\t"+correlationMatrix[i][j]);
			}
		}
		
		if (outputFile != null) {
			log.debug("Writing to output file");
			try (BufferedWriter writer = Files.newBufferedWriter(outputFile, Charset.defaultCharset())) {
				writer.write(output.toString());
			}
		} else {
			System.out.println(output.toString());
		}
	}
	
	/**
	 * Get a collapsed data vector from w
	 * This is the "binning" part of ACT
	 * @param w a WigFile
	 * @return a vector of binned values from w in a consistent order
	 * @throws IOException 
	 * @throws WigFileException 
	 */
	private float[] getDataVector(WigFile w) {
		float[] values = new float[totalNumBins];
		int[] counts = new int[totalNumBins];
		
		int binOffset = 0;
		for (int i = 0; i < chromosomes.size(); i++) {
			String chr = chromosomes.get(i);
			try {
				Iterator<WigItem> result = w.query(chr, w.getChrStart(chr), w.getChrStop(chr));
				while (result.hasNext()) {
					WigItem item = result.next();
					for (int bp = item.getStartBase(); bp <= item.getEndBase(); bp++) {
						int bin = (bp-chrStarts[i]) / windowSize;
						values[bin+binOffset] += item.getWigValue();
						counts[bin+binOffset]++;
					}
				}
			} catch (WigFileException | IOException e) {
				log.error("Error getting data from wig file: "+w.getPath());
				throw new CommandLineToolException(e.getMessage());
			}
			
			binOffset += nBins[i];
		}
		
		// Compute the average
		for (int i = 0; i < totalNumBins; i++) {
			if (counts[i] > 0) {
				values[i] /= counts[i];
			} else {
				values[i] = Float.NaN;
			}
		}
		
		return values;
	}
	
	/**
	 * @param args
	 * @throws WigFileException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, WigFileException {
		new WigCorrelate().instanceMain(args);
	}
	
	public enum Correlation {
		PEARSON("pearson"),
		SPEARMAN("spearman");
		
		private String name;
		
		Correlation(final String name) {
			this.name = name;
		}
		
		public static Correlation fromName(final String name) {
			for (Correlation c : Correlation.values()) {
				if (c.getName().equalsIgnoreCase(name)) {
					return c;
				}
			}
			
			return null;
		}

		/**
		 * @return the name
		 */
		public String getName() {
			return name;
		}
	}

}
