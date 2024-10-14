package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;

public class ParallelStream {
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static final ThreadLocal<Series> sigma70_pattern =
            ThreadLocal.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];

    static
    {
        complement['C'] = 'G'; complement['c'] = 'g';
        complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a';
        complement['A'] = 'T'; complement['a'] = 't';
    }
    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true)
        {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    private static boolean Homologous(PeptideSequence A, PeptideSequence B)
    {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene)
    {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location-1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
        else
        {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i=0; i<upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart-i]];
            return new NucleotideSequence(result);
        }
    }

    private static Match PredictPromoter(NucleotideSequence upStreamRegion)
    {
        return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir)
    {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    private static List<String> ListGenbankFiles(String dir)
    {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    private static GenbankRecord Parse(String file) throws IOException
    {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    // added
    public synchronized void addConsensus(String name, Match prediction) {
        consensus.get(name).addMatch(prediction);
        consensus.get("all").addMatch(prediction);
    }

    public void run3rd(String referenceFile, String dir, int threadNum) throws IOException
    {
        System.out.println("Number of Threads: " + threadNum);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",Integer.toString(threadNum));
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);

        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            for (Gene referenceGene :referenceGenes) {
                System.out.println(referenceGene.name);
                // Parse the GenBank file
                try {
                    GenbankRecord record = Parse(filename);

                    // Parallel stream to process each gene in the record
                    record.genes.parallelStream().forEach(gene -> {
                        if (Homologous(gene.sequence, referenceGene.sequence)) {
                            NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                            Match prediction = PredictPromoter(upStreamRegion);
                            if (prediction != null) {
                                addConsensus(referenceGene.name, prediction);
                            }
                        }
                    });
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }

    public void runWithPrep(String referenceFile, String dir, int threadNum) throws IOException {
        // Preparing Data and store in List<TaskHandler>
        List<TaskHandler> taskHandlers = new ArrayList<>();
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    taskHandlers.add(new TaskHandler(referenceGene, gene, record));
                }
            }
        }
        System.out.println("parallelStream - Preparing data finished!\nComputing...\n");

        // Initialize Threads with ParallelStream() and compute
        // System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(Runtime.getRuntime().availableProcessors()));
        System.out.println("Now run on " + threadNum + " threads.");
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(threadNum));
        taskHandlers.parallelStream()
                .filter(task -> Homologous(task.getGene().sequence, task.getReferenceGene().sequence))
                .forEach(task -> {
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(task.getRecord().nucleotides, task.getGene());
                    Match prediction = PredictPromoter(upStreamRegion);
                    if (prediction != null) {
                        addConsensus(task.getReferenceGene().name, prediction);
                    }
                });

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }

    public static void main(String[] args) throws IOException{
        long startTime = System.currentTimeMillis();
        //new ParallelStream().run3rd("./referenceGenes.list", "./Ecoli", 16);
        new ParallelStream().runWithPrep("./referenceGenes.list", "./Ecoli", 16);
        long timeLapsed = System.currentTimeMillis() - startTime;
        System.out.println("\nTime: " + timeLapsed/1000.0 + "s");
    }
}