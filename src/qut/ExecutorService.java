package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class ExecutorService {
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

    public class RunnableTask implements Runnable {
        private final Gene referenceGene;
        private final Gene gene;
        private final GenbankRecord record;

        public RunnableTask(Gene referenceGene, Gene gene, GenbankRecord record) {
            this.referenceGene = referenceGene;
            this.gene = gene;
            this.record = record;
        }
        @Override
        public void run() {
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                Match prediction = PredictPromoter(upStreamRegion);
                if (prediction != null) {
                    consensus.get(referenceGene.name).addMatch(prediction);
                    consensus.get("all").addMatch(prediction);
                }
            }
        }
    }
    public void run(String referenceFile, String dir, int threadNum) throws FileNotFoundException, IOException, ExecutionException, InterruptedException {
        // Create a fixed thread pool based on the provided thread number
        java.util.concurrent.ExecutorService executorService = Executors.newFixedThreadPool(threadNum);
        System.out.println("Number of Threads: " + threadNum);

        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<String> genBankFiles = ListGenbankFiles(dir);

        // List to store future tasks
        List<Future<?>> futureTasks = new ArrayList<>();

        // Iterate over reference genes and files to submit tasks
        for (Gene referenceGene : referenceGenes) {
            for (String filename : genBankFiles) {
                GenbankRecord record = Parse(filename);

                // For each gene in the record, submit a task to the executor
                for (Gene gene : record.genes) {
                    Future<?> futureTask = executorService.submit(new RunnableTask(referenceGene, gene, record));
                    futureTasks.add(futureTask);
                }
            }
        }
        System.out.println("ExecutorService\n...Executing...\n");

        // shut down the executor service after submitting all tasks
        executorService.shutdown();

        // Wait for all tasks to complete
        for (Future<?> futureTask : futureTasks) {
            futureTask.get();  // Ensures that the task is completed
        }
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet()) {
            System.out.println(entry.getKey() + " " + entry.getValue());
        }
    }

    public static void main(String[] args) throws FileNotFoundException, IOException, ExecutionException, InterruptedException {
        long startTime = System.currentTimeMillis();
        new ExecutorService().run("./referenceGenes.list", "./Ecoli", 16);
        long timeLapsed = System.currentTimeMillis() - startTime;
        System.out.println("\nTime: " + timeLapsed/1000 + "s");
    }
}
