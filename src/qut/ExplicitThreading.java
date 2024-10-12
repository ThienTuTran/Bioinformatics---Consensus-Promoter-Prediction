package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.locks.ReentrantLock;

class MyRunnable implements Runnable
{
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];
    private final String gbkFile;
    private static String referenceFile;
    static ReentrantLock lock = new ReentrantLock();

    public MyRunnable(String gbkFile, String referenceFile) {
        this.gbkFile = gbkFile;
        this.referenceFile = referenceFile;
    }

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
        return BioPatterns.getBestMatch(sigma70_pattern, upStreamRegion.toString());
    }

    private static GenbankRecord Parse(String file) throws IOException
    {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    public void run()
    {
        List<Gene> referenceGenes = null;
        GenbankRecord record = null;

        try {
            referenceGenes = ParseReferenceGenes(referenceFile);
            record = Parse(gbkFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(gbkFile);

        for (Gene referenceGene : referenceGenes) {
            System.out.println(referenceGene.name);
            for (Gene gene : record.genes) {
                if (Homologous(gene.sequence, referenceGene.sequence)) {
                    lock.lock();
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                    Match prediction = PredictPromoter(upStreamRegion);
                    if (prediction != null) {
                        consensus.get(referenceGene.name).addMatch(prediction);
                        consensus.get("all").addMatch(prediction);
                    }
                    lock.unlock();
                }
            }
        }
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }
}

class ExplicitThreading {
    private static final HashMap<String, Sigma70Consensus> consensus = new HashMap<>();
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
    public static void main(String[] args) throws InterruptedException {

        long startTime = System.currentTimeMillis();
        List<Thread> threads = new ArrayList<>();
        List<String> listGenBankFiles = ListGenbankFiles("./Ecoli");

        //create and start 4 threads for 4 files
        for (int i = 0; i < listGenBankFiles.size(); i++) {
            String gbkFile = listGenBankFiles.get(i);
            Thread thread = new Thread(new MyRunnable(gbkFile, "./referenceGenes.list"));
            thread.setName("Thread-" + (i + 1));
            threads.add(thread);
            thread.start();
        }

        //wait for all threads to complete
        for (Thread thread : threads)  thread.join();

        long timeLapsed = System.currentTimeMillis() - startTime;
        System.out.println("\nTime: " + timeLapsed/1000 + " s");
    }
}

