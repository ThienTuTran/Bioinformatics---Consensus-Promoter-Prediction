package qut;

public class ThreadTaskForPS {
    private final Gene referenceGene;
    private final Gene gene;
    private final GenbankRecord record;
    public ThreadTaskForPS(Gene referenceGene, Gene gene, GenbankRecord record) {
        this.referenceGene = referenceGene;
        this.gene = gene;
        this.record = record;
    }
    public Gene getReferenceGene() {
        return referenceGene;
    }
    public Gene getGene() {
        return gene;
    }
    public GenbankRecord getRecord() {
        return record;
    }
}

