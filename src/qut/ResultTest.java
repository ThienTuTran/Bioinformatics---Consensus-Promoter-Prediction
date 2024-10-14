package qut;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class ResultTest {
    private HashMap<String, String> defaultConsensus;

    @BeforeEach
    void declareDefaultConsensus() {
        defaultConsensus = new HashMap<>();
        defaultConsensus.put("all", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (5430 matches)");
        defaultConsensus.put("fixB", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (965 matches)");
        defaultConsensus.put("carA", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (1079 matches)");
        defaultConsensus.put("fixA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (896 matches)");
        defaultConsensus.put("caiF", " Consensus: -35: T T C A A A gap: 18.0 -10: T A T A A T  (11 matches)");
        defaultConsensus.put("caiD", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (550 matches)");
        defaultConsensus.put("yaaY", " Consensus: -35: T T G T C G gap: 18.0 -10: T A T A C T  (4 matches)");
        defaultConsensus.put("nhaA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (1879 matches)");
        defaultConsensus.put("folA", " Consensus: -35: T T G A C A gap: 17.5 -10: T A T A A T  (46 matches)");
    }
    @Test
    void ExplicitThreadingTest() throws InterruptedException {
        ExplicitThreading.main(null);
        for (Map.Entry<String, Sigma70Consensus> entry : ExplicitThreading.getConsensus().entrySet()) {
            assertEquals(defaultConsensus.get(entry.getKey()), entry.getValue().toString());
        }
    }
    @Test
    void ExecutorServiceTest() throws IOException, InterruptedException, ExecutionException {
        ExecutorService.main(null);
        for (Map.Entry<String, Sigma70Consensus> entry : ExecutorService.getConsensus().entrySet()) {
            assertEquals(defaultConsensus.get(entry.getKey()), entry.getValue().toString());
        }
    }
    @Test
    void ParallelStreamTest() throws IOException {
        ParallelStream.main(null);
        for (Map.Entry<String, Sigma70Consensus> entry : ParallelStream.getConsensus().entrySet()) {
            assertEquals(defaultConsensus.get(entry.getKey()), entry.getValue().toString());
        }
    }
}
