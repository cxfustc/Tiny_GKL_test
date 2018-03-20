import java.io.*;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeBinding;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import com.intel.gkl.pairhmm.IntelPairHmm;
import com.intel.gkl.pairhmm.IntelPairHmmOMP;

public class PairHMMTest {
    public static void main (String[] args) {
        PairHMMRun pairHMMRun = new PairHMMRun();
        pairHMMRun.run_pairHMM();

        NewThread thread1 = new NewThread();
        NewThread thread2 = new NewThread();
        NewThread thread3 = new NewThread();
        NewThread thread4 = new NewThread();
        thread1.start();
        thread2.start();
        thread3.start();
        thread4.start();

        try {
            thread1.join();
        } catch (Exception e) {

        }

        try {
            thread2.join();
        } catch (Exception e) {

        }

        try {
            thread3.join();
        } catch (Exception e) {

        }

        try {
            thread4.join();
        } catch (Exception e) {

        }
    }
}

class NewThread extends Thread {
    public void run() {
        PairHMMRun pairHMMRun = new PairHMMRun();
        pairHMMRun.run_pairHMM();
    }
}

class PairHMMRun {
    private final int numReadRepeat = 1000;
    private final int numHapRepeat = 100;

    private ReadDataHolder[] readDataArray;
    private HaplotypeDataHolder[] mHaplotypeDataArray;
    private double[] mLogLikelihoodArray;

    private PairHMMNativeBinding piarhmm;

    private int FAIL_TO_START_PAIRHMM_ENGINE = 1;
    private int NUM_THRED_FOR_PAIRHMM = 2;

    private byte[] read = {
            'A', 'G', 'C', 'T', 'C', 'G', 'A', 'G', 'G', 'G',
            'A', 'G', 'A', 'G', 'A', 'G', 'T', 'T', 'T', 'T',
            'A', 'T', 'T', 'C', 'A', 'G', 'A', 'A', 'T', 'G',
            'G', 'G', 'A', 'A', 'A', 'C', 'T', 'T', 'A', 'C',
            'C', 'A', 'A', 'G', 'G', 'A', 'G', 'A', 'C', 'T',
            'A', 'C', 'T', 'A', 'T', 'G', 'A', 'A', 'T', 'C',
            'A', 'C', 'G', 'A', 'T', 'A', 'C', 'T', 'A', 'C',
            'G', 'A', 'T', 'G', 'A', 'T', 'C', 'C', 'T', 'C',
            'G', 'G', 'G', 'A', 'A', 'T', 'A', 'C', 'A', 'G',
            'G', 'G', 'A', 'T', 'T', 'A', 'C'
    };

    private byte[] q = {
            28, 30, 31, 31, 30, 23, 30, 32, 31, 31,
            30, 32, 29, 32, 30, 32, 28, 30, 30, 30,
            29, 29, 30, 30, 30, 32, 30, 30, 29, 31,
            30, 30, 30, 31, 31, 30, 31, 30, 29, 31,
            29, 30, 31, 32, 31, 30, 33, 30, 31, 31,
            29, 31, 31, 29, 30, 32, 30, 31, 30, 31,
            31, 31, 24, 31, 30, 30, 32, 32, 29, 33,
            25, 29, 29, 33, 29, 31, 31, 31, 32, 31,
            26, 32, 32, 30, 31, 29, 28, 33, 30, 31,
            30, 29, 25, 24, 25, 24, 29
    };

    private byte[] i = {
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 39, 39, 39, 37, 37,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 38, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 39, 39,
            38, 39, 40, 40, 40, 40, 42
    };

    private byte[] d = {
            40, 40, 40, 40, 40, 40, 40, 40, 40, 39,
            40, 40, 40, 40, 40, 39, 39, 39, 33, 33,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 34, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 41
    };

    private byte[] c = {
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10
    };

    private byte[] hap = {
            'A', 'G', 'C', 'T', 'C', 'G', 'A', 'G', 'G', 'G',
            'A', 'G', 'A', 'G', 'A', 'G', 'T', 'T', 'T', 'T',
            'A', 'T', 'T', 'C', 'A', 'G', 'A', 'A', 'T', 'G',
            'G', 'G', 'A', 'A', 'A', 'C', 'T', 'T', 'A', 'C',
            'C', 'A', 'A', 'G', 'G', 'A', 'G', 'A', 'C', 'T',
            'A', 'C', 'T', 'A', 'T', 'G', 'A', 'A', 'T', 'C',
            'A', 'C', 'G', 'A', 'T', 'A', 'C', 'T', 'A', 'C',
            'G', 'A', 'T', 'G', 'A', 'T', 'C', 'C', 'T', 'C',
            'G', 'G', 'G', 'A', 'A', 'T', 'A', 'C', 'A', 'G',
            'G', 'G', 'A', 'T', 'T', 'A', 'C', 'A', 'G', 'G',
            'A', 'A', 'T', 'G', 'A', 'T', 'C', 'C', 'T', 'T',
            'A', 'T', 'G', 'A', 'A', 'C', 'A', 'A', 'G', 'A',
            'T', 'A', 'T', 'T', 'A', 'G', 'G', 'G', 'A', 'A',
            'T', 'A', 'T', 'A', 'G', 'T', 'T', 'A', 'C', 'A',
            'G', 'G', 'C', 'A', 'A', 'A', 'A', 'A', 'A', 'T',
            'A', 'A', 'C', 'A', 'A', 'A', 'G', 'A', 'G', 'A',
            'A', 'C', 'G', 'T', 'G', 'A', 'A', 'A', 'G', 'A',
            'T', 'T', 'T', 'G', 'A', 'G', 'T', 'C', 'T', 'G',
            'A', 'C', 'C', 'G', 'G', 'G', 'A', 'C', 'A', 'G',
            'A', 'G', 'A', 'C', 'C', 'A', 'T', 'G', 'A', 'G',
            'A', 'G', 'G', 'A', 'G', 'G', 'C', 'C', 'G', 'A',
            'T', 'T', 'G', 'A', 'A', 'C', 'G', 'A', 'A', 'G',
            'T', 'C', 'A', 'A', 'A', 'G', 'T'
    };

    public void run_pairHMM() {
        readDataArray = new ReadDataHolder[numReadRepeat];
        mHaplotypeDataArray = new HaplotypeDataHolder[numHapRepeat];
        mLogLikelihoodArray = new double[numReadRepeat * numHapRepeat];

        piarhmm = new IntelPairHmmOMP();
        if (!piarhmm.load(null)) {
            System.exit(FAIL_TO_START_PAIRHMM_ENGINE);
        }

        PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.maxNumberOfThreads = NUM_THRED_FOR_PAIRHMM;
        args.useDoublePrecision = false;
        piarhmm.initialize(args);

        for (int j = 0; j < numReadRepeat; ++j) {
            readDataArray[j] = new ReadDataHolder();
            readDataArray[j].readBases = read;
            readDataArray[j].readQuals = q;
            readDataArray[j].insertionGOP = i;
            readDataArray[j].deletionGOP = d;
            readDataArray[j].overallGCP = c;
        }

        for (int j = 0; j < numHapRepeat; ++j) {
            mHaplotypeDataArray[j] = new HaplotypeDataHolder();
            mHaplotypeDataArray[j].haplotypeBases = hap;
        }

        long beg = System.nanoTime();
        piarhmm.computeLikelihoods(readDataArray, mHaplotypeDataArray, mLogLikelihoodArray);
        System.out.printf("Cost %fs\n", (double) (System.nanoTime() - beg) / 1000000000.0);
    }
}
