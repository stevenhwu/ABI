package test.srp.haplotypes.operator;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.SingleBaseFrequencyOperator;
import srp.likelihood.ShortReadLikelihood;
import test.TestUtils;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.Trace;

public class SingleBaseFrequencyOperatorTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testGetOperatorName() {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, 3);

		SimpleMCMCOperator operator = new SingleBaseFrequencyOperator(haplotypeModel, null);
    	assertEquals(operator.getOperatorName(), "SingleBaseFrequencyOperator");
    	assertEquals(operator.getPerformanceSuggestion(), "");
    	
	}

	
	@Test
	public void testDoOperation() throws OperatorFailedException {
		String[] seqs = new String[]{
				"GGGGGGGGGGGGG.....",
				".....CCCCCCCCCCCCCCCCCCCCCCCCCCC....",
				"..........GGGGGGGGGGGGGGCGCGGGGGGGGG",
//				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		
		String[] haps = new String[]{
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA"
//				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
				"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
				};
		
		HaplotypeModel haplotypeModel = AlignmentUtils.createHaplotypeModel(seqs, haps);
		Parameter freqs = new Parameter.Default("frequency", new double[] {0.0,0.25,0.75,0.0});
		SimpleMCMCOperator operator = new SingleBaseFrequencyOperator(haplotypeModel, freqs);
    	
		int[] count = new int['T'+1]; 
    	for (int i = 0; i < 100; i++) {
    		operator.doOperation();
    		
    		int newChar = haplotypeModel.getSwapInfo().getSwapInfoSWAPBASE()[2];
    		count[newChar]++;
    		
			String newHap = haplotypeModel.getHaplotypeString(0);
			assertTrue( newHap.contains("C") || newHap.contains("G"));
			haplotypeModel.reject();
			newHap = haplotypeModel.getHaplotypeString(0);
			assertEquals(haps[0], newHap);

		}
    	TestUtils.assertExpectationRange(count['C'], 25, 10);
    	TestUtils.assertExpectationRange(count['G'], 75, 10);
    	
    	
    	freqs.setParameterValue(1, 0);
    	freqs.setParameterValue(3, 0.25);
		count = new int['T'+1]; 
		int ite = 1000;
    	for (int i = 0; i < ite; i++) {
    		operator.doOperation();
    		
    		int newChar = haplotypeModel.getSwapInfo().getSwapInfoSWAPBASE()[2];
    		count[newChar]++;
    		
			String newHap = haplotypeModel.getHaplotypeString(0);
			assertTrue( newHap.contains("G") || newHap.contains("T"));
			haplotypeModel.reject();
			newHap = haplotypeModel.getHaplotypeString(0);
			assertEquals(haps[0], newHap);

		}
    	TestUtils.assertExpectationRange(count['G'], 0.75*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['T'], 0.25*ite, 0.05*ite);


	}
	@Test
	public void testDoOperationMCMC() {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		String[] haps = new String[]{
				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		Parameter freqs = new Parameter.Default("frequency", haplotypeModel.getStateFrequencies());
		
		SimpleMCMCOperator operator = new SingleBaseFrequencyOperator(haplotypeModel, freqs);
		MCMCOperator[] operators = new MCMCOperator[]{operator};
    	// Operators
    	OperatorSchedule schedule = new SimpleOperatorSchedule();
    	schedule.addOperator(operator);
    	
    	//CompoundLikelihood
        ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(haplotypeModel);
		double initLikelihood = srpLikelihood.getLogLikelihood();
    	
		ArrayLogFormatter formatter = new ArrayLogFormatter(false);
    	MCLogger[] loggers = new MCLogger[1];
    	loggers[0] = new MCLogger(formatter, 1, false);
    	loggers[0].add(srpLikelihood );

//    	// MCMC
//    	MCMCOptions options = new MCMCOptions();
//    	options.setChainLength(100);
////    	mcmc.setShowOperatorAnalysis(true);
    	
    	MCMC mcmc = new MCMC("mcmc1");
		mcmc.init(100, srpLikelihood, operators , loggers);
    	mcmc.setShowOperatorAnalysis(false);
    	mcmc.run();
    	
        for (int i = 0; i < schedule.getOperatorCount(); i++) {
        	MCMCOperator op = schedule.getOperator(i);
        	double acceptanceProb = MCMCOperator.Utils.getAcceptanceProbability(op);
        	assertTrue("~0 AcceptanceProb:"+acceptanceProb, acceptanceProb > 0.01);
        	assertTrue("~1 AcceptanceProb:"+acceptanceProb, acceptanceProb < 0.99);
        }
//        System.out.println(MCMCOperator.Utils.getAcceptanceProbability(operator));
        
        
    	List<Trace> traces = formatter.getTraces();
		for (Trace<?> trace : traces) {
			if (trace.getName().equals("ShortReadLikelihood")) {
				double startValue = (Double) trace.getValue(0);
				double endValue = (Double) trace.getValue(trace.getValuesSize() - 1);
				assertEquals(initLikelihood , startValue, 0);
				assertTrue(endValue > startValue);
			}
		}

	}



}