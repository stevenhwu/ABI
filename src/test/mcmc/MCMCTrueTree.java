package test.mcmc;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeLoggerWithTrueHaplotype;
import srp.haplotypes.HaplotypeModel;
import srp.likelihood.ShortReadLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.ext.TreeLikelihoodExt;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundLikelihoodParser;

public class MCMCTrueTree {

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
	public void testMCMCFixTree() throws Exception {
		String dataDir = "/home/sw167/workspace/ABI/data/";
//		String truePhylogenyFile = "H6_005_true_tree.trees";
//		String shortReadFile = "H6_srp.fasta";
//		
//		String dataDir = "/home/sw167/workspace/ABI/data/H7/";
		String truePhylogenyFile = "H7Srp.tree";
		String shortReadFile = "H7Srp.fasta";
		String trueHaplotypeFile = "H7Srp_fullHaplotype.fasta";
		
		String prefix = "FixTree_H7";
		String logTracerName = prefix+".log";
		String logTreeName = prefix+".trees";
		String logHaplotypeName = prefix+"_haplatype.hap";
		String operatorAnalysisFile = prefix+"_operatorAnalysisFile";
		
		int numberOfHaplotype = 7;
		int logInterval = 1000;
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);

		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
		HaplotypeModel haplotypeModel = new HaplotypeModel(alignmentMapping, numberOfHaplotype);

		// coalescent
		Parameter popSize = new Parameter.Default(
				ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);
		ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.DAYS);

		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), startingTree);
		coalescent.setId("coalescent");

		
		// clock model
		Parameter rateParameter = new Parameter.Default(StrictClockBranchRates.RATE, 1e-5, 0, 1);
		StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
		// sub model
		Parameter freqs = new Parameter.Default("frequency", haplotypeModel.getStateFrequencies());
		
		// treeLikelihood
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);
		TreeLikelihoodExt treeLikelihood = MCMCUtils.setupTreeLikelihood(kappa, freqs,
				haplotypeModel, treeModel, branchRateModel);

		// ShortReadLikelihood
		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(haplotypeModel);

		// CompoundLikelihood
		HashMap<String, Likelihood> compoundlikelihoods = MCMCUtils.setupCompoundLikelihood(
				popSize, kappa, coalescent, treeLikelihood, srpLikelihood);
		
		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
		Likelihood shortReadLikelihood = compoundlikelihoods.get(ShortReadLikelihood.SHORT_READ_LIKELIHOOD);
		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);
		
		// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		schedule.addOperators(MCMCUtils.defalutOperators(haplotypeModel, freqs, kappa, popSize));
//		schedule = defalutTreeOperators(schedule, treeModel);
		Parameter rootHeight = treeModel.getRootHeightParameter();
		
		
		double expectedInit = shortReadLikelihood.getLogLikelihood();
		assertEquals(expectedInit, srpLikelihood.getLogLikelihood(), 0);
		
		MCLogger[] loggers = new MCLogger[4];
//		loggers[0] = new MCLogger(formatter, logInterval, false);
		
		loggers[0] = new MCLogger(logTracerName, logInterval, false, 0);
		MCMCUtils.addToLogger(loggers[0], posterior, prior, likelihood, shortReadLikelihood,
				rootHeight, 
				//rateParameter,
				popSize, kappa, coalescent,
				freqs
				);

		loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
		MCMCUtils.addToLogger(loggers[1], posterior, prior, likelihood, shortReadLikelihood,
				popSize, kappa, coalescent, 
				rootHeight
				);
		
		TabDelimitedFormatter treeFormatter = new TabDelimitedFormatter(
				new PrintWriter(new FileOutputStream(new File(logTreeName))));

		loggers[2] = new TreeLogger(treeModel, branchRateModel, null, null,
				treeFormatter, logInterval, true, true, true, null, null);

		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		loggers[3] = new HaplotypeLoggerWithTrueHaplotype(haplotypeModel, trueAlignment, logHaplotypeName, logInterval*10);
		
		// MCMC
		
		MCMCOptions options = setMCMCOptions(logInterval);
		
		MCMC mcmc = new MCMC("mcmc1");
		mcmc.setShowOperatorAnalysis(true);
		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
		
		mcmc.init(options, posterior, schedule, loggers);
		mcmc.run();

		System.out.println(mcmc.getTimer().toString());
		// Tracer
		// List<Trace> traces = formatter.getTraces();
		// ArrayTraceList traceList = new ArrayTraceList("test", traces, 0);
		//
		//
		//
		// // Trace trace = traces.get(0);
		// for (Trace trace : traces) {
		// if (trace.getName().equals("ShortReadLikelihood")) {
		//
		// double startValue = (Double) trace.getValue(0);
		// double endValue = (Double) trace
		// .getValue(trace.getValuesSize() - 1);
		// assertEquals(expectedInit , startValue,0);
		// assertTrue(endValue > startValue);
		// // System.out.println(trace.getName());
		// // break;
		// }
		// }

		// for (int j = 0; j < trace.getValuesSize(); j++) {
		// System.out.print(trace.getValue(j) +"\t");
		// }
		// System.out.println();
		// System.out.println(Arrays.toString(trace.getRange()));
		// System.out.println(trace.getTraceT9ype());

	}

	private static MCMCOptions setMCMCOptions(int logInterval) {
		MCMCOptions options = new MCMCOptions(logInterval * 300,
				logInterval * 0, 1, MarkovChain.EVALUATION_TEST_THRESHOLD,
				false, logInterval * 5, 1.0);
//		options.setChainLength(logInterval * 300);;
//		options.setUseCoercion(true); // autoOptimize = true
//		options.setCoercionDelay(logInterval * 2);
//		options.setTemperature(1.0);
//		options.setFullEvaluationCount(logInterval*0);

		return options;
	}

}
