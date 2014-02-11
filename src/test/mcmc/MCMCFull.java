package test.mcmc;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

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
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
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
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundLikelihoodParser;

public class MCMCFull {

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
	public void testMCMCFull() throws Exception {

		String dataDir = "/home/sw167/workspace/ABI/data/";

		String shortReadFile = "H7Srp.fasta";
		String trueHaplotypeFile = "H7Srp_fullHaplotype.fasta";
		
		String prefix = "FullTree_H7";
		String logTracerName = prefix+".log";
		String logTreeName = prefix+".trees";
		String logHaplotypeName = prefix+"_haplatype.hap";
		String operatorAnalysisFile = prefix+"_operatorAnalysisFile.txt";
		
		int numberOfHaplotype = 7;
		int logInterval = 1000;

		DataImporter dataImporter = new DataImporter(dataDir);

		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
		HaplotypeModel haplotypeModel = new HaplotypeModel(alignmentMapping, numberOfHaplotype);

		// coalescent
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);

//		 Random treeModel
		ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);
		ConstantPopulation constant = (ConstantPopulation) startingTree.getDemographicFunction();
		CoalescentSimulator simulator = new CoalescentSimulator();
		Tree tree = simulator.simulateTree(haplotypeModel, constant);
		TreeModel treeModel = new TreeModel(tree);// treeModel

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
		ArrayList<MCMCOperator> defalutOperatorsList = MCMCUtils.defalutOperators(haplotypeModel, freqs, kappa, popSize);
		schedule.addOperators(defalutOperatorsList);
		schedule.addOperators(MCMCUtils.defalutTreeOperators(treeModel));
		Parameter rootHeight = treeModel.getRootHeightParameter();
		
		int total = 0;
		for (int i = 0; i < schedule.getOperatorCount(); i++) {
			MCMCOperator operator = schedule.getOperator(i);
			total += operator.getWeight() ;
		}
		System.out.println("totalWeight: "+total);
		
		
		// MCLogger
		MCLogger[] loggers = new MCLogger[4];
		
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
		

	}

	private static MCMCOptions setMCMCOptions(int logInterval) {
//		MCMCOptions options = new MCMCOptions();
		MCMCOptions options = new MCMCOptions(logInterval * 500,
				logInterval * 0, 1, MarkovChain.EVALUATION_TEST_THRESHOLD,
				false, logInterval * 2, 1.0);

//		options.setChainLength(logInterval * 500);;
//		options.setUseCoercion(false); // autoOptimize = true
//		options.setCoercionDelay(logInterval * 2);
//		options.setTemperature(1.0);
//		options.setFullEvaluationCount(logInterval*0);

		return options;
	}

}
