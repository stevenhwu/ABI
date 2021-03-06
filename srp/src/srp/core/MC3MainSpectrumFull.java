package srp.core;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.SpectraParameter.SpectraType;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.likelihood.spectrum.AbstractShortReadsSpectrumLikelihood.DistType;
import srp.likelihood.spectrum.ShortReadsSpectrumLikelihood;
import srp.likelihood.spectrum.treelikelihood.SpectrumTreeLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.mcmcmc.MCMCMC;
import dr.inference.mcmcmc.MCMCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundLikelihoodParser;


public class MC3MainSpectrumFull {

	public static void main(String[] args) throws Exception {
//
//0.30933196101825705	0.30249512435225206	0.2825730502398382	0.3122140438344868	0.31071974332425467	0.30928379887532603	0.31231278315168154	
//0.4313886921245247	0.42251157301604264	0.36601314864409784	0.4287955518793434	0.43182972173904316	0.431444897184788	0.43294820478526	
//0.31212654047288 17	0.30638171409639386	0.2880006593359247	0.3099814108119972	0.3123797592007158	0.3114106319137973	0.31290882634473327	
//0.37250962641044366	0.36800711790778085	0.3565109918690321	0.37268010551554437	0.37269166754477395	0.3730513482348595	0.370813206645306	
//0.540636490287451	0.5263861349379844	0.4270550208705616	0.538498800782783	0.5416501177059991	0.5405796169959256	0.5397434743780452	
//0.48127893944007205	0.4702431976343528	0.38967231908010985	0.481094963242507	0.4825315431397799	0.4821172766781288	0.48092511228341606	
//0.4281713946916407	0.4198173796641415	0.3524266561999545	0.42895690329065106	0.42928715725961364	0.42590787930341933	0.427904513786816	
//
//

		String dataDir;
		int runIndex;
		int totalSamples;
		int logInterval;
		int noOfTrueHaplotype;
		int noOfRecoveredHaplotype;
		boolean randomTree = true;
		boolean randomSpectrum = true;
		SpectraType randomSpectrumType = SpectraType.RANDOM;
		DistType distTypeCode = DistType.flat;
		
		boolean commandLine = true;
		commandLine = false;
		if(commandLine){
			dataDir = args[0];
			runIndex = Integer.parseInt(args[1]);
			totalSamples = Integer.parseInt(args[2]);
			logInterval = Integer.parseInt(args[3]);
			noOfTrueHaplotype = Integer.parseInt(args[4]);
			noOfRecoveredHaplotype= Integer.parseInt(args[5]);
		}
		else{	
			dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/testData/";
			runIndex = 54;
			dataDir += "H7_"+runIndex+"/";
			
			totalSamples = 10	 ;
			logInterval = 100;
			
			randomTree = true;
//			randomTree = false;
			
			randomSpectrum = true;
			randomSpectrumType = SpectraType.DOMINANT;
			randomSpectrumType = SpectraType.RANDOM ;
			
//			randomSpectrumType = SpectrumAlignmentModel.SpectrumType.EQUAL;
//			randomSpectrum = false;
			
			distTypeCode = DistType.betaMean;
			distTypeCode = DistType.flat;
			
			noOfTrueHaplotype = 7;
			noOfRecoveredHaplotype=7;
		}
		
		String hapRunIndex = "H"+noOfTrueHaplotype+"_"+runIndex;
		String prefix = dataDir+"FullTree_"+hapRunIndex;
		
		String shortReadFile = hapRunIndex +"_Srp.fasta";
		String trueHaplotypeFile = hapRunIndex +"_Srp_fullHaplotype.fasta";
//		String trueHaplotypeFile = "FullTree_H7_"+runIndex+".haplatypepartial";
//		String trueHaplotypeFile = prefix+".haplatypepartial";
		
		String logTracerName = prefix+".log";
		String logTreeName = prefix+".trees";
		String logHaplotypeName = prefix+".haplatype";
		String operatorAnalysisFile = prefix+"_operatorAnalysisFile.txt";
		
		String partialSpectrumName = hapRunIndex+".haplatypepartial";
		String partialTreeName = "FullTree_"+hapRunIndex+".treespartial";
//		String partialTreeName = hapRunIndex + "_Srp.tree";
		
		DataImporter dataImporter = new DataImporter(dataDir);

		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		ShortReadMapping srpMap = new ShortReadMapping(shortReads);
		int spectrumLength = srpMap.getLength();
//		 SpectrumModel and ShortReadLikelihood;
		
//		SpectrumAlignmentModel spectrumModel ;
//		ShortReadsSpectrumLikelihood srpLikelihood ;
//		
////		boolean redo = true;
////		int c = 0;
//		if (randomSpectrum) {
////			do {
//				spectrumModel = new SpectrumAlignmentModel(spectrumLength, noOfRecoveredHaplotype, randomSpectrumType);
//				srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap, distTypeCode);
////				redo = (srpLikelihood.getLogLikelihood() == Double.NEGATIVE_INFINITY);
////				c++;
////				if(c==100){
////					System.err.println("After 100 try, ShortreadLikelihood= "+ srpLikelihood.getLogLikelihood());
////					System.exit(-1);
////				}
////			} while (redo);
//		}		
//		else{
//			spectrumModel = dataImporter.importPartialSpectrumFile(partialSpectrumName );
//			srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap, distTypeCode);
//		}
////		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, noOfRecoveredHaplotype, randomSpectrumType);
////		srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap, distTypeCode);
//		// coalescent
//		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);
//		
//		// Random treeModel
//		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
//		TreeModel treeModel;
//		if(randomTree){
//			treeModel = MCMCSetupHelper.setupRandomTreeModel(popModel, spectrumModel, Units.Type.YEARS);
//		}
//		else{
//			Tree partialPhylogeny = dataImporter.importTree(partialTreeName);
//			treeModel = new TreeModel(TreeModel.TREE_MODEL, partialPhylogeny, false);
//		}
//
//		// Coalescent likelihood
//		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), popModel);
//		coalescent.setId("coalescent");
//
//		
//		// Simulate  treeLikelihood
//		HashMap<String, Object> parameterList = MCMCSetupHelperSpectrum.setupSpectrumTreeLikelihoodSpectrumModel(treeModel, spectrumModel);
//		Parameter kappa = (Parameter) parameterList.get("kappa");
//		Parameter freqs = (Parameter) parameterList.get("freqs");
//		Parameter rate = (Parameter) parameterList.get("rate");
//		StrictClockBranchRates branchRateModel = (StrictClockBranchRates) parameterList.get("branchRateModel");
//		SpectrumTreeLikelihood treeLikelihood = (SpectrumTreeLikelihood) parameterList.get("treeLikelihood");
//		
//		
//		
//		// CompoundLikelihood
//		HashMap<String, Likelihood> compoundlikelihoods = MCMCSetupHelper.setupCompoundLikelihood(
//				popSize, kappa, coalescent, treeLikelihood, srpLikelihood);
//		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
//		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
//		Likelihood shortReadLikelihood = compoundlikelihoods.get(ShortReadsSpectrumLikelihood.SHORT_READ_LIKELIHOOD);
//		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);
//
//		// Operators
//		OperatorSchedule schedule = new SimpleOperatorSchedule();
////		ArrayList<MCMCOperator> defalutOperatorsList = 
//		MCMCSetupHelperSpectrum.defalutSpectrumOperators(schedule, spectrumModel, freqs, popSize, kappa);
//		MCMCSetupHelperSpectrum.defalutTreeOperators(schedule, treeModel);
//		
//		
////		MCMCOperator operator;
////		operator = new RJTreeOperator(spectrumModel, treeModel);
////		operator.setWeight(100);
////		schedule.addOperator(operator);
//		
//		Parameter rootHeight = treeModel.getRootHeightParameter();
//		rootHeight.setId("rootHeight");
//		
////		double total = 0;
////		for (int i = 0; i < schedule.getOperatorCount(); i++) {
////			MCMCOperator op= schedule.getOperator(i);
//////			System.out.println(op.getOperatorName());
////			total += op.getWeight() ;
////		}
////		System.out.println("totalWeight: "+total);
//		
//
//		// MCLogger
//		MCLogger[] loggers = new MCLogger[4];
//		// log tracer
//		loggers[0] = new MCLogger(logTracerName, logInterval, false, 0);
//		MCMCSetupHelper.addToLogger(loggers[0], posterior, prior, likelihood, shortReadLikelihood,
//				rootHeight, 
//				rate,
//				popSize, kappa, coalescent,
//				freqs
//				);
//		// System.out
//		loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
//		MCMCSetupHelper.addToLogger(loggers[1],
////				freqs
//				posterior, prior, likelihood, shortReadLikelihood,
//				rate,
//				popSize, kappa, coalescent, rootHeight
//				);
//		// log Tree
//		TabDelimitedFormatter treeFormatter = new TabDelimitedFormatter(
//				new PrintWriter(new FileOutputStream(new File(logTreeName))));
//
//		loggers[2] = new TreeLogger(treeModel, branchRateModel, null, null,
//				treeFormatter, logInterval, true, true, true, null, null);
//
////		 log spectrum??
////		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
////		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
////		ShortReadLikelihood trueSrp = new ShortReadLikelihood(HaplotypeModelUtils.factory(shortReads, trueAlignment));
////		System.err.println("\'trueShortReadLikelihood\': "+trueSrp.getLogLikelihood());
//		loggers[3] = new SpectrumLogger(spectrumModel, trueAlignment, logHaplotypeName, logInterval*totalSamples/10);
//		
		
		
		// MCMC
//		MCMCOptions options = MCMCSetupHelper.setMCMCOptions(logInterval, totalSamples);
		
//		MCMC mcmc = new MCMC("mcmc1");
//		mcmc.setShowOperatorAnalysis(true);
//		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
//		
//		mcmc.init(options, posterior, schedule, loggers);
////		mcmc.run();
//		System.out.println(mcmc.getTimer().toString());

		
		//MCMCMC
		double DEFAULT_DELTA = 1.0;
	    int DEFAULT_SWAP_CHAIN_EVERY = 100;
	    int swapChainsEvery = 100;
		double delta = DEFAULT_DELTA;

		int chainCount = 2;
	    double[] chainTemperatures = new double[chainCount];
	    chainTemperatures[0] = 1.0;

		for (int i = 1; i < chainCount; i++) {
			chainTemperatures[i] = 1.0 / (1.0 + (delta * i));
		} 
		
		//

         MCMC[] chains = new MCMC[chainCount];
         MCMCMCOptions optionsMC3 = new MCMCMCOptions(chainTemperatures, swapChainsEvery);
         
         for (int i = 0; i < chainCount; i++) {
             // parse the file once for each hot chain
        	 chains[i] = new MCMC("mcmc"+i);
        	 chains[i].setShowOperatorAnalysis(true);
        	 chains[i].setOperatorAnalysisFile(new File(operatorAnalysisFile));
        	 
////////////
        	 SpectrumAlignmentModel spectrumModel ;
     		ShortReadsSpectrumLikelihood srpLikelihood ;
     		
//     		boolean redo = true;
//     		int c = 0;
     		if (randomSpectrum) {
//     			do {
     				spectrumModel = new SpectrumAlignmentModel(spectrumLength, noOfRecoveredHaplotype, randomSpectrumType);
     				srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap, distTypeCode);
//     				redo = (srpLikelihood.getLogLikelihood() == Double.NEGATIVE_INFINITY);
//     				c++;
//     				if(c==100){
//     					System.err.println("After 100 try, ShortreadLikelihood= "+ srpLikelihood.getLogLikelihood());
//     					System.exit(-1);
//     				}
//     			} while (redo);
     		}		
     		else{
     			spectrumModel = dataImporter.importPartialSpectrumFile(partialSpectrumName );
     			srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap, distTypeCode);
     		}
//     		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, noOfRecoveredHaplotype, randomSpectrumType);
//     		srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap, distTypeCode);
     		// coalescent
     		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);
     		
     		// Random treeModel
     		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
     		TreeModel treeModel;
     		if(randomTree){
     			treeModel = MCMCSetupHelper.setupRandomTreeModel(popModel, spectrumModel, Units.Type.YEARS);
     		}
     		else{
     			Tree partialPhylogeny = dataImporter.importTree(partialTreeName);
     			treeModel = new TreeModel(TreeModel.TREE_MODEL, partialPhylogeny, false);
     		}
     		// Coalescent likelihood
     		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), popModel);
     		coalescent.setId("coalescent");

     		
     		// Simulate  treeLikelihood
     		HashMap<String, Object> parameterList = MCMCSetupHelperSpectrum.setupSpectrumTreeLikelihoodSpectrumModel(treeModel, spectrumModel);
     		Parameter kappa = (Parameter) parameterList.get("kappa");
     		Parameter freqs = (Parameter) parameterList.get("freqs");
     		Parameter rate = (Parameter) parameterList.get("rate");
     		StrictClockBranchRates branchRateModel = (StrictClockBranchRates) parameterList.get("branchRateModel");
     		SpectrumTreeLikelihood treeLikelihood = (SpectrumTreeLikelihood) parameterList.get("treeLikelihood");
     		
     		
     		
     		// CompoundLikelihood
     		HashMap<String, Likelihood> compoundlikelihoods = MCMCSetupHelper.setupCompoundLikelihood(
     				popSize, kappa, coalescent, treeLikelihood, srpLikelihood);
     		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
     		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
     		Likelihood shortReadLikelihood = compoundlikelihoods.get(ShortReadsSpectrumLikelihood.SHORT_READ_LIKELIHOOD);
     		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);

     		// Operators
     		OperatorSchedule schedule = new SimpleOperatorSchedule();
//     		ArrayList<MCMCOperator> defalutOperatorsList = 
     		MCMCSetupHelperSpectrum.defalutSpectrumOperators(schedule, spectrumModel, freqs, popSize, kappa);
     		MCMCSetupHelperSpectrum.defalutTreeOperators(schedule, treeModel);
     		
     		
//     		MCMCOperator operator;
//     		operator = new RJTreeOperator(spectrumModel, treeModel);
//     		operator.setWeight(100);
//     		schedule.addOperator(operator);
     		
     		Parameter rootHeight = treeModel.getRootHeightParameter();
     		rootHeight.setId("rootHeight");
        	 MCMCOptions options = MCMCSetupHelper.setMCMCOptions(logInterval, totalSamples);
			// MCLogger
			MCLogger[] loggers = new MCLogger[2];
			// log tracer
			loggers[0] = new MCLogger(logTracerName+i, logInterval, false, 0);
			MCMCSetupHelper.addToLogger(loggers[0], posterior, prior,
					likelihood, shortReadLikelihood, rootHeight, rate, popSize,
					kappa, coalescent, freqs);
			// System.out
			loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out),
					logInterval, true, logInterval * 2);
			MCMCSetupHelper.addToLogger(loggers[1],
					// freqs
					posterior, prior, likelihood, shortReadLikelihood, rate,
					popSize, kappa, coalescent, rootHeight);
        	 chains[i].init(options, posterior, schedule, loggers);

         }
	            
	            
		  MCMCMC mc3 = new MCMCMC(chains, optionsMC3);
		  Thread thread = new Thread(mc3);
//          thread.start();
		  mc3.run();
          System.out.println(mc3.getColdChain());
		  System.out.println(mc3.getCurrentState());
          System.out.println(mc3.getChainLength());
          System.out.println(mc3.getTimer().toString());
	}


}


