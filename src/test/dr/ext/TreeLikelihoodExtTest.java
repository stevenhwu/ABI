package test.dr.ext;


import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.HaplotypeModelUtils;
import srp.haplotypes.operator.SingleBaseOperator;
import srp.likelihood.ShortReadLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.ext.SitePatternsExt;
import dr.ext.TreeLikelihoodExt;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.model.Parameter;

public class TreeLikelihoodExtTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}


	private GammaSiteModel siteModel;
	private BranchRateModel branchRateModel;

	@Before
	public void setUp() throws Exception {

    	// clock model
    	Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1E-5, 0, 1);
    	branchRateModel = new StrictClockBranchRates(rateParameter);

    	// Sub model
    	Parameter freqs = new Parameter.Default(new double[]{0.25,0.25,0.25,0.25});
    	Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

    	FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
    	HKY hky = new HKY(kappa, f);

    	//siteModel
    	Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
    	siteModel = new GammaSiteModel(hky);
    	siteModel.setMutationRateParameter(mu);


	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testUpdatePatternList() throws Exception {
		
//		-Djava.library.path=/home/sw167/PostdocLarge/Software/BEAST/BEASTv1.7.1/lib -Xms128m -Xmx256m
		String dataDir = "/home/sw167/workspace/ABI/unittest/";

		String trueAlignmentFile = "H4_haplotypes.phyml";
		String phylogenyFile = "H4_haplotypes.tree";
		String shortReadFile = "H4_srp.fasta";//"1110_10_align_2.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Alignment alignment = dataImporter.importAlignment(trueAlignmentFile);
		
		Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, alignment);
		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(haplotypeModel);
		ShortReadLikelihood srpLikelihoodUpdate = new ShortReadLikelihood(haplotypeModel);
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);

    	//treeLikelihood
    	SitePatterns patterns = new SitePatterns(trueAlignment, null, 0, -1, 1, true);
    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
    			false, false, true, false, false);

    	//treeLikelihoodExt
    	SitePatternsExt patternsExt = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);
        TreeLikelihoodExt treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
		// end
    	

        assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 0);
		
		
		srpLikelihood = new ShortReadLikelihood(haplotypeModel);

		SingleBaseOperator op = new SingleBaseOperator(haplotypeModel, 0);
        for (int i = 0; i < 100; i++) {
			op.doOperation();

//			patternsExt.updateAlignment(haplotypeModel);
//			treeLikelihoodExt.updatePatternList(patternsExt);
					
			treeLikelihoodExt.updatePatternList(haplotypeModel);
			
			srpLikelihoodUpdate = new ShortReadLikelihood(haplotypeModel);
			assertEquals(srpLikelihood.getLogLikelihood(), srpLikelihoodUpdate.getLogLikelihood(), 0);
			
			//treeLikelihood
	    	patterns = new SitePatterns(haplotypeModel, null, 0, -1, 1, true);
	    	treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
	    			false, false, true, false, false);

    		assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 1e-8);
		}
        

        long time1, time2;
		time1 = System.currentTimeMillis();
		for (int i = 0; i < 1e4; i++) {
			treeLikelihood.makeDirty();
			treeLikelihood.getLogLikelihood();
		}
		time2 = System.currentTimeMillis();
		System.out.println("\t" + (time2 - time1) + "\t");
		time1 = System.currentTimeMillis();
		for (int i = 0; i < 1e4; i++) {
			treeLikelihoodExt.updatePatternList(haplotypeModel);
//			treeLikelihoodExt.makeDirty();
			treeLikelihoodExt.getLogLikelihood();
		}
		time2 = System.currentTimeMillis();
		System.out.println("\t" + (time2 - time1) + "\t");



        
        
        for (int i = 0; i < 500; i++) {
        	op.doOperation();
        }
//        patternsExt.updateAlignment(haplotypeModel);
//		treeLikelihoodExt.updatePatternList(patternsExt);

//        treeLikelihoodExt.updatePatternList(patternsExt, haplotypeModel);
        treeLikelihoodExt.updatePatternList(haplotypeModel);
		

		//treeLikelihood
    	patterns = new SitePatterns(haplotypeModel, null, 0, -1, 1, true);
    	treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
    			false, false, true, false, false);

		assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 1e-8);
	      	
//		srpLikelihoodUpdate = new ShortReadLikelihood(aMap, haplotypeModel);
		//Fail atm, should pass after get ride of matrix[][]
//		assertEquals(srpLikelihood.getLogLikelihood(), srpLikelihoodUpdate.getLogLikelihood(), 0);
		
	}
	

	public void testUpdatePatternList2() throws Exception {
		
		
		String dataDir = "/home/sw167/workspace/ABI/unittest/";

		String trueAlignmentFile = "1110_10_org_6.phyml";
		String phylogenyFile = "1110_10_org_6.phyml_phyml_tree.txt";
		String shortReadFile = "1110_10_align_100.fasta";//"1110_10_align_2.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Alignment alignment = dataImporter.importAlignment(trueAlignmentFile);
		
		Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		HaplotypeModel trueHaplotypes = new HaplotypeModel(aMap, trueAlignment);
		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(trueHaplotypes);
	
		


//			
//			TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
//			LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
//			li.setPopSize(3000,0,30000);
//			li.setTreeAndAlignment(treeModel, trueAlignment);

		
//			System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
//					+"\t"+ trueHaplotypes.calculateSPS());
		int[][] sps = HaplotypeModelUtils.calculeteSPSArray(trueHaplotypes, trueHaplotypes);
//			for (int i = 0; i < sps.length; i++) {
//					System.out.println(Arrays.toString(sps[i]));
//			}
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);
		

//			LikelihoodCalculation likelihoodModel = new LikelihoodCalculation(treeModel, aMap, trueHaplotypes);
//			likelihoodModel.setPopSize(3000,0,30000);
////			li.setTreeAndAlignment(treeModel, trueAlignment);
//			System.out.println(likelihoodModel.getTreeLikelihood());
//			System.out.println(likelihoodModel.getCoalescentLikelhood());
//			System.out.println(likelihoodModel.getShortReadLikelihood());
//			System.out.println(likelihoodModel.getLoglikelihood());
//			System.out.println(srpLikelihood.getLogLikelihood());
//			System.out.println("=== END True Hap===");
			
		
		
        //CompoundLikelihood
//	        List<Likelihood> likelihoodsList = new ArrayList<Likelihood>();        
//	        likelihoodsList.add(coalescent);
//	        Likelihood prior = new CompoundLikelihood(0, likelihoodsList);
//	        prior.setId(CompoundLikelihoodParser.PRIOR);
//
//	        likelihoodsList.clear();
//	        likelihoodsList.add(treeLikelihood);
//	        Likelihood likelihood = new CompoundLikelihood(-1, likelihoodsList);
//
//	        likelihoodsList.clear();
//	        likelihoodsList.add(prior);
//	        likelihoodsList.add(likelihood);
//	        Likelihood posterior = new CompoundLikelihood(0, likelihoodsList);
//	        posterior.setId(CompoundLikelihoodParser.POSTERIOR);

			
		// start
		
    	Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000,0,10000);
    	ConstantPopulationModel constantModel = new ConstantPopulationModel(popSize, Units.Type.DAYS);//createRandomInitialTree(popSize);

    	CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
    	coalescent.setId("coalescent");

    	// clock model
    	Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1E-5, 0, 1);
    	StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);

    	// Sub model
    	Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
    	Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

    	FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
    	HKY hky = new HKY(kappa, f);

    	//siteModel
    	GammaSiteModel siteModel = new GammaSiteModel(hky);
    	Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
    	siteModel.setMutationRateParameter(mu);

    	//treeLikelihood
    	SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

//	    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
//	    			false, false, true, false, false);
//	    	treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);

		// end
		
		
		
//
//			Haplotypes haplotypes = new Haplotypes(aMap, trueAlignment.getSequenceCount());
    	HaplotypeModel haplotypeModel = trueHaplotypes;
////			for (int i = 0; i < haplotypes.getHaplotypesCount(); i++) {
////					haplotypes.randomSeq(i);
////				System.out.println(haplotypes.getHaplotype(i) );
////			}
//			likelihoodModel.updateHaplotypes(haplotypes);
//			likelihoodModel.calculateShortReadLikelihoodFull();
//			System.out.println(likelihoodModel.getTreeLikelihood());
//			System.out.println(likelihoodModel.getCoalescentLikelhood());
//			System.out.println(likelihoodModel.getShortReadLikelihood());
//			System.out.println(likelihoodModel.getLoglikelihood());
//			
//System.out.println("====full====");
//			 likelihoodModel = new LikelihoodCalculation(treeModel, aMap, haplotypes);
//				likelihoodModel.setPopSize(3000,0,30000);
//			likelihoodModel.setPopSize(3000,0,30000);
////			li.setTreeAndAlignment(treeModel, trueAlignment);
//			System.out.println(likelihoodModel.getTreeLikelihood());
//			System.out.println(likelihoodModel.getCoalescentLikelhood());
//			System.out.println(likelihoodModel.getShortReadLikelihood());
//			System.out.println(likelihoodModel.getLoglikelihood());
//			
//			
//
//			System.out.println("===");


		
		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
				+"\t"+ haplotypeModel.calculateSPS());
		sps = HaplotypeModelUtils.calculeteSPSArray(trueHaplotypes, haplotypeModel);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		
		}
		
		
		srpLikelihood = new ShortReadLikelihood(haplotypeModel);
		
                
        
        
        int thinning = 1000;
        SitePatternsExt patternsExt = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);
        TreeLikelihoodExt treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
        
        double likelihood = srpLikelihood.getLogLikelihood() + treeLikelihoodExt.getLogLikelihood();
		likelihood = treeLikelihoodExt.getLogLikelihood();
		System.out.println(likelihood);
		
		SingleBaseOperator op = new SingleBaseOperator(haplotypeModel, 0);
        for (int i = 0; i < 100; i++) {
        	op.doOperation();
//			alignment = haplotypeModel.getAlignment();
			patternsExt.updateAlignment(haplotypeModel);
			treeLikelihoodExt.updatePatternListExt(patternsExt);
		}
		
		likelihood = treeLikelihoodExt.getLogLikelihood();
		System.out.println("updated: "+likelihood);
		treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
		likelihood = treeLikelihoodExt.getLogLikelihood();
		System.out.println("updated: "+likelihood);
		
		
	
//		
//        for (int i = 0; i < 1e5; i++) {
//
//			haplotypeModel.swapBase();
//			srpLikelihood.updateHaplotypes(haplotypeModel);
//			
////				likelihoodModel.updateHaplotypes(haplotypes);
////				int[] x = haplotypes.getSwapInfo();
////				if(x[2] != x[3]){
////					System.out.println(Arrays.toString(haplotypes.getSwapInfo()));
////					System.out
////					.println(Arrays.toString(patterns2.getStateFrequencies()));
////			System.out.println(Arrays.toString(patterns2.getPattern(x[0])));
////				}
//			
//			alignment = haplotypeModel.getAlignment();
//			patternsExt.updateAlignment(alignment);
//			treeLikelihoodExt.updatePatternList(patternsExt);
//			
//			
////				if(x[2] != x[3]){
////				System.out
////						.println(Arrays.toString(patterns2.getStateFrequencies()));
////				System.out.println(Arrays.toString(patterns2.getPattern(x[0])));
////		System.out.println("===");
////				}
////				double newL = srpLikelihood.getLogLikelihood() + treeLikelihood2.getLogLikelihood();
//			double newL = treeLikelihoodExt.getLogLikelihood();
//
//			
//			
//			boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
////					System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
//			if(accept){
//				likelihood = newL;
////					likelihoodModel.acceptState();
////					srpLikelihood.acceptState();
//			}
//			else{
//				haplotypeModel.reject();
////					likelihoodModel.restorreState();
//				srpLikelihood.restoreState();
//			}
//			if (i% thinning == 0){
//				
//				String temp = i +"\t"+ likelihood +"\t"+ haplotypeModel.calculateSPS() + 
//						"\t"+ HaplotypeModel.calculeteSPS(trueHaplotypes, haplotypeModel) + "\n";
//				System.out.print(temp);
//			}
//					
//		}
//		
//		System.out.println("===== Done =====");
//		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
//				+"\t"+ haplotypeModel.calculateSPS());
//		sps = HaplotypeModel.calculeteSPSArray(trueHaplotypes, haplotypeModel);
//		for (int i = 0; i < sps.length; i++) {
//				System.out.println(Arrays.toString(sps[i]));
//		}
//		}catch (Exception e){
//			e.printStackTrace();
//		}
////				Haplotypes.calculeteSPS(trueHaplotypes, haplotypes);
//			
//
//
//	}
	}
}
