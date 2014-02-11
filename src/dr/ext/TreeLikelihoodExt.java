package dr.ext;

import javax.swing.text.TabableView;

import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.GeneralLikelihoodCore;
import dr.evomodel.treelikelihood.NativeNucleotideLikelihoodCore;
import dr.evomodel.treelikelihood.NucleotideLikelihoodCore;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Model;

public class TreeLikelihoodExt extends TreeLikelihood {

	private HaplotypeModel haplotypeModel;
	private SitePatternsExt sitePatternExt;

	public TreeLikelihoodExt(HaplotypeModel haplotypeModel, TreeModel treeModel,
			SiteModel siteModel, BranchRateModel branchRateModel,
			TipStatesModel tipStatesModel, boolean useAmbiguities,
			boolean allowMissingTaxa, boolean storePartials,
			boolean forceJavaCore, boolean forceRescaling) {
		
		
		super(new SitePatternsExt (haplotypeModel, null, 0, -1, 1, true),
				treeModel, siteModel, branchRateModel,
				tipStatesModel, useAmbiguities, allowMissingTaxa,
				storePartials, forceJavaCore, forceRescaling);
		
		this.sitePatternExt = (SitePatternsExt) getPatternList(); 
		this.haplotypeModel = haplotypeModel;
		addModel(this.haplotypeModel);
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 6348244745369755426L;


    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    	if (model == haplotypeModel){
//    		System.out.println("GOOD here");
    		sitePatternExt.updateAlignment(haplotypeModel);
    		updatePatternListExt(sitePatternExt);
    		likelihoodKnown = false;
    	}
    	else{
    		super.handleModelChangedEvent(model, object, index);
    	}
    }
    protected void restoreState() {

//        if (storePartials) {
//            likelihoodCore.restoreState();
//        } else {
//            updateAllNodes();
//        }
    	sitePatternExt.updateAlignment(haplotypeModel);
		updatePatternListExt(sitePatternExt);
        super.restoreState();

    }


	
	public void updatePatternListExt(PatternList patternList) {
		

        this.patternList = patternList;
//        this.dataType = patternList.getDataType();
        patternCount = patternList.getPatternCount();
//        stateCount = dataType.getStateCount();

        patternWeights = patternList.getPatternWeights();

//        this.treeModel = treeModel;
//        addModel(treeModel);

//        nodeCount = treeModel.getNodeCount();

//        updateNode = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }

        likelihoodKnown = false;
		
        patternLogLikelihoods = new double[patternCount];
//        resetRootPartials();//TODO This will FAIL!! need to change BEAST code to make this method work
//        rootPartials = new double[patternCount * stateCount];

//	}
//
//	public void tTreeLikelihood(PatternList patternList, TreeModel treeModel,
//			SiteModel siteModel, BranchRateModel branchRateModel,
//			TipStatesModel tipStatesModel, 
//			boolean useAmbiguities,
//			boolean allowMissingTaxa, boolean storePartials,
//			boolean forceJavaCore, boolean forceRescaling) {

//		super(TreeLikelihoodParser.TREE_LIKELIHOOD, patternList, treeModel);
		 

//		        super(TreeLikelihoodParser.TREE_LIKELIHOOD, patternList, treeModel);

//		        this.storePartials = storePartials;
		boolean useAmbiguities = false;
		boolean allowMissingTaxa = false;
//		boolean storePartials,
//		boolean forceJavaCore, 
		boolean forceRescaling = false;

		        try {
//		            this.siteModel = siteModel;
//		            addModel(siteModel);

//		            this.frequencyModel = siteModel.getFrequencyModel();
//		            addModel(frequencyModel);

//		            this.tipStatesModel = tipStatesModel;

//		            integrateAcrossCategories = siteModel.integrateAcrossCategories();

		            categoryCount = siteModel.getCategoryCount();

//		            final Logger logger = Logger.getLogger("dr.evomodel");
//		            String coreName = "Java general";
//		            boolean forceJavaCore = false;
//		            if (integrateAcrossCategories) {
//
//		                final DataType dataType = patternList.getDataType();
//
//		                if (dataType instanceof dr.evolution.datatype.Nucleotides) {
//
//		                    if (!forceJavaCore && NativeNucleotideLikelihoodCore.isAvailable()) {
//		                        coreName = "native nucleotide";
//		                        likelihoodCore = new NativeNucleotideLikelihoodCore();
//		                    } else {
//		                        coreName = "Java nucleotide";
//		                        likelihoodCore = new NucleotideLikelihoodCore();
//		                    }
//
//		                } 
//		            } else {
//		                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
//		            }

		            probabilities = new double[stateCount * stateCount];

		            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);

		            int extNodeCount = treeModel.getExternalNodeCount();
		            int intNodeCount = treeModel.getInternalNodeCount();
		           
		            {
		                for (int i = 0; i < extNodeCount; i++) {
		                    // Find the id of tip i in the patternList
		                    String id = treeModel.getTaxonId(i);
		                    int index = patternList.getTaxonIndex(id);

		                    if (index == -1) {
		                        if (!allowMissingTaxa) {
		                            throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
		                                    ", is not found in patternList, " + patternList.getId());
		                        }
		                        if (useAmbiguities) {
		                            setMissingPartials(likelihoodCore, i);
		                        } else {
		                            setMissingStates(likelihoodCore, i);
		                        }
		                    } else {
		                        if (useAmbiguities) {
		                            setPartials(likelihoodCore, patternList, categoryCount, index, i);
		                        } else {
		                            setStates(likelihoodCore, patternList, index, i);
		                        }
		                    }
		                }
		            }
		            for (int i = 0; i < intNodeCount; i++) {
		                likelihoodCore.createNodePartials(extNodeCount + i);
		            }



		        } catch (TaxonList.MissingTaxonException mte) {
		            throw new RuntimeException(mte.toString());
		        }

//		        addStatistic(new SiteLikelihoodsStatistic());
//		        System.out.println(getStatisticCount());
//		        System.out.println(getStatistic(0).getDimension()+"\t"+ getStatistic(0).getStatisticValue(10));
	}

	public void updatePatternList(SitePatternsExt patterns, HaplotypeModel haplotypeModel) {
//		TODO: more test required
//		Alignment alignment = haplotypeModel.getAlignment();
        
		patterns.updateAlignment(haplotypeModel);
		updatePatternListExt(patterns);
	}

	public void updatePatternList(HaplotypeModel haplotypeModel) {
//		TODO: more test required
//		sitePatternExt.updateAlignment(haplotypeModel);
		updatePatternListExt(haplotypeModel);
	}


}
