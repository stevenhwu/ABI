package srp.spectrum.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.google.common.primitives.Doubles;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.likelihood.LikelihoodScaler;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.SpectrumOperationRecord;
import dr.app.beauti.util.NumberUtil;
import dr.evolution.datatype.DataType;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

/*

Three methods to go through multiple sites
Method 1: Go through them, might be faster for small srp size/small number of bases change

	for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
		mapToSrp = aMap.getMapToSrp(k);
		for (int i : mapToSrp) {
			for (int r = 0; r < TWO; r++) {
				j = twoSpectrums[r];
				storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
				storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
				count2++;
			}
		}
	}


Method 2:
get all unique sites with HashSet

	Set<Integer> allSrpPos = new HashSet<Integer>();
	allSrpPos.clear();
	for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
		mapToSrp = aMap.getMapToSrp(i);
		allSrpPos.addAll(mapToSrp);
	}
	for (int i : allSrpPos) {
		for (int r = 0; r < TWO; r++) {
			for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
				j = twoSpectrums[r];
				if(allLogLikelihood[i][j][k] != 0){
					storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
					storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
				}
			}
		}
	}
Method 3: boolean index array
				
	int srpCount = aMap.getSrpCount();
	boolean[] srpSwitch = new boolean[srpCount];
	Arrays.fill(srpSwitch, false);
	for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
		mapToSrp = aMap.getMapToSrp(k);
		for (int i : mapToSrp) {
			srpSwitch[i] = true;
		}
	}
				
	for (int i = 0; i < srpSwitch.length; i++) {
		if(srpSwitch[i]){
			for (int r = 0; r < TWO; r++) {
				j = twoSpectrums[r];
				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
					if(allLogLikelihood[i][j][k] != 0){
						storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
						storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
					}
				}
			}
		}
	}
				

 */


public class ShortReadsSpectrumLikelihood  extends AbstractModelLikelihood {

	/**
	 * 
	 */
	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadSpectrumLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double NOT_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
//	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	private static final int TWO = 2;

	private static final int GAP_STATE = 17;  
	private final int AMBIGUOUS_STATE_COUNT;
	private final int STATE_COUNT;
	
	protected boolean likelihoodKnown;
	
	private int spectrumLength;
	private int spectrumCount;
	private int srpCount;

	private double logLikelihood;
	private double storedLogLikelihood;

	private double[] eachSrpLogLikelihood;
	private double[] storedEachSrpLogLikelihood;

	private double[] sumScaledSrpLogLikelihood;
	private double[] storedSumSrpLogLikelihood;

	
	private double[][] spectrumLogLikelihood;
	private double[][] storedSpectrumLogLikelihood;
	
//	private double[][] spectrumScaledLogLikelihood;
//	private double[][] storedSpectrumScaledLogLikelihood;
	
	private double[][][] allLogLikelihood;
	private double[][][] storedAllLogLikelihood;


	private AlignmentMapping aMap;
	private LikelihoodScaler liS;
	
	private SpectrumAlignmentModel spectrumModel;
	private DataType dataType;
	
	
	private boolean[] srpSwitch;
	private Set<Integer> allSrpPos;

	private boolean useSrpArray = true; //true: sue booleanArray, false: use HashSet
//	private boolean useSrpArray = false; //true: sue booleanArray, false: use HashSet

	
		
	public ShortReadsSpectrumLikelihood(SpectrumAlignmentModel spectrumModel){
		super(SHORT_READ_LIKELIHOOD);
		this.spectrumModel = spectrumModel;
		this.aMap = this.spectrumModel.getAlignmentMapping();
		
		this.dataType = this.spectrumModel.getDataType();
		STATE_COUNT = dataType.getStateCount();//4
		AMBIGUOUS_STATE_COUNT = dataType.getAmbiguousStateCount();//18

		likelihoodKnown = false;
		
		
		addModel(this.spectrumModel);
		
		preprocessLikelihoodAlignmentMap();
		getLogLikelihood();
		storeEverything();
		
	}
	

	private void preprocessLikelihoodAlignmentMap() {
//		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = aMap.getSrpCount();
		spectrumCount = spectrumModel.getSpectrumCount();
		spectrumLength = spectrumModel.getSpectrumLength();
		
		this.srpSwitch = new boolean[srpCount];
		this.allSrpPos = new HashSet<Integer>();
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;
		
		allLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];
		storedAllLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];

		spectrumLogLikelihood = new double[srpCount][spectrumCount];
		storedSpectrumLogLikelihood = new double[srpCount][spectrumCount];
		
//		spectrumScaledLogLikelihood = new double[srpCount][spectrumCount];
//		storedSpectrumScaledLogLikelihood = new double[srpCount][spectrumCount];
		
		sumScaledSrpLogLikelihood = new double[srpCount];
		storedSumSrpLogLikelihood = new double[srpCount];

		
		eachSrpLogLikelihood = new double[srpCount];
		storedEachSrpLogLikelihood = new double[srpCount];

		
		
	}


    
	@Override
	public double getLogLikelihood(){


//		long time1 = System.currentTimeMillis();
		
		
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
        
        return logLikelihood;


		
		

	}
	
	protected double calculateLogLikelihood() {
//		SpectrumOperationRecord operationReocrd = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumModel.getSpectrumOperation();
		double logLikelihood = Double.NEGATIVE_INFINITY;
//System.err.println("calculateLikelihood\t"+operation);

		switch (operation) {
			case NONE:
			case FULL:
				if(DEBUG){
					System.out.println("Calculate ShortReadLikelihood:\t"+operation);
				}
				logLikelihood = calculateSrpLikelihoodFull();
				break;
//			case SWAPMULTI:
//				logLikelihood = calculateSrpLikelihoodMultiBasesSwap();
//				break;

			case DELTA_SINGLE:
			case SWAP_SINGLE:
//				System.out.println("single");
//				logLikelihood = calculateSrpLikelihoodFull();				
				logLikelihood = calculateSrpLikelihoodSingle();
				break;
			case DELTA_COLUMN:
			case SWAP_COLUMN:
				logLikelihood = calculateSrpLikelihoodColumn();
//				logLikelihood = calculateSrpLikelihoodFull();
				break;
			case DELTA_MULTI:
			case SWAP_MULTI:
//				logLikelihood = calculateSrpLikelihoodFull();
				logLikelihood = calculateSrpLikelihoodMulti();
				
				break;
			case SWAP_SUBCOLUMN:
				logLikelihood = calculateSrpLikelihoodSwapSubColumn();
				break;
			case RECOMBINATION:
				logLikelihood = calculateSrpLikelihoodRecombination();
				break;

//			case PASS:
//				logLikelihood = storedLogLikelihood;
//				break;
			default:
				throw new IllegalArgumentException("Unknown operation type: "+operation);
	
			}
//	    double logLikelihood = calculateShoreReadLikelihood4();
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
//		storeState();
//		System.out.println("likelihood\t"+ logLikelihood);
//		if( (logLikelihood-this.logLikelihood)==0 ){
//		System.err.println("Delta: "+ (logLikelihood-this.logLikelihood) +"\t"+ logLikelihood +"\t"+ this.logLikelihood +"\t"+ operation);
////		System.out.println("Delta: "+ (logLikelihood-this.logLikelihood) +"\t"+ logLikelihood +"\t"+ this.logLikelihood +"\t"+ operation);
//		}
//		else{
//			System.out.println("Delta: "+ (logLikelihood-this.logLikelihood) +"\t"+ logLikelihood +"\t"+ this.logLikelihood +"\t"+ operation);
//		}
		return logLikelihood;
	}
	
	private double calculateSrpLikelihoodFull() {

//		System.out.println("calculateSrpLikelihoodFull");
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex(); 
//		int siteIndex = record.getAllSiteIndexs()[0];
//		System.out.println("In fullCaluclation\t"+spectrumIndex +"\t"+ siteIndex +"\t"+ 
//				Arrays.toString(spectrumModel.getSpectrum(spectrumIndex)
//				.getFrequencies(siteIndex)));
		for (int i = 0; i < srpCount; i++) {

//			String srp = aMap.getSrpFragment(i);
			String fullSrp = aMap.getSrpFull(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			
			double logLikelihood;
			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {

				Spectrum spectrum = spectrumModel.getSpectrum(j);
				spectrumLogLikelihood[i][j] = 0;
				for (int k = start; k < end; k++) {
					double[] frequencies = spectrum.getFrequenciesAt(k);
					char srpChar = fullSrp.charAt(k);
					int state = dataType.getState(srpChar);
					if(state<STATE_COUNT){
						double likelihood = frequencies[state] * NOT_ERROR_RATE
								+ (1 - frequencies[state]) * ERROR_RATE;
						logLikelihood = Math.log(likelihood);
					}
					else{
						logLikelihood = LOG_ERROR_RATE;
					}
//					if(spectrumIndex>0 && (allLogLikelihood[i][j][k] != logLikelihood)){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ logLikelihood + " "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
					allLogLikelihood[i][j][k] = logLikelihood;
					spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
					
				}
//				spectrumScaledLogLikelihood[i][j] = liS.scale(spectrumLogLikelihood[i][j]);
//				liS.addScaledLogProb(spectrumScaledLogLikelihood[i][j]);
				liS.addScaleLogProb(spectrumLogLikelihood[i][j]);
				

			}	
			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
		}
		double logLikelihood = liS.sumLogLikelihood(sumScaledSrpLogLikelihood);
		double logLikelihood2 = StatUtils.sum(eachSrpLogLikelihood);
		
//for (int i = 0; i < srpCount; i++) {
//	System.out.println(Arrays.toString(allDists[i]));
//}
//System.out.println("==");

//		double logLikelihood2 = calculateSrpLikelihoodFull();
//		if(logLikelihood != logLikelihood2){
//			System.out.println(logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
//		
//		System.out.println(logLikelihood);
		return logLikelihood;
	}

	
	private double calculateSrpLikelihoodSingle() {


		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int j = record.getSpectrumIndex(); 
		int k = record.getAllSiteIndexs()[0];
//
		Spectrum spectrum = spectrumModel.getSpectrum(j);
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
		 
		double[] statesLogLikelihood = calculateStatesLogLikelihood(spectrum, k);
		double totalLogLikelihood = this.logLikelihood;
//		System.out.println(j +"\t"+ k);
//		System.out.println(totalLogLikelihood);
		for (int i : mapToSrp) {
			
			String fullSrp = aMap.getSrpFull(i);
			totalLogLikelihood = updateLikelihoodAtIJK(i, j, k, fullSrp, statesLogLikelihood, totalLogLikelihood);
//			updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);
//			updateEachSrpAt(i);
		}
//		System.out.println(totalLogLikelihood);
//		System.out.println(Arrays.toString(eachSrpLogLikelihood));
//		for (int i : mapToSrp) {
//			
//			String fullSrp = aMap.getSrpFull(i);
//			updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);
//			updateEachSrpAt(i);
//		}
//		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
//		System.out.println(logLikelihood);
//		System.out.println(Arrays.toString(eachSrpLogLikelihood));
//		System.out.println();
		return totalLogLikelihood;
	}
	

	@Deprecated
	private double calculateSrpLikelihoodSingle_Old() {


		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int j = record.getSpectrumIndex(); 
		int k = record.getAllSiteIndexs()[0];

		Spectrum spectrum = spectrumModel.getSpectrum(j);
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
		
		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);
			updateEachSrpAt(i);
		}

		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
	
		return logLikelihood;
	}

	private double calculateSrpLikelihoodMulti() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();

		int[] siteIndexs = record.getAllSiteIndexs();
		int j= record.getSpectrumIndex(); 
		Spectrum spectrum = spectrumModel.getSpectrum(j);

		double[][] allStateLogLikelihood = new double[siteIndexs.length][]; 
		for (int s = 0; s < siteIndexs.length; s++) {
			int k = siteIndexs[s];
			allStateLogLikelihood[s] = calculateStatesLogLikelihood(spectrum, k);
		}
		
		if(useSrpArray){
			for (int s : siteIndexs) {
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
				for (int i : mapToSrp) {
					srpSwitch[i] = true;
				}
			}
			
			for (int i = 0; i < srpSwitch.length; i++) {
				if(srpSwitch[i]){
					String fullSrp = aMap.getSrpFull(i);
				
					for (int s = 0; s < siteIndexs.length; s++) {
						int k = siteIndexs[s];
						if(allLogLikelihood[i][j][k]!=0){
//							updateLikelihoodAtIJK(i, j, k,spectrum, fullSrp);
							updateLikelihoodAtIJK(i, j, k, fullSrp, allStateLogLikelihood[s], logLikelihood);
						}
					}
					updateEachSrpAt(i);
				}

			}
		}
		else{
			allSrpPos.clear();
			for (int s : siteIndexs) {
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
				allSrpPos.addAll(mapToSrp);
			}
			for (int i : allSrpPos) {
				String fullSrp = aMap.getSrpFull(i);
			
				for (int s = 0; s < siteIndexs.length; s++) {
					int k = siteIndexs[s];
					if(allLogLikelihood[i][j][k]!=0){
//							updateLikelihoodAtIJK(i, j, k,spectrum, fullSrp);
						updateLikelihoodAtIJK(i, j, k, fullSrp, allStateLogLikelihood[s], logLikelihood);
					}
				}
				updateEachSrpAt(i);
			}
		}
		
		

		
		
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);

		return logLikelihood;

	}

	@Deprecated
	private double calculateSrpLikelihoodMulti2() {
	
			SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
	
			int[] siteIndexs = record.getAllSiteIndexs();
			int j= record.getSpectrumIndex(); 
			Spectrum spectrum = spectrumModel.getSpectrum(j);
			
	//		Arrays.fill(srpSwitch, true);
			for (int s : siteIndexs) {
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
				for (int i : mapToSrp) {
					srpSwitch[i] = true;
				}
			}
	
		
	//	Set<Integer> allSrpPos = new HashSet<Integer>();
	//	allSrpPos.clear();
	//	for (int s : siteIndexs) {
	//
	//		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
	//		allSrpPos.addAll(mapToSrp);
	//	}
	//	System.out.println("Total:\t"+allSrpPos.size());
	
	
			for (int i = 0; i < srpSwitch.length; i++) {
				
				if(srpSwitch[i]){
					String fullSrp = aMap.getSrpFull(i);
				
					for (int s = 0; s < siteIndexs.length; s++) {
						int k = siteIndexs[s];
						if(allLogLikelihood[i][j][k]!=0){
							updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);
						}
					}
					updateEachSrpAt(i);
				}
	
			}
			double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
	
			return logLikelihood;
	
		}

	@Deprecated
	private double calculateSrpLikelihoodMulti3() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();

		int[] siteIndexs = record.getAllSiteIndexs();
		int j= record.getSpectrumIndex(); 
		Spectrum spectrum = spectrumModel.getSpectrum(j);
		
		Set<Integer> allSrpPos = new HashSet<Integer>();
		for (int i = 0; i < siteIndexs.length; i++) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(siteIndexs[i]);
			allSrpPos.addAll(mapToSrp);
		}
		
		for (Integer integer : allSrpPos) {
			int i = integer;
			String fullSrp = aMap.getSrpFull(i);
		
			for (int s = 0; s < siteIndexs.length; s++) {
				int k = siteIndexs[s];
				if(allLogLikelihood[i][j][k]!=0){
					updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);
				}
			}
			updateEachSrpAt(i);

		}
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);

		return logLikelihood;


	}
	
	private double calculateSrpLikelihoodColumn() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int k = record.getColumnIndex();
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);

		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			
			for (int j = 0; j < spectrumCount; j++) {
				Spectrum spectrum = spectrumModel.getSpectrum(j);
				updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);

			}
			updateEachSrpAt(i);

		}
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;
	}


	private double calculateSrpLikelihoodSwapSubColumn() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int k = record.getColumnIndex();
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
		int[] allSpectrumIndexs = record.getAllSpectrumIndexs();
		
		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			
//			for (int j = 0; j < spectrumCount; j++) {
			for (int j : allSpectrumIndexs) {
				Spectrum spectrum = spectrumModel.getSpectrum(j);
				updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);

			}
			updateEachSrpAt(i);

		}
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;
	}

	
	private double calculateSrpLikelihoodRecombination() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int[] twoSpectrums = record.getRecombinationSpectrumIndex();
		int[] twoPositions = record.getRecombinationPositionIndex();

//		Spectrum[] spectrums = new Spectrum[] {
//				spectrumModel.getSpectrum(twoSpectrums[0]),
//				spectrumModel.getSpectrum(twoSpectrums[1]) };
	
		int j0 = twoSpectrums[0];
		int j1 = twoSpectrums[1];
		
		for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
			
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
//			System.out.println("Site: "+k +"\t"+ mapToSrp.size());
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}
		
		Set<Integer> allSrpPos = new HashSet<Integer>();
		allSrpPos.clear();
		for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(i);
			allSrpPos.addAll(mapToSrp);
		}


		for (int i = 0; i < srpSwitch.length; i++) {
			if(srpSwitch[i]){
	
//				String fullSrp = aMap.getSrpFull(i);
//				for (int r = 0; r < TWO; r++) {
				double LL0 = storedSpectrumLogLikelihood[i][j0];
				double LL1 = storedSpectrumLogLikelihood[i][j1];
				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
					if(storedAllLogLikelihood[i][j0][k]!=0){

						double L0 = storedAllLogLikelihood[i][j0][k];
						double L1 = storedAllLogLikelihood[i][j1][k];

						LL0 += (-L0 + L1);
						LL1 += (-L1 + L0);

						allLogLikelihood[i][j1][k] = L0;
						allLogLikelihood[i][j0][k] = L1;
					}

				}
				spectrumLogLikelihood[i][j0] = LL0;
				spectrumLogLikelihood[i][j1] = LL1;
				updateEachSrpAt(i);
			}

		}
	
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;

	}

	@Deprecated
	private double calculateSrpLikelihoodRecombination_full() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int[] twoSpectrums = record.getRecombinationSpectrumIndex();
		int[] twoPositions = record.getRecombinationPositionIndex();

		Spectrum[] spectrums = new Spectrum[] {
				spectrumModel.getSpectrum(twoSpectrums[0]),
				spectrumModel.getSpectrum(twoSpectrums[1]) };

		for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}
		for (int i = 0; i < srpSwitch.length; i++) {
			if(srpSwitch[i]){
				
				String fullSrp = aMap.getSrpFull(i);
//				for (int r = 0; r < TWO; r++) {
//					int j=twoSpectrums[r];	
					for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
						
						int j=twoSpectrums[0];	
						if(allLogLikelihood[i][j][k]!=0){
							
							updateLikelihoodAtIJK(i,j,k, spectrums[0], fullSrp);
							j=twoSpectrums[1];	
							updateLikelihoodAtIJK(i,j,k, spectrums[1], fullSrp);
						}
					}
//				}
				updateEachSrpAt(i);
			}

		}

		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;

	}

	private double[] calculateStatesLogLikelihood(Spectrum spectrum, int k) {
		SpectraParameter spectra = spectrum.getSpectra(k);
		double[] statesLogLikelihood = new double[AMBIGUOUS_STATE_COUNT];
		// Arrays.fill(stateLogLikelihood, LOG_ERROR_RATE);
		for (int state = 0; state < STATE_COUNT; state++) {
			double frequency = spectra.getFrequency(state);
			statesLogLikelihood[state] = Math.log(frequency * NOT_ERROR_RATE
					+ (1 - frequency) * ERROR_RATE);
		}
		statesLogLikelihood[GAP_STATE] = LOG_ERROR_RATE;
		return statesLogLikelihood;
	}

	private void updateEachSrpAt(int i) {
		double temp = eachSrpLogLikelihood[i];
		liS.reset();
		for (int j = 0; j < spectrumCount; j++) {
			liS.addScaleLogProb(spectrumLogLikelihood[i][j]);
		}
//		System.out.print("eachsrp" +"\t"+ eachSrpLikelihood[i]);
		eachSrpLogLikelihood[i] = liS.getLogLikelihood();
		if(eachSrpLogLikelihood[i] != temp){
			System.out.println("diff eachSrp: "+i +"\t"+ temp +"\t"+ eachSrpLogLikelihood[i]);
		}
//		System.out.println(eachSrpLikelihood[i]);
	}
	private void updateLikelihoodAtIJK(int i, int j, int k, Spectrum spectrum, String fullSrp) {
		
		char srpChar = fullSrp.charAt(k);
		int state = dataType.getState(srpChar);
		double frequency = 0;
		double logLikelihood;
		if(state<STATE_COUNT){
			 frequency = spectrum.getFrequency(k, state);
			double likelihood = frequency * NOT_ERROR_RATE
					+ (1 - frequency) * ERROR_RATE;
			logLikelihood = Math.log(likelihood);
		}
		else{
			logLikelihood = LOG_ERROR_RATE;
		}
//		System.out.println(i +"\t"+ j +"\t"+ k +"\t"+ state +"\t"+ frequency +"\t"+ logLikelihood);
//		System.out.print(spectrumLogLikelihood[i][j] +"\t" );
		if(allLogLikelihood[i][j][k] != logLikelihood){
			System.out.println("should be equal\t"+ allLogLikelihood[i][j][k] +"\t"+ logLikelihood);
		}
		
		spectrumLogLikelihood[i][j] -= allLogLikelihood[i][j][k]; 
		allLogLikelihood[i][j][k] = logLikelihood;
		spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
//		System.out.println(spectrumLogLikelihood[i][j] +"\t"+ logLikelihood);
	}


	private double updateLikelihoodAtIJK(int i, int j, int k, String fullSrp, double[] statesLogLikelihood, double totalLogLikelihood){
		char srpChar = fullSrp.charAt(k);
		int state = dataType.getState(srpChar);
//		int state = 0;
		
		double oneStateLogLikelihood= statesLogLikelihood[state];
//		double totalLogLikelihood = logLikelihood;		
		
		if(allLogLikelihood[i][j][k] != oneStateLogLikelihood){
			
//			System.out.println("Diff StateLog:\t"+state +"\t"+ oneStateLogLikelihood+"\t"+ allLogLikelihood[i][j][k]);
//			System.out.println(StatUtils.sum(eachSrpLogLikelihood) +"\t"+ this.logLikelihood);
			totalLogLikelihood -= eachSrpLogLikelihood[i];

			liS.reset();
			liS.addScaleLogProb(eachSrpLogLikelihood[i]);
//			liS.setsumScaledLikelihood(sumScaledSrpLogLikelihood[i]);
			liS.minusScaleLogProb( spectrumLogLikelihood[i][j]);
			
			spectrumLogLikelihood[i][j] -= allLogLikelihood[i][j][k]; 
			allLogLikelihood[i][j][k] = oneStateLogLikelihood;
			spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];

			liS.addScaleLogProb(spectrumLogLikelihood[i][j]);
//			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
			totalLogLikelihood += eachSrpLogLikelihood[i];

		}
		return totalLogLikelihood;
	}
	
	@Deprecated
	private void updateLikelihoodAtIJK2(int i, int j, int k, Spectrum spectrum, String fullSrp) {
		
		char srpChar = fullSrp.charAt(k);
		int state = dataType.getState(srpChar);
		double logLikelihood = 0;
		if(state<STATE_COUNT){
			double[] frequencies = spectrum.getFrequenciesAt(k);
			double likelihood = frequencies[state] * NOT_ERROR_RATE
					+ (1 - frequencies[state]) * ERROR_RATE;
			logLikelihood = Math.log(likelihood);
		}
		else{
			logLikelihood = LOG_ERROR_RATE;
		}

		spectrumLogLikelihood[i][j] -= allLogLikelihood[i][j][k]; 
		allLogLikelihood[i][j][k] = logLikelihood;
		spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
		
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == spectrumModel) {
            // treeModel has changed so recalculate the intervals
//            eventsKnown = false;
        }
        else{
        	System.err.println("Call handleModelChangedEvent in ShortReadSpectrumLikelihood" +"\t"+ model.getModelName());
        }
        likelihoodKnown = false;
//        System.err.println("Call handleModelChangedEvent in ShortReadSpectrumLikelihood" +"\t"+ model.getModelName());
//        makeDirty();
        
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent in SpectrumAlignmentModel");
	}

	@Override
	protected void storeState() {
//long time1 = System.currentTimeMillis();

//		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		storeEverything();
		SpectrumOperationRecord spectrumOperationRecord = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex; = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;// = aMap.getMapToSrp(siteIndex);

		int j;
		int k;
		switch (operation) {
		case NONE:
			if(DEBUG){
				System.out.println("StoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("StoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			storeEverything();
			break;

		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
			k= spectrumOperationRecord.getColumnIndex();
			mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
//				storedEachSrpLikelihood[i] = eachSrpLikelihood[i];
				for (int sp = 0; sp < spectrumCount; sp++) {
					storedAllLogLikelihood[i][sp][k] = allLogLikelihood[i][sp][k];
					storedSpectrumLogLikelihood[i][sp] = spectrumLogLikelihood[i][sp];
				}
			}

			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getAllSiteIndexs()[0];
			mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
				storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
				storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
				storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
			}
			
		case DELTA_MULTI:
		case SWAP_MULTI:
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			
			if(useSrpArray){

				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
						for (int s = 0; s < siteIndexs.length; s++) {
							k = siteIndexs[s];
							if (allLogLikelihood[i][j][k] != 0) {
									storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
									storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
							}
						}
						srpSwitch[i] = false;
					}
				}
			}
			else{
				for (int i : allSrpPos) {
					storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
					for (int s = 0; s < siteIndexs.length; s++) {
						k = siteIndexs[s];
						if (allLogLikelihood[i][j][k] != 0) {
								storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
								storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
						}
					}
				}
			}
			
			
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();

			for (int i = 0; i < srpSwitch.length; i++) {
				if (srpSwitch[i]) {
//					storedEachSrpLikelihood[i] = eachSrpLikelihood[i];
					for (int r = 0; r < TWO; r++) {
						j = twoSpectrums[r];
						for ( k = twoPositions[0]; k < twoPositions[1]; k++) {
							if (allLogLikelihood[i][j][k] != 0) {
								storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
								storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
							}
						}
					}
					srpSwitch[i] = false;
				}
			}
			
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsSpectrumLikelihood.class.getSimpleName() );
//			storeEverything();
//			break;
			
		}


//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
		
		
//		System.out.println(spectrumIndex +"\t"+ siteIndex);
//		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(siteIndex);
//		for (int i : mapToSrp) {
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
//			storedSpectrumLogLikelihood[i][spectrumIndex] = spectrumLogLikelihood[i][spectrumIndex];
//		}
		
//		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
//			for (int j = 0; j < spectrumCount; j++) {
//				for (int k = 0; k < spectrumLength; k++) {
//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
//				}
//				System.arraycopy(allLogLikelihood[i][j], 0, storedAllLogLikelihood[i][j], 0, spectrumLength);
//			}
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
//			
//		}
		
		
//		System.err.println("SR likelihood store: " + logLikelihood +"\t"+ storedLogLikelihood);
//		System.err.println(allLogLikelihood);
	}
//	private void storeMultiDelta(SpectrumOperationRecord spectrumOperationRecord) {
//		
//		
//	}

	@Override
	protected void restoreState() {
//		long time1 = System.currentTimeMillis();
		
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
//		restoreEverything();
//		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		
		SpectrumOperationRecord spectrumOperationRecord = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;// = aMap.getMapToSrp(siteIndex);
//		int[] siteIndexs;
		int j;
		int k;
		switch (operation) {
		case NONE:
			if(DEBUG){
				System.out.println("RestoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("RestoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			restoreEverything();
			break;
		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
			k = spectrumOperationRecord.getColumnIndex();
			mapToSrp = aMap.getMapToSrp(k);
			
			for (j = 0; j < spectrumCount; j++) {
				for (int i : mapToSrp) {
					allLogLikelihood[i][j][k] = storedAllLogLikelihood[i][j][k];
					spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
				}
			}

			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getAllSiteIndexs()[0];
			mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
				sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
				spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
				allLogLikelihood[i][j][k] = storedAllLogLikelihood[i][j][k];
			}
			
			
		case DELTA_MULTI:
		case SWAP_MULTI:
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			
			if(useSrpArray){
				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
						for (int s = 0; s < siteIndexs.length; s++) {
							k = siteIndexs[s];
							if (allLogLikelihood[i][j][k] != 0) {
								allLogLikelihood[i][j][k] = storedAllLogLikelihood[i][j][k];
								spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
							}
						}
					}
				}
			}
			else{
				for (int i : allSrpPos) {
					storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
					for (int s = 0; s < siteIndexs.length; s++) {
						k = siteIndexs[s];
						if (allLogLikelihood[i][j][k] != 0) {
							allLogLikelihood[i][j][k] = storedAllLogLikelihood[i][j][k];
							spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
						}
					}
					
				}
				
			}
			
			
//			
//			for (int s = 0; s < siteIndexs.length; s++) {
//				mapToSrp = aMap.getMapToSrp(siteIndexs[s]);
//				for (int i : mapToSrp) {
//					allLogLikelihood[i][spectrumIndex][siteIndexs[s]] = storedAllLogLikelihood[i][spectrumIndex][siteIndexs[s]];
//					spectrumLogLikelihood[i][spectrumIndex] = storedSpectrumLogLikelihood[i][spectrumIndex];
//				}
//			}
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();

			for (int i = 0; i < srpSwitch.length; i++) {
				if (srpSwitch[i]) {
//					storedEachSrpLikelihood[i] = eachSrpLikelihood[i];
					for (int r = 0; r < TWO; r++) {
						j = twoSpectrums[r];
						for ( k = twoPositions[0]; k < twoPositions[1]; k++) {
							if (allLogLikelihood[i][j][k] != 0) {
								allLogLikelihood[i][j][k] = storedAllLogLikelihood[i][j][k];
								spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
							}
						}
					}
				}
			}

			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsSpectrumLikelihood.class.getSimpleName() );

		}
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
//
//		
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex();
//		int siteIndex = record.getSiteIndex();
//
//		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(siteIndex);
//		for (int i : mapToSrp) {
//			allLogLikelihood[i][spectrumIndex][siteIndex] = storedAllLogLikelihood[i][spectrumIndex][siteIndex];
//			spectrumLogLikelihood[i][spectrumIndex] = storedSpectrumLogLikelihood[i][spectrumIndex];
//		}
//		
		
//		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			allLogLikelihood[i][spectrumIndex][siteIndex] = storedAllLogLikelihood[i][spectrumIndex][siteIndex];
//			for (int j = 0; j < spectrumCount; j++) {
//				
//				System.arraycopy(storedAllLogLikelihood[i][j], 0, allLogLikelihood[i][j], 0, spectrumLength);
//			}
//			
//		}
//		System.err.println(
//				spectrumIndex +"\t"+ siteIndex +"\t"
//						+Arrays.toString(spectrumModel.getSpectrum(spectrumIndex)
//						.getFrequencies(siteIndex)));
//
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		
	}

	@Override
	protected void acceptState() {
		//Do nothing
	}
	private void restoreEverything(){
		
		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
			System.arraycopy(storedSpectrumLogLikelihood[i],0, spectrumLogLikelihood[i], 0, spectrumCount);
			for (int j = 0; j < spectrumCount; j++) {
//				for (int k = 0; k < spectrumLength; k++) {
//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
//				}
				System.arraycopy(storedAllLogLikelihood[i][j], 0, allLogLikelihood[i][j], 0, spectrumLength);
			}
			
		}
	}
	private void storeEverything() {

		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
//		storedLogLikelihood = logLikelihood;
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex();
//		int siteIndex = record.getSiteIndex();
//		System.out.println(spectrumIndex +"\t"+ siteIndex);
		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
			System.arraycopy(spectrumLogLikelihood[i],0, storedSpectrumLogLikelihood[i], 0, spectrumCount);
			for (int j = 0; j < spectrumCount; j++) {
//				for (int k = 0; k < spectrumLength; k++) {
//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
//				}
				System.arraycopy(allLogLikelihood[i][j], 0, storedAllLogLikelihood[i][j], 0, spectrumLength);
			}
			
		}
	}

	public SpectrumOperation getOperation(){
		return spectrumModel.getSpectrumOperation();
	}
	
	@Override
	public Model getModel() {
		return this;
		
	}


	@Override
	public void makeDirty() {
//		System.err.println("make dirty");
		spectrumModel.resetSpectrumOperation();
        likelihoodKnown = false;
		
	}

	public double[] unittestMethodGetEachLikelihood() {
		double[] copyOfValues = new double[eachSrpLogLikelihood.length];
        System.arraycopy(eachSrpLogLikelihood, 0, copyOfValues, 0, copyOfValues.length);
		return copyOfValues;
	}


	@Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }


}