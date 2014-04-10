package srp.likelihood.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;

import srp.evolution.OperationRecord;
import srp.evolution.OperationType;
import srp.evolution.shortreads.ShortRead;
//import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.evolution.spectrum.SpectraParameter;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import srp.likelihood.AbstractShortReadsLikelihood;
import srp.likelihood.LikelihoodScaler;
//import srp.likelihood.spectrum.AbstractShortReadsSpectrumLikelihood;








import srp.likelihood.stateLikelihood.StateLikelihood;

import com.carrotsearch.hppc.BitSet;
//import java.util.BitSet;



public class ShortReadsHaplotypeLikelihood  extends AbstractShortReadsLikelihood {

	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadHaplotypeLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;

	private final double MIN_LOG_LIKELIHOOD;

	

	private LikelihoodScaler liS;



	
	protected HaplotypeModel alignmentModel;
	
	@Deprecated protected double[] allStateLogLikelihood;
	@Deprecated protected double[] allStoredStateLogLikelihood;

//	protected HaplotypeModel haplotypeModel;
	
	private HashMap<Integer, double[]> scaledLogBinomialDesnity;
	private int[][] allDists;
	private int[][] storedAllDists;
	
	@Deprecated protected StateLikelihood stateLikelihood;

	

	public ShortReadsHaplotypeLikelihood(HaplotypeModel haplotypeModel, ShortReadMapping srpMap){
		super(SHORT_READ_LIKELIHOOD, srpMap);
		this.alignmentModel = haplotypeModel;

		operationRecord = alignmentModel.getOperationRecord();
//		multiType = MultiType.Array;
		multiType = MultiType.BitSet;
		
//		type = MultiType.Hash;
//		type = MultiType.All;
//		distTypeCode = "flat";//"betaMean"  "betaMode" "gTest"
//		setDistType(distType);
		MIN_LOG_LIKELIHOOD = 0;//stateLikelihood.caluclateStateLogLikelihood(SpectraParameter.MIN_FREQ);
		
		

		likelihoodKnown = false;
		
		addModel(this.alignmentModel);
		
		preprocessLikelihoodAlignmentMap();

//		calculateSrpLikelihoodFull();//TODO FIX this? shouldn't needed
		getLogLikelihood();
				
		storeEverything();
		
		for (int j = 0; j < sequenceCount; j++) {
			Haplotype haplotype = haplotypeModel.getHaplotype(j);
			haplotype.storeState();
		}
		
	}
	
	

	private void preprocessLikelihoodAlignmentMap() {
//		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
//		srpCount = srpMap.getSrpCount();
		sequenceCount = alignmentModel.getHaplotypeCount();
//		sequenceLength = alignmentModel.getHaplotypeLength();
		operationRecord = alignmentModel.getOperationRecord();
//		this.srpSwitch = new boolean[srpCount];
//		this.allSrpPos = new HashSet<Integer>();
//		this.srpIndex = new int[srpCount];
		
//		logLikelihood = Double.NEGATIVE_INFINITY;
//		storedLogLikelihood = Double.NEGATIVE_INFINITY;

//		spectrumLogLikelihood = new double[srpCount*sequenceCount];
//		storedSpectrumLogLikelihood = new double[srpCount*sequenceCount];
//
//		scaledSpectrumLogLikelihood = new double[srpCount*sequenceCount];
//		storedScaledSpectrumLogLikelihood = new double[srpCount*sequenceCount];
//		
//		sumScaledSrpLogLikelihood = new double[srpCount];
//		storedSumSrpLogLikelihood = new double[srpCount];
//
//		
//		eachSrpLogLikelihood = new double[srpCount];
//		storedEachSrpLogLikelihood = new double[srpCount];

//		allStateLogLikelihood = new double[sequenceLength*STATE_COUNT];
//		allStoredStateLogLikelihood = new double[sequenceLength*STATE_COUNT];
		
//		mapToSrpArray = srpMap.getMapToSrpArray();
//		
//		
//		String[] srpArray = srpMap.getSrpArray();
//		
//		allSrpState2D = new int[srpArray.length][sequenceLength];
//		allSrpChar2D = new char[srpArray.length][sequenceLength];
//		
//		for (int i = 0; i < srpArray.length; i++) {
//			String srp = srpArray[i];
//			for (int j = 0; j < sequenceLength; j++) {
//				allSrpState2D[i][j] = getStateAtK(srp, j);
//				allSrpChar2D[i][j] = srp.charAt(j);
//			}
//		}
//		
//		
//		bitSet = new BitSet(srpCount);
		
		
		scaledLogBinomialDesnity = new HashMap<Integer, double[]>();
		
		allDists = new int[srpCount][sequenceCount];
		storedAllDists = new int[srpCount][sequenceCount];


		
		int maxDist=0;
		for (int s = 0; s < srpCount; s++) {
			String srp = srpMap.getSrpFragment(s);

			int srLength = srp.length();
			int srLength1 = srLength+1;
			maxDist = Math.max(maxDist, srLength1);
			double[] logBinomD = new double[srLength1];
			double[] scaledBinomD = new double[srLength1];
			for (int i = 0; i < logBinomD.length; i++) {
	
				logBinomD[i] = i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
	//			logBinomD[i] = i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
				scaledBinomD[i] = liS.scale(logBinomD[i]); 
			}
	//		System.out.println(Arrays.toString(logBinomD));
//			logBinomialDesnity.put(srLength, logBinomD);
			scaledLogBinomialDesnity.put(srLength, scaledBinomD);
		}
		
//		for (int i = 0; i < srpCount; i++) {
//			String srp = srpMap.getSrpFragment(i);//srpArray[i]
//			int start = srpMap.getSrpStart(i);
//			int end = srpMap.getSrpEnd(i);
//			for (int j = 0; j < sequenceCount; j++) {
//				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getSequenceString());
//				allDists[i][j]=dist;
//			}
//		}
	}


    
	protected double calculateSrpLikelihoodFull() {


//		System.out.println("calculateSrpLikelihoodFull");
		for (int i = 0; i < srpCount; i++) {
			
			String srp = srpMap.getSrpFragment(i);//srpArray[i]
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));

			liS.reset();
			for (int j = 0; j < sequenceCount; j++) {

				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getSequenceString());
				allDists[i][j]=dist;
				liS.add(logPD[dist]);
//				liS.scaleLogProb(logPD[dist]);
			}	
			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
//			sumScaledSrpLogLikelihood[i] = liS.scale(eachSrpLogLikelihood[i], LOG_C);

		}
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;
//		System.out.println("calculateSrpLikelihoodFull");
//
//		double[][] allStateLogLikelihoodFull2D = new double[sequenceCount][sequenceLength*STATE_COUNT];
//		
////		for (int j = 0; j < sequenceCount; j++) {
////			Haplotype haplotype = haplotypeModel.getHaplotype(j);
////			for (int k = 0; k < sequenceLength; k++) {
////				SpectraParameter spectra = spectrum.getSpectra(k);
////				int kOffset = k*STATE_COUNT;
////				stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
////			}
////			System.arraycopy(allStateLogLikelihood, 0  , allStateLogLikelihoodFull2D[j], 0, allStateLogLikelihood.length );
////		}
//
//		double stateLogLikelihood;
//		for (int i = 0; i < srpCount; i++) {
//
////			String fullSrp = srpMap.getSrpFull(i);
//			int start = srpMap.getSrpStart(i);
//			int end = srpMap.getSrpEnd(i);
//			
//			liS.reset();
//			for (int j = 0; j < sequenceCount; j++) {
//				int offset = i*sequenceCount+j;
////				Haplotype haplotype = spectrumModel.getHaplotype(j);
//
////				spectrumLogLikelihood[offset] = 0;
//				for (int k = start; k < end; k++) {
////					int state = getStateAtK(fullSrp, k);
//					int state = allSrpState2D[i][k];
//					
//					if(state<STATE_COUNT){
//						int kOffset = k * STATE_COUNT+state;
//						 stateLogLikelihood = allStateLogLikelihoodFull2D[j][kOffset];
//					}
//					else{
//						stateLogLikelihood = MIN_LOG_LIKELIHOOD;
//					}
////					spectrumLogLikelihood[offset] += stateLogLikelihood;
//					
//				}
//				
////				scaledSpectrumLogLikelihood[offset] = liS.scale(spectrumLogLikelihood[offset]);
//				
////				liS.add(scaledSpectrumLogLikelihood[offset]);
//
//			}	
////			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
//			
//			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
////			System.out.println(i +"\t"+ eachSrpLogLikelihood[i]);
//		}
//		
//		double logLikelihood = liS.sumLogLikelihood(sumScaledSrpLogLikelihood);
//		double totalLogLikelihood = StatUtils.sum(eachSrpLogLikelihood);
////		totalLogLikelihood = calculateSrpLikelihoodFullMaster();
//		return totalLogLikelihood;
	}



	
	protected double calculateSrpLikelihoodSingle() {


//		OperationRecord record = alignmentModel.getOperationRecord();
		int j = operationRecord.getSpectrumIndex(); 
		int k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
		double currentLogLikelihood = getStoredLogLikelihood();
		
//		System.out.println("StartSingle" +"\t"+ j +"\t"+ k +"\t"+ currentLogLikelihood);
//		SpectraParameter spectra = haplotypeModel.getHaplotype(j).getSpectra(k);
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);

//		stateLikelihood.calculateStatesLogLikelihood2D(spectra, allStateLogLikelihood2D[0]);
//		stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra, allStoredStateLogLikelihood2D[0]);

//		stateLikelihood.calculateStatesLogLikelihood(spectra, 0, allStateLogLikelihood);
//		stateLikelihood.calculateStoredStatesLogLikelihood(spectra, 0, allStoredStateLogLikelihood);

////		for (int i : mapToSrp) {
		Haplotype haplotype = alignmentModel.getHaplotype(j);
		char oldChar = haplotype.getStoredChar(k);
		char newChar = haplotype.getChar(k);
//		char 
		

	
		if(newChar!= oldChar){ // if(newChar!= oldChar && isHapEqualNew)
	//		getChar(k);
			for (int i : mapToSrpArray[k]){
				int srpLength = srpMap.getSrpLength(i);
	//			int state = getStateAtK(fullSrp, k);
	//			int state = allSrpState2D[i][k];
				
				char srpChar = allSrpChar2D[i][k];
				int deltaDist = 0;
				
				if (srpChar==newChar){
					deltaDist = -1;
				}
				else if(srpChar==oldChar){
					deltaDist = 1;
				}
				
				if (deltaDist != 0) {
//					currentLogLikelihood = updateLikelihoodAtIJK(deltaDist
//	//						allStateLogLikelihood, allStoredStateLogLikelihood,
//							currentLogLikelihood);

					int srpIndex = i; int hapIndex = j; 
//					int swapPos = 0; 
//					int newChar = 0; 
//					int oldChar = 0;
////					int srpIndex, int hapIndex, int swapPos, int newChar, int oldChar
//					ShortRead srp = srpMap.getShortRead(srpIndex);
//					int srpChar = srp.getFullSrpCharAt(swapPos);
//					
//					int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);

					double[] logPD = scaledLogBinomialDesnity.get(srpLength);
					int oldDist = storedAllDists[srpIndex][hapIndex];
					int newDist = storedAllDists[srpIndex][hapIndex] + deltaDist;
					allDists[srpIndex][hapIndex] = newDist;
					
					currentLogLikelihood -= eachSrpLogLikelihood[i];
					
					sumScaledSrpLogLikelihood[i] -= logPD[oldDist];//scaledSpectrumLogLikelihood[offset];
			
//						spectrumLogLikelihood[offset] -= storedStateLn; 
//						spectrumLogLikelihood[offset] += stateLn;
//						scaledSpectrumLogLikelihood[offset] = LikelihoodScaler.scale(spectrumLogLikelihood[offset], LOG_C);
					
					sumScaledSrpLogLikelihood[i] += logPD[newDist];//scaledSpectrumLogLikelihood[offset];

					eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
					
					currentLogLikelihood += eachSrpLogLikelihood[i];

				}
			}
		}
		double trueL = calculateSrpLikelihoodFullMaster();
		if((trueL - currentLogLikelihood)> 1e-6) {
			System.out.println("DEBUG "+trueL +"\t"+ currentLogLikelihood);
			
			for (int i: mapToSrpArray[k]) {
				
				String srp = srpMap.getSrpFragment(i);//srpArray[i]
				int start = srpMap.getSrpStart(i);
				int end = srpMap.getSrpEnd(i);
				
				double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));

				liS.reset();
				double oldPD = 0;
				double newPD = -1;
				for (int hj = 0; hj < sequenceCount; hj++) {

					int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(hj).getSequenceString());
//					System.out.println(hj +" "+ j +":\t"+ dist +"\t"+ storedAllDists[i][hj] +"\t"+ (dist == storedAllDists[i][hj]));
					liS.add(logPD[dist]);
					
					if(hj == j){
						int srpLength = srpMap.getSrpLength(i);
						char srpChar = allSrpChar2D[i][k];
						int deltaDist = 0;
						if (srpChar==newChar){
							deltaDist = -1;
						}
						else if(srpChar==oldChar){
							deltaDist = 1;
						}
						if (deltaDist != 0) {
		
							int srpIndex = i; int hapIndex = j; 
			
							logPD = scaledLogBinomialDesnity.get(srpLength);
							int oldDist = storedAllDists[srpIndex][hapIndex];
							int newDist = storedAllDists[srpIndex][hapIndex] + deltaDist;
//							System.out.println();
							oldPD = logPD[oldDist];
							newPD = logPD[newDist];
							System.out.println(hj + " " + j + ":\t" + dist
									+ "\t" + storedAllDists[i][hj] + "\t" +
									logPD[oldDist] +"\t"+ logPD[newDist]
									+ (dist == storedAllDists[i][hj]));
						}
					}
//					liS.scaleLogProb(logPD[dist]);
				}	
				
				double e1 = liS.getSumScaledLikelihood(); 
				double e2 = liS.getLogLikelihood();
				
				System.out.println((e1 == sumScaledSrpLogLikelihood[i]));
				double d1 = e1 - newPD;
				double d2 = storedSumScaledSrpLogLikelihood[i] - oldPD;
				double d3 = sumScaledSrpLogLikelihood[i] - newPD;
				System.out.println((d1==d2) +"\t"+ d1 +"\t"+ d2);
				System.out.println((d2==d3) +"\t"+ d2 +"\t"+ d3);
				if(e2 != 	eachSrpLogLikelihood[i]){
					System.out.println(i +"\t"+ eachSrpLogLikelihood[i] + "\t" + e2 +"\t"+ 
							(eachSrpLogLikelihood[i] -e2) +"\t"+ 
							e1 +"\t"+  
							sumScaledSrpLogLikelihood[i]
							
							);
					

					
				}
				else{
					System.out.println(i);
				}
					
			}
		}
		return currentLogLikelihood;
	}
	
	
	protected double calculateSrpLikelihoodMulti() {
		
		OperationRecord record = alignmentModel.getOperationRecord();

		int[] siteIndexs = record.getAllSiteIndexs();
		int j= record.getSpectrumIndex(); 
		Haplotype haplotype = alignmentModel.getHaplotype(j);
		
//		stateLikelihood.calculateStatesLogLikelihood(spectra, allStateLogLikelihood2D[k]);
//		stateLikelihood.calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood2D[k]);
		
//		for (int k : siteIndexs) {
//			SpectraParameter spectra = spectrum.getSpectra(k);
//			int kOffset = k*STATE_COUNT;
//			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
//		}
		
		int multihere;
		
		double currentLogLikelihood = getStoredLogLikelihood();
		if(multiType == MultiType.BitSet){
			recalculateBitSet(siteIndexs);
			
			int count = 0;
			for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i+1)) {
				srpIndex[count++] = i;

				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j,
						siteIndexs, allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);

			}
			srpIndexCount = count;

			
		}
		else if(multiType==MultiType.Array){
			recalculateArray(siteIndexs);

			for (int i = 0; i < srpSwitch.length; i++) {
				if(srpSwitch[i]){

					currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j, siteIndexs,
							allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);

				}
			}
	
		}
		else if(multiType==MultiType.Hash){
			recalculateHashSet(siteIndexs);

			for (int i : allSrpPos) {
				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j, siteIndexs,
						allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);
			}
		}
		
		
		

		return currentLogLikelihood;

	}


	protected double calculateSrpLikelihoodColumn() {

		OperationRecord record = alignmentModel.getOperationRecord();
		int k = record.getSingleIndex();
		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
		int[] allSpectrumIndexs = record.getAllSpectrumIndexs();
		double currentLogLikelihood = getStoredLogLikelihood();
		
//		for (int j : allSpectrumIndexs) {
//			SpectraParameter spectra = haplotypeModel.getHaplotype(j).getSpectra(k);
//			
//			int kOffset = j*STATE_COUNT;
//			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
//		}
		for (int i : mapToSrp) {
			int state = allSrpState2D[i][k];
			for (int j : allSpectrumIndexs) {
//				currentLogLikelihood = updateLikelihoodAtIJK(i, j, state, allStateLogLikelihood2D[j],
//						allStoredStateLogLikelihood2D[j], currentLogLikelihood);
				if (state < STATE_COUNT) {
					currentLogLikelihood = updateLikelihoodAtIJK(i, j, j*STATE_COUNT+state,
							currentLogLikelihood);
				}
			}
		}
		return currentLogLikelihood;
	}
	
	protected double calculateSrpLikelihoodSubColumn() {

		OperationRecord record = alignmentModel.getOperationRecord();
		int k = record.getSingleIndex();
		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
		int[] allSpectrumIndexs = record.getAllSpectrumIndexs();
		
		double currentLogLikelihood = getStoredLogLikelihood();

//		for (int j = 0; j < spectrumCount; j++) {
//		for (int j : allSpectrumIndexs) {
//			SpectraParameter spectra = haplotypeModel.getHaplotype(j).getSpectra(k);
//			
//			int kOffset = j*sequenceCount;
//			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
//	
//		}
		
		for (int i : mapToSrp) {
//			String fullSrp = srpMap.getSrpFull(i);
//			int state = getStateAtK(fullSrp, k);
			int state = allSrpState2D[i][k];
//			for (int j = 0; j < spectrumCount; j++) {
			for (int j : allSpectrumIndexs) {
				if (state < STATE_COUNT) {
					currentLogLikelihood = updateLikelihoodAtIJK(i, j, j*STATE_COUNT+state,
							currentLogLikelihood);
				}
			}
			

		}
		return currentLogLikelihood;
	}

	

	protected double calculateSrpLikelihoodRecombination() {
		
		OperationRecord record = alignmentModel.getOperationRecord();
		int[] twoSpectrums = record.getRecombinationSpectrumIndex();
		int[] twoPositions = record.getRecombinationPositionIndex();

//		Spectrum[] spectrums = new Spectrum[] {
//				spectrumModel.getHaplotype(twoSpectrums[0]),
//				spectrumModel.getHaplotype(twoSpectrums[1]) };
	
		int j0 = twoSpectrums[0];
		int j1 = twoSpectrums[1];
		int length = twoPositions[1] - twoPositions[0];
		
		Haplotype haplotype0 = alignmentModel.getHaplotype(j0); 
		Haplotype haplotype1 = alignmentModel.getHaplotype(j1);
		
		int[] siteIndexs = new int[length];
		for (int k = twoPositions[0], s=0; k < twoPositions[1]; k++, s++) {
			siteIndexs[s] = k;
			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}

//		for (int k : siteIndexs) {
//			
//			SpectraParameter spectra0 = spectrum0.getSpectra(k);
//			SpectraParameter spectra1 = spectrum1.getSpectra(k);
//			int kOffset = k*STATE_COUNT;
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra0, kOffset, allStateLogLikelihood);
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra1, kOffset, allStoredStateLogLikelihood);
//			
////			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra0, allStateLogLikelihood2D[s]);
////			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra1, allStoredStateLogLikelihood2D[s]);
//			spectra0.setStateLikelihood(allStoredStateLogLikelihood, kOffset);
//			spectra1.setStateLikelihood(allStateLogLikelihood, kOffset);
//		
//		}

		int multihere;
		double currentLogLikelihood = getStoredLogLikelihood();
		if(multiType == MultiType.BitSet){

			bitSet.clear();
//			BitSet bitSet = new BitSet(srpCount);
			for (int s : siteIndexs) {
				BitSet tempSet = srpMap.getBitSet(s);
				bitSet.or(tempSet);
			}
			
			int count = 0;
			for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i+1)) {
				srpIndex[count] = i;
				count++;

				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j0,
						siteIndexs, allStoredStateLogLikelihood, allStateLogLikelihood, currentLogLikelihood);
				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j1,
						siteIndexs, allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);


			}
			srpIndexCount = count;

			
		}
		else if(multiType==MultiType.Array){

			for (int s : siteIndexs) {
				for (int i : mapToSrpArray[s]){
					srpSwitch[i] = true;
				}
			}

			for (int i = 0; i < srpSwitch.length; i++) {
				if(srpSwitch[i]){

					currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j0,
							siteIndexs, allStoredStateLogLikelihood, allStateLogLikelihood, currentLogLikelihood);
					currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j1,
							siteIndexs, allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);

				}
			}
	
		}
		

		
	
//		totalLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return currentLogLikelihood;

	}




	private double updateLikelihoodAtIJK(int i, int j, int state, 
//			double[] statesLogLikelihood, double[] storedStatesLogLikelihood, 
			double currentLogLikelihood) {
	
		int srpIndex = i; int hapIndex = j; 
		int swapPos = 0; 
		int newChar = 0; 
		int oldChar = 0;
//		int srpIndex, int hapIndex, int swapPos, int newChar, int oldChar
		ShortRead srp = srpMap.getShortRead(srpIndex);
		int srpChar = srp.getFullSrpCharAt(swapPos);
		
		int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);

		if (deltaDist!= 0){
			double[] logPD = scaledLogBinomialDesnity.get(srp.getLength());
			int oldDist = storedAllDists[srpIndex][hapIndex];
			int newDist = storedAllDists[srpIndex][hapIndex] + deltaDist;
			allDists[srpIndex][hapIndex] = newDist;
	
			liS.reset();		
			for (int s = 0; s < sequenceCount ; s++) {
				liS.add(logPD[allDists[srpIndex][s]]);
			}
			eachSrpLogLikelihood[srpIndex] = liS.getLogLikelihood();
			
			
			currentLogLikelihood -= eachSrpLogLikelihood[i];
			sumScaledSrpLogLikelihood[i] -= logPD[oldDist];//scaledSpectrumLogLikelihood[offset];
	
//			spectrumLogLikelihood[offset] -= storedStateLn; 
//			spectrumLogLikelihood[offset] += stateLn;
//			scaledSpectrumLogLikelihood[offset] = LikelihoodScaler.scale(spectrumLogLikelihood[offset], LOG_C);
			
			sumScaledSrpLogLikelihood[i] += logPD[newDist];//scaledSpectrumLogLikelihood[offset];

			eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
			currentLogLikelihood += eachSrpLogLikelihood[i];
		}

	//////////////
		
//			double stateLn= allStateLogLikelihood[state];
//			double storedStateLn = allStoredStateLogLikelihood[state];
//			int offset = i*sequenceCount+j;
//			if(storedStateLn != stateLn){
//				currentLogLikelihood -= eachSrpLogLikelihood[i];
//				sumScaledSrpLogLikelihood[i] -= scaledSpectrumLogLikelihood[offset];
//		
//				spectrumLogLikelihood[offset] -= storedStateLn; 
//				spectrumLogLikelihood[offset] += stateLn;
//
//				scaledSpectrumLogLikelihood[offset] = LikelihoodScaler.scale(spectrumLogLikelihood[offset], LOG_C);
//				
//				sumScaledSrpLogLikelihood[i] += scaledSpectrumLogLikelihood[offset];
//
//				eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
//				currentLogLikelihood += eachSrpLogLikelihood[i];
//			}

		return currentLogLikelihood;
	}

	private static int calculateDeltaDist(int srpChar, int newChar, int oldChar){//, boolean isHapEqualNew){
		
		int deltaDist = 0;
	
		if(newChar!= oldChar){ // if(newChar!= oldChar && isHapEqualNew)
			if (srpChar==newChar){
				deltaDist = -1;
			}
			else if(srpChar==oldChar){
				deltaDist = 1;
			}
		}
		return deltaDist;
		
	}


	private double updateLikelihoodAtIJ_Local1DArray(int i, int j, int[] siteIndexs, 
			double[] allStateLogLikelihood, double[] allStoredStateLogLikelihood, double currentLogLikelihood) {

		int offsetIJ = i*sequenceCount+j;
		currentLogLikelihood -= eachSrpLogLikelihood[i];
//		sumScaledSrpLogLikelihood[i] -= scaledSpectrumLogLikelihood[offsetIJ];

//		double localSpectrum = spectrumLogLikelihood[offsetIJ];
		int[] thisSrpStata = allSrpState2D[i];
		for (int k : siteIndexs) {
			int state = thisSrpStata[k];
			if (state < STATE_COUNT) {
				int offset = k*STATE_COUNT+state;
				
//				localSpectrum -= allStoredStateLogLikelihood[offset];
//				localSpectrum += allStateLogLikelihood[offset];
				

			}
		}
		
//		spectrumLogLikelihood[offsetIJ] = localSpectrum;
//		scaledSpectrumLogLikelihood[offsetIJ] = LikelihoodScaler.scale(localSpectrum, LOG_C);
//		sumScaledSrpLogLikelihood[i] += scaledSpectrumLogLikelihood[offsetIJ];

//		eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
		currentLogLikelihood += eachSrpLogLikelihood[i];

		return currentLogLikelihood;
	}
	

	@Override
	protected void storeState() {

//		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		storeEverything();
		OperationRecord spectrumOperationRecord = alignmentModel.getOperationRecord();
		OperationType operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex; = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;

		int j;
		int k;
		if(DEBUG){
			System.out.println("StoreState in ShortReadsSpectrumLikelihood:\t"+operation);
		}
		switch (operation) {
		case NONE:
			break;
		case FULL:
			storeEverything();
			break;

		case COLUMN:
		
		case SWAP_SUBCOLUMN:
			k= spectrumOperationRecord.getSingleIndex();
			mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				storeI(i);
				for (j = 0; j < sequenceCount; j++) {
					storeIJ(i, j);
				}
			}

			break;
		case SINGLE:
		
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getSingleIndex();//AllSiteIndexs()[0];
//			mapToSrp = srpMap.getMapToSrp(k);
//			spectrumModel.getSpectrum(j).getSpectra(k).storeState();
//			for (int i : mapToSrp) {
			
			for (int i : mapToSrpArray[k]){
				storeI(i);
				storeIJ(i, j);
				
//				storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
//				storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];
//				storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
//				storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];

			}
			break;
			
		case MULTI:
		
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			int storeMulti;
			if(multiType==MultiType.BitSet){
				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					storeI(i);
					storeIJ(i, j);
				}
			}
			
			else if(multiType==MultiType.Array){
				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						storeI(i);
						storeIJ(i, j);
						srpSwitch[i] = false;
					}
				}
			}
			else if(multiType==MultiType.Hash){
				for (int i : allSrpPos) {
					storeI(i);
					storeIJ(i, j);
				}
			}
			else if(multiType==MultiType.All){
				for (int kk : siteIndexs) {
//				for (int s = 0; s < siteIndexs.length; s++) {
//					k = siteIndexs[s];
					mapToSrp = srpMap.getMapToSrp(kk);
					for (int i : mapToSrp) {
						storeI(i);
						storeIJ(i, j);
					}
				}
			}
						
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
//			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();

//			for (int i = 0; i < srpSwitch.length; i++) {
//				if (srpSwitch[i]) {
//					storeI(i);
//					storeIJ(i, twoSpectrums[0]);
//					storeIJ(i, twoSpectrums[1]);
//					srpSwitch[i] = false;
//				}
//			}
			
			if(multiType==MultiType.BitSet){
				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					storeI(i);
					storeIJ(i, twoSpectrums[0]);
					storeIJ(i, twoSpectrums[1]);

				}
			}
			
			else if(multiType==MultiType.Array){
				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						storeI(i);
						storeIJ(i, twoSpectrums[0]);
						storeIJ(i, twoSpectrums[1]);
						srpSwitch[i] = false;
					}
				}
			}
			
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsHaplotypeLikelihood.class.getSimpleName() );
//			storeEverything();
//			break;
			
		}


//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
		
		
//		System.out.println(spectrumIndex +"\t"+ siteIndex);
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(siteIndex);
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


	private void storeI(int i) {
		storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
		storedSumScaledSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
		
	}


	private void storeIJ(int i, int j) {
		storedAllDists[i][j] = allDists[i][j];
		
	}


	//	private void storeMultiDelta(SpectrumOperationRecord spectrumOperationRecord) {
	//		
	//		
	//	}
	
	
	private void storeEverything() {

		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(sumScaledSrpLogLikelihood, 0, storedSumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
		for (int i = 0; i < allDists.length; i++) {
			System.arraycopy(allDists[i], 0, storedAllDists[i], 0, sequenceCount);
		}
		
	}


	@Override
	protected void restoreState() {
//		long time1 = System.currentTimeMillis();
		
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
//		restoreEverything();
//		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		
		OperationRecord spectrumOperationRecord = alignmentModel.getOperationRecord();
		OperationType operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;
//		int[] siteIndexs;
		int j;
		int k;
		if(DEBUG){
			System.out.println("RestoreState in ShortReadsSpectrumLikelihood:\t"+operation);
		}
		switch (operation) {
		case NONE:
			
			break;
		case FULL:
			
			restoreEverything();
			break;
		case COLUMN:
		
			k = spectrumOperationRecord.getSingleIndex();
			mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				restoreI(i);
				for (j = 0; j < sequenceCount; j++) {
					restoreIJ(i, j);
				}
			}

			break;
		case SINGLE:

			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getSingleIndex();//AllSiteIndexs()[0];
//			mapToSrp = srpMap.getMapToSrp(k);
//			spectrumModel.getSpectrum(j).getSpectra(k).restoreState();
//			for (int i : mapToSrp) {
			for (int i : mapToSrpArray[k]){
				restoreI(i);
				restoreIJ(i, j);
			}
			
			
		case MULTI:
		
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			int restoreMulti;
			if(multiType==MultiType.BitSet){

				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					restoreI(i);
					restoreIJ(i, j);
				}
			}

			else if(multiType==MultiType.Array){

				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						restoreI(i);
						restoreIJ(i, j);
					}
				}
			}

			else if(multiType==MultiType.Hash){
				for (int i : allSrpPos) {
					restoreI(i);
					restoreIJ(i, j);
				}
				
			}
			else if(multiType==MultiType.All){
				for (int kk : siteIndexs) {
					mapToSrp = srpMap.getMapToSrp(kk);
					for (int i : mapToSrp) {
						restoreI(i);
						restoreIJ(i, j);
					}
				}

			}
			
			break;
		case RECOMBINATION:
			int[] twoSpectrumsIndex = spectrumOperationRecord.getRecombinationSpectrumIndex();

			if(multiType==MultiType.BitSet){
				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					restoreI(i);
					restoreIJ(i, twoSpectrumsIndex[0]);
					restoreIJ(i, twoSpectrumsIndex[1]);
				}
			}

			else if(multiType==MultiType.Array){

				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						restoreI(i);
						restoreIJ(i, twoSpectrumsIndex[0]);
						restoreIJ(i, twoSpectrumsIndex[1]);
					}
				}
			}

			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsHaplotypeLikelihood.class.getSimpleName() );

		}
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
//
//		
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex();
//		int siteIndex = record.getSiteIndex();
//
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(siteIndex);
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


	private void restoreI(int i) {
		eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
		sumScaledSrpLogLikelihood[i] = storedSumScaledSrpLogLikelihood[i]; 
	}


	private void restoreIJ(int i, int j) {
		allDists[i][j] = storedAllDists[i][j];
		
	}


	private void restoreEverything() {
		
		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(storedSumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);

		for (int i = 0; i < allDists.length; i++) {
			System.arraycopy(storedAllDists[i], 0, allDists[i], 0, sequenceCount);
		}
		
	}


	public double[] unittestMethodGetEachLikelihood() {
		double[] copyOfValues = new double[eachSrpLogLikelihood.length];
		System.arraycopy(eachSrpLogLikelihood, 0, copyOfValues, 0,
				copyOfValues.length);
		return copyOfValues;
	}
	
	public double calculateSrpLikelihoodFullMaster() {


//		System.out.println("calculateSrpLikelihoodMaster");
		double logLikelihood = 0;
		for (int i = 0; i < srpCount; i++) {
			
			String srp = srpMap.getSrpFragment(i);//srpArray[i]
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));

			liS.reset();
			for (int j = 0; j < sequenceCount; j++) {

				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getSequenceString());
				liS.add(logPD[dist]);
//				liS.scaleLogProb(logPD[dist]);
			}	
			
			double eachSrpLogLikelihood = liS.getLogLikelihood();
			logLikelihood += eachSrpLogLikelihood;

		}
//		 = StatUtils.sum(eachSrpLogLikelihood);
		
//for (int i = 0; i < srpCount; i++) {
//	System.out.println(Arrays.toString(allDists[i]));
//}
//System.out.println("==");

//		double logLikelihood2 = calculateSrpLikelihoodFull();
//		if(logLikelihood != logLikelihood2){
//			System.out.println(logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
//		
		
		return logLikelihood;
	}

	@Override
	public void makeDirty() {
		alignmentModel.resetOperation();
		likelihoodKnown = false;
	
	}

//
//		double logLikelihood = 0;
//		double spectrumLogLikelihood = 0;
//		double stateLogLikelihood = 0;
//		
//		for (int i = 0; i < srpCount; i++) {
//
//			String fullSrp = srpMap.getSrpFull(i);
//			int start = srpMap.getSrpStart(i);
//			int end = srpMap.getSrpEnd(i);
//			
//			liS.reset();
//			for (int j = 0; j < spectrumCount; j++) {
//
//				Haplotype haplotype = haplotypeModel.getHaplotype(j);
//				spectrumLogLikelihood = 0;
//				for (int k = start; k < end; k++) {
//					double[] frequencies = spectrum.getFrequenciesAt(k);
//					int state = getStateAtK(fullSrp, k);
//					if(state<STATE_COUNT){
////						stateLogLikelihood = caluclateStateLogLikelihood(frequencies[state]);
//						stateLogLikelihood = stateLikelihood.caluclateStateLogLikelihood(frequencies[state]);
////						double likelihood = frequencies[state] * NOT_ERROR_RATE
////								+ (1 - frequencies[state]) * ERROR_RATE;
////						stateLogLikelihood = Math.log(likelihood);
//					}
//					else{
//						stateLogLikelihood = MIN_LOG_LIKELIHOOD;
//					}
//					
//					spectrumLogLikelihood += stateLogLikelihood;
//					
//				}
//				liS.addLogProb(spectrumLogLikelihood);
//				
//
//			}	
//			logLikelihood += liS.getLogLikelihood();
//		}
//
////		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
//		if(DEBUG){
//			if(logLikelihood != this.logLikelihood){
//				System.out.println(logLikelihood +"\t"+ this.logLikelihood +"\t"+ getStoredLogLikelihood());
//			
//			
//			
//				OperationRecord record = haplotypeModel.getOperationRecord();
////				int[] siteIndexs = record.getAllSiteIndexs();
//				int j= record.getSpectrumIndex(); 
//				Haplotype haplotype = haplotypeModel.getHaplotype(j);
//				
////				stateLikelihood.calculateStatesLogLikelihood(spectra, allStateLogLikelihood[k]);
////				stateLikelihood.calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[k]);
//				
//				for (int s = 0; s < spectrumLength; s++) {
//					int k = s;
//					SpectraParameter spectra = spectrum.getSpectra(k);
////					int kOffset = k*STATE_COUNT;
//					double[] stateLn = spectra.getStateLikelihood();
////					stateLikelihood.calculateStatesLogLikelihood(spectra, 0, stateLn);
////					stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
////					System.out.println(Arrays.toString());
//					
//					double[] frequencies = spectrum.getFrequenciesAt(k);
//					double[] stateLn2 = new double[4];
//					for (int state = 0; state < 4; state++) {
//						stateLn2[state] = stateLikelihood.caluclateStateLogLikelihood(frequencies[state]);
//						if(stateLn2[state] != stateLn[state]){
//							System.out.println(s);
//							System.out.println(Arrays.toString(stateLn));
//							System.out.println(Arrays.toString(stateLn2));
//						}
//					}
////					System.out.println(Arrays.toString(stateLn2));	
//				}
////				System.exit(-1);
//			
//			
//			
//			
//			
//			
//			
//			}
//		}
//		return logLikelihood;
//	}
//

}
