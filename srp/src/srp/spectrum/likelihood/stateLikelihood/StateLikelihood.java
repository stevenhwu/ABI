package srp.spectrum.likelihood.stateLikelihood;

import srp.dr.evolution.datatype.ShortReads;
import srp.spectrum.SpectraParameter;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;

public abstract class StateLikelihood {

	private static final int STATE_COUNT = ShortReads.INSTANCE.getStateCount();
	protected static final double ERROR_RATE = ShortReadsSpectrumLikelihood.ERROR_RATE;
	protected static final double NOT_ERROR_RATE = ShortReadsSpectrumLikelihood.NOT_ERROR_RATE;
	
	public StateLikelihood() {
	}
	
	public double[] calculateStatesLogLikelihood(SpectraParameter spectra, 
			double[] statesLogLikelihood){
		for (int state = 0; state < STATE_COUNT; state++) {
			double frequency = spectra.getFrequency(state);
			statesLogLikelihood[state] = caluclateStateLogLikelihood(frequency);
		}
		
		return statesLogLikelihood;
	
	}
	
	
	public double[] calculateStoredStatesLogLikelihood(SpectraParameter spectra, 
			double[] statesLogLikelihood){
		for (int state = 0; state < STATE_COUNT; state++) {
			double frequency = spectra.getStoredFrequency(state);
			statesLogLikelihood[state] = caluclateStateLogLikelihood(frequency);
		}
		return statesLogLikelihood;
	}
	
	public abstract double caluclateStateLogLikelihood(double frequency);


	
	//XXX: beta: don't think this will work
	//XXX:: beta work!! with /300. new BetaDistribution(296.79, 3.21);
//	static BetaDistribution betaD = new BetaDistribution(296.79, 3.21);
//	static double high = betaD.logPdf(0.98)/300;
////	static double low = LOG_ERROR_RATE;// BAD
//	static double low = betaD.logPdf(0.02)/300;
	
	//XXX: MLE at mode
	//  alpha = 1.9893, beta = 1.0107
	
//	static double high = betaD.logPdf(0.98);
//	static double low = betaD.logPdf(0.02);
	
	//XXX G-test
	
//	static double high = betaD.logPdf(0.98);
//	static double low = betaD.logPdf(0.02);

}
