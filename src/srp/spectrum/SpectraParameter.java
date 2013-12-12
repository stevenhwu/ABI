package srp.spectrum;

import java.util.Arrays;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

//public class Spectra implements Parameter{
public class SpectraParameter extends Parameter.Default{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3519136356577040343L;
	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
	
	public SpectraParameter(boolean random){
		this();
		if(random){
			double[] freq = new double[4];
			for (int i = 0; i < freq.length; i++) {
				freq[i] = MathUtils.nextInt(100);
			}
			double sum = getSumOfFrequencies(freq);
			for (int i = 0; i < freq.length; i++) {
				freq[i] /= sum;
			}
			System.out.println(Arrays.toString(freq));
//			super(freq);
			setFrequenciesQuietly(freq);
			addBounds(new DefaultBounds(1.0, 0.0, getDimension()));
			
		}
	}
	
	public SpectraParameter(){
		this(EQUAL_FREQ);
		
	}

    public SpectraParameter(double[] frequencies) {
    	super(frequencies);
    	setId("spectraSSS");
    	

        double sum = getSumOfFrequencies(frequencies);
    	if(getDimension()!=4){
    		throw new IllegalArgumentException("Frequencies should have 4 elements, frequencies.length= "+getDimension());
    	}
    	//TODO make sure it's ON!!!
        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
        }
    	
		addBounds(new DefaultBounds(1.0, 0.0, getDimension()));
		if(!isWithinBounds()){
			throw new IllegalArgumentException("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies)); 
		}

    }

    private static double getSumOfFrequencies(double[] frequencies) {
        double total = 0.0;
        for (int i = 0; i < frequencies.length; i++) {
            total += frequencies[i];
        }
        return total;
    }

    public void setFrequency(int i, double value) {
    	setParameterValue(i, value);
    }
    
    //Sets the value of the parameter without firing a changed event.
    protected void setFrequenciesQuietly(double[] values){
    	for (int i = 0; i < values.length; i++) {
    		setParameterValueQuietly(i, values[i]);
		}
//    	fireParameterChangedEvent(i, Parameter.ChangeType.VALUE_CHANGED);
//    	System.arraycopy(values, 0, this.values, 0, values.length);
    }
    
    public double getFrequency(int i) {
        return getParameterValue(i);
    }

    public int getFrequencyCount() {
        return getDimension();
    }

//    public Parameter getFrequencyParameter() {
//        return spectra;
//    }

    public double[] getFrequencies() {
    	
//        double[] frequencies = new double[getFrequencyCount()];
//        for (int i = 0; i < frequencies.length; i++) {
//            frequencies[i] = getFrequency(i);
//        }
        return getParameterValues();
    }

    public double[] getCumulativeFrequencies() {
        double[] frequencies = getFrequencies();
        for (int i = 1; i < frequencies.length; i++) {
            frequencies[i] += frequencies[i - 1];
        }
        return frequencies;
    }
    
    protected void storeState() {
//    	System.err.println("storeState in Spectra");
    	super.storeValues();
	}
    protected void restoreState() {
//    	System.err.println("restoreState in Spectra");
    	super.restoreValues();
	}
    public String diagnostic(){
    	return super.diagnostic();
    }

}
