package srp.operator.haplotypes;

import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class BasesMultiOperator extends AbstractMultiOperator {


	public final static String OPERATOR_NAME = BasesMultiOperator.class.getSimpleName();
//	public final static Operation OP = Operation.SWAPMULTI;

	

	
	public BasesMultiOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		
	}


	@Override
	public String getOperatorName() {

		return OPERATOR_NAME;
	}


	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
//	    
//	    
		int[] siteIndexs = 
				generateUniqueSites(basesCount);

		for (int i = 0; i < basesCount; i++) {
			
//			SpectraParameter spectra = spectrum.getSpectra(siteIndexs[i]);
//			swapFrequency(spectra);
			char newChar = getNextBase();
			haplotype.setCharAt(i, newChar);
	        
		}
        // symmetrical move so return a zero hasting ratio
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndexs);
	
		haplotypeModel.endAlignmentModelOperation();

		return 0.0;
	}

}
