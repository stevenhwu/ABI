package test.srp.spectrum.likelihood;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.SpectrumOperationRecord;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.operator.AbstractSpectrumOperator;
import srp.spectrum.operator.ColumnSpectrumDeltaExchangeOperator;
import srp.spectrum.operator.MultiSpectrumDeltaExchangeOperator;
import srp.spectrum.operator.RecombinationSpectrumOperator;
import srp.spectrum.operator.SingleSpectrumDeltaExchangeOperator;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class ShortReadsSpectrumLikelihoodTest {

	public static final double ERROR = ShortReadsSpectrumLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsSpectrumLikelihood.NOT_ERROR_RATE;
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
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
//		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double logLikelihood = likelihood.getLogLikelihood();
		double expected = Math.log(1*NOT_ERROR+0*ERROR)*8;
//		System.out.println((0.25*NOT_ERROR+0.75*ERROR) +"\t"+ Math.log(1*NOT_ERROR+0*ERROR) );
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
	}
	
	@Test
	public void testCalculateLikelihoodFixedSpectrum() {
		String[] seqs = new String[]{
				".AA",
				".AC",
				".GT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		seqs = new String[]{"AAA"};
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double[] eachLikelihood = likelihood.getEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds = new double[]{ 
				0+Math.log(1*NOT_ERROR+0*ERROR)*2,
				0+Math.log(1*NOT_ERROR+0*ERROR)*1+Math.log(0*NOT_ERROR+1*ERROR)*1,
				0+Math.log(1*NOT_ERROR+0*ERROR)*0+Math.log(0*NOT_ERROR+1*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
	}

	
	@Test
	public void testCalculateLikelihoodSpectrum() {
		String[] seqs = new String[]{
				".AA",
				".AC",
				".GT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double[] eachLikelihood = likelihood.getEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds = new double[]{ 
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		
//		logLikelihood = likelihood .getLogLikelihood();
//		expected = -0.086061253223681313806*4; //dbinom(0,8,E,log=T)
//		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);

//		double[] expecteds = new double[]{ 
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*2),
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*1)+Math.log((0.75*NOT_ERROR+0.25*ERROR)*1),
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*1)+Math.log((0.75*NOT_ERROR+0.25*ERROR)*1)
//			};

	}



	@Test
	public void testCalculateLikelihoodCustomSpectrum() {
		String[] seqs = new String[]{
				"AAAC",
				"AACT",
				"ACGT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1);
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < spectrum.getLength(); i++) {
			double[] freqs = new double[]{1-(0.1*i*3), 0.1*i, 0.1*i, 0.1*i};
			spectrum.resetFrequencies(i, freqs);
//			System.out.println("SITE: "+i +"\t"+  Arrays.toString(spectrum.getFrequencies(i)));
		}
//		spectrumModel.setSpectrum(0, spectrum);
		
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double[] eachLikelihood = likelihood.getEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		//Site1: 1, 0, 0, 0
		//Site2: 0.7, 0.1, 0.1, 0.1
		//Site3: 0.4, 0.2, 0.2, 0.2
		//Site4: 0.1, 0.3, 0.3, 0.3
		double[] expecteds = new double[] {
				// MMMD
				Math.log((1 * NOT_ERROR + 0 * ERROR)
						* (0.7 * NOT_ERROR + 0.3 * ERROR)
						* (0.4 * NOT_ERROR + 0.6 * ERROR)
						* (0.3 * NOT_ERROR + 0.7 * ERROR)),
				// MMDD
				Math.log((1 * NOT_ERROR + 0 * ERROR)
						* (0.7 * NOT_ERROR + 0.3 * ERROR)
						* (0.2 * NOT_ERROR + 0.8 * ERROR)
						* (0.3 * NOT_ERROR + 0.7 * ERROR)),

				// MDDD
				Math.log((1 * NOT_ERROR + 0 * ERROR)
						* (0.1 * NOT_ERROR + 0.9 * ERROR)
						* (0.2 * NOT_ERROR + 0.8 * ERROR)
						* (0.3 * NOT_ERROR + 0.7 * ERROR)) 
		};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		

	}
	
	@Test
	public void testFullvsSingle() throws Exception {
	
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "HaplotypeModelTest_10_srp.fasta");
		AlignmentMapping aMap = new AlignmentMapping(alignment);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		SingleSpectrumDeltaExchangeOperator op = new SingleSpectrumDeltaExchangeOperator(spectrumModel, 0.25, null);
		
		for (int i = 0; i < 1e4; i++) {
			try {
				op.doOperation();
				double logLikelihoodSingle = likelihood.getLogLikelihood();
				assertEquals(SpectrumOperation.SINGLE_DELTA, likelihood.getOperation());
				record.setOperation(SpectrumOperation.NONE);
				likelihood.makeDirty();
				double logLikelihoodFull = likelihood.getLogLikelihood();
				assertEquals(SpectrumOperation.NONE, likelihood.getOperation());
				assertEquals(logLikelihoodFull, logLikelihoodSingle, 1e-8);
				
			} catch (Exception e) {
			}
		}
	}

	private static void assertLikelihoodOperator(SpectrumAlignmentModel spectrumModel,
			AbstractSpectrumOperator op) {
		
		int ite = (int) 10;
		SpectrumOperation expectedSpectrumOperation = op.getSpectrumOperation();
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		double logLikelihoodOperator;
		double logLikelihoodFull;
		
		for (int i = 0; i < ite; i++) {
			try {
				likelihood.storeModelState();
				op.doOperation();
				likelihood.makeDirty();
				logLikelihoodOperator = likelihood.getLogLikelihood();
				assertEquals(expectedSpectrumOperation, likelihood.getOperation());
				
				SpectrumAlignmentModel spectrumModelFull = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
				ShortReadsSpectrumLikelihood likelihoodFull = new ShortReadsSpectrumLikelihood(spectrumModelFull);
				logLikelihoodFull = likelihoodFull.getLogLikelihood();
				assertEquals(SpectrumOperation.NONE, likelihoodFull.getOperation());
				assertEquals(logLikelihoodFull, logLikelihoodOperator, 1e-8);
	
				double rand = MathUtils.nextDouble();
				if(rand>0.5){
					likelihood.acceptModelState();
				}
				else{
					likelihood.restoreModelState();
				}
			} catch (OperatorFailedException e) {
			}
		}
	}

	@Test
	public void testFullvsSingleStoreRestore() throws Exception {
	
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "HaplotypeModelTest_10_srp.fasta");
		AlignmentMapping aMap = new AlignmentMapping(alignment);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
		
		SingleSpectrumDeltaExchangeOperator op = new SingleSpectrumDeltaExchangeOperator(spectrumModel, 0.25, null);
		assertLikelihoodOperator(spectrumModel, op);
	}
	

	@Test
	public void testFullvsMultiStoreRestore() throws Exception {
	
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "HaplotypeModelTest_10_srp.fasta");
		AlignmentMapping aMap = new AlignmentMapping(alignment);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
		
		MultiSpectrumDeltaExchangeOperator op = new MultiSpectrumDeltaExchangeOperator(spectrumModel, 0.1, null);
		assertLikelihoodOperator(spectrumModel, op);
	}


	@Test
	public void testFullvsColumnStoreRestore() throws Exception {
	
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "HaplotypeModelTest_10_srp.fasta");
		AlignmentMapping aMap = new AlignmentMapping(alignment);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
		ColumnSpectrumDeltaExchangeOperator op = new ColumnSpectrumDeltaExchangeOperator(spectrumModel, 0.1, null);
		assertLikelihoodOperator(spectrumModel, op);
		
	}

	

	@Test
	public void testFullvsRecombinationStoreRestore() throws Exception {
	
//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "HaplotypeModelTest_10_srp.fasta");
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "H4_srp.fasta");
		
		AlignmentMapping aMap = new AlignmentMapping(alignment);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
		
		RecombinationSpectrumOperator op = new RecombinationSpectrumOperator(spectrumModel, 100, null);
		assertLikelihoodOperator(spectrumModel, op);
	}
}
