package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
//	test.dr.ext.TreeLikelihoodExtTest.class,
	test.srp.haplotypes.AAllTestsHaplotypes.class,
	test.srp.haplotypes.operator.AllTestsHaplotyesOperator.class,
	test.srp.likelihood.AllTestsLikelihood.class,
	
	})
public class AllTests {

}