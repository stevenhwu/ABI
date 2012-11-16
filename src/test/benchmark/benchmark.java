package test.benchmark;

import static org.junit.Assert.*;

import java.io.OutputStream;
import java.io.PrintStream;

import likelihood.LikelihoodCalculation;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import core.DataImporter;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;

public class benchmark {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		System.setErr(new PrintStream(new OutputStream() {
		    @Override
			public void write(int b) {
		    }
		}));
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
	public void tesAllLikelihoodCustomData(){
		//Take a while to run
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/benchMark/";
		
		for (int i = 0; i < 2; i++) {
			String index = "B"+i;
			String prefix = dataDir+index+"/"+index;
			System.out.print("Bencmarking"+prefix +"\n"+index +"\t");
			String trueAlignmentFile = "_true_seqs.fasta";
			String truePhylogenyFile = "_true_tree.newick";
			String shortReadFile = "_short_reads.fasta";

			DataImporter dataImporter = new DataImporter(prefix);
			SimpleAlignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
			Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
//			Sequences shortReads = dataImporter.importSequence(shortReadFile);
//			Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
			
			TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
			LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment);
			System.out.print(li.getTreeLikelihood() +"\t");
			
			li.setPopSize(3000, 0, 30000);
			System.out.print(li.getCoalescentLikelhood() +"\t");
			
			li.setPopSize(10000, 0, 1000000);
			System.out.print(li.getCoalescentLikelhood() +"\t");
			
			
			shortReadFile = "_short_reads_1.fasta";
			Sequences shortReads = dataImporter.importSequence(shortReadFile);
			li.setShortReads(shortReads);
			System.out.print(li.getShortReadLikelihood() +"\t");
			
			shortReadFile = "_short_reads_10.fasta";
			shortReads = dataImporter.importSequence(shortReadFile);
			li.setShortReads(shortReads);
			System.out.print(li.getShortReadLikelihood() +"\t");
			
			shortReadFile = "_short_reads_100.fasta";
			shortReads = dataImporter.importSequence(shortReadFile);
			li.setShortReads(shortReads);
			System.out.print(li.getShortReadLikelihood() +"\t");

			shortReadFile = "_short_reads_all.fasta";
			shortReads = dataImporter.importSequence(shortReadFile);
			li.setShortReads(shortReads);
			System.out.print(li.getShortReadLikelihood() +"\t");
			System.out.println();
		}
	}

}
