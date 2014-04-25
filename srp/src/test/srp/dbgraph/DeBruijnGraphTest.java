package test.srp.dbgraph;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.dbgraph.DeBruijnGraph;

public class DeBruijnGraphTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	private String dataDir;

	@Before
	public void setUp() throws Exception {
		dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/DBGraphData";
	}


	@After
	public void tearDown() throws Exception {
	}
	
	
	@Test
	public void testDBG() throws Exception {
//		DeBruijnImporter dbi = new DeBruijnImporter(dataDir);
//		dbi.importDeBruijnGraph("H7_56-Paired.bfast.60.cond.graph");
		
		DeBruijnGraph dbg = new DeBruijnGraph();
		dbg.addNode(0, "AA");
		dbg.addEdge(0, 1);
		dbg.addEdge(0, 2);
		dbg.addNode(1, "CCC");
		dbg.addEdge(1, 2);
		dbg.addNode(2, "GGGG");
		dbg.addEdge(2, 2);

		dbg.preprocess();
		ArrayList<String> allNodes = dbg.getAllNodes();
		ArrayList<Integer> allLength = dbg.getAllLength();
		ArrayList<ArrayList<Integer>> allEdges = dbg.getAllEdges();
		String[] expectedStrings = new String[]{"AA", "CCC", "GGGG"};
		for (int i = 0; i < dbg.getSize(); i++) {
			assertEquals(expectedStrings[i], allNodes.get(i));
			assertEquals(expectedStrings[i].length(), (int) allLength.get(i) );
		}
		assertEquals(1, (int) allEdges.get(0).get(0));
		assertEquals(2, (int) allEdges.get(0).get(1));
		assertEquals(2, (int) allEdges.get(1).get(0));
		assertEquals(2, (int) allEdges.get(2).get(0));
	}
	
}


//
//private String[] expectedsString = new String[]{
//	"TCGCACGCGGGGTAGTTGCCGTCAGCGCCTCGAGAGGTTGCCCACCCCTCTGTTTTCGCGATTTGTCTAGCGAGGCTCAATTCGTCCTTCCTATCGGATTATCCCGAACTCCGCGCGGCTAGGGTCAGGCCGCACGAGCTCT",
//	"ACGGGCACCGCTATAAAAAATCGACTCGCGTTTCGTAATTGGACCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTTGGCGCACATAATCAGGCGGCTAAAGCCGTGACCGGAAAATGGCCGGGACGGAGTTGTCCAGCCCATAAGACGGCGTATCTGAGGCTCGTAC",
//	"TCGCACGTGGGGTAGTTGCCGTCAGCGCCTCGAGAGGTTGACCACCCCCCTTTTTTCGTGGTGGGTCTAGCGAGGCTCAATTTGTCCTTCCTATCGGATTATCCCGAACCCCGCGCGGATAGGGTTAGGCCGCACGAGCTCTTTAGGTAAAAGGGGTCCCACACAACAAGTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGTCGA",
//	"ACGGGCACCGCTGTAAAGAATCGACTCGCGTTTCGTAATTGGACCACTCTTTTCTGCTGTAGTCACCTGAGGAGCTCGGCGCACATAATCAGACGGCTAAAGCCGTGACCGGAAAATTGC",
//	"CGGGCACCGCTGTAAAGAATCGACTCGCGTTTCGTAATTGGACCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTCGGCGCACATAATCA",
//	"CGGGCACCGCTTTATATAATTGACACACGTTCCGTAACTGGACCACTCTCTTATGCTGTTATCACCTGAGGAGCTTGGCGCACATAATCAGACGGCTATAGCCGTGACCGGAAAATGGCCGGGCCGGACTTGCCCAGCCTATAGGACGGCGTAACTGAGGCTCATTCCTGCTTCGCGCACGCGGGCAATTGGAGAACGTTGTAGTAGGGCGCTAGCTGGTCACGACGCGACCGAAATGTCTACCGCTATTAGTCACTCCCGCACCAATTCTAATTGGCACTAGACCAAAGCTCAAGCATAGTCTTAGGACACTCCACCGTCTTAATATACATCAGGAAGAGTGTAGCGTTCCATGGAACTGCGTTCTCAGTACGAATCGCCGTAAGTGGAATTGTAGACGTGCGCCCTCTAGATGTATACTTCTAATCGAGGGACCAAAGACAGTTGGAAGTAGTAATCAGCATCCTATTGGATAACGCGCTTGTCGATGAAGGTTTCTGACTATCTTGAAAGAGTCCACAACGTAAGTACGTGTTCCCTAGTGTGGTTACCGGTATTCTAGTTAATGGGCACGAAACTGTAGACTTATACCACAGTGTTCTCATCATGGGAGACAAAACGACTTGGGCGTTCCTGGACCATTACTTCAACAATAGATACGTTGAAATGAGTCTCGCACCATAAACCGCGGCGAACTTAGGACGAGCAAGCTGGCTTTTGGATAAATCTAGCCGGCTAAGGAGAACGATAAAATACTCAATTTACTACATGATCGGGGCCACGCGTTGTCCGCGAGCAAGGATCCGTCTTGCCTCCGGCAGCAGCCAAGTCCGGGTTATCCTGTCGTCGCGAATCTAGTGACTGTTGCAGAGTACGGCGCACGGCTGAAAGTGTCCCCTCGAACTGAGCTAGGAAGTTCGGGGTGGGCTACCCATTTGCGGGGCCTGCGCCCGCGAAGCTAACCGACAAGATTATATTGACTAGCAGCTCGACAGCGCCTCTGCTGAACGTCATTTTCATTCATCGCCACTTGTTGTGTGGGACCGCCTTGACCTAAAGAGCTCGTGCGGCCTTACCCTAGCCGTGCGGAGTTCGGGATAGTCCGATAGGAAGGAGAGATTGAGCCTCGCTAGACCAACCGCGAAAAAGGGGGTGTGGGCAACCTCTTGAGCCGTTGACAGCAACTTCCCTGCATACGA",
//	"TGGCACGCGGGGTAGATGCCGTCAGCGCCTCGAGAGGTTGACCACCCCTCCGTTTTCGCGGTTGGTCTAGCGAGGCTCAATTTGTCCTTCCTATCGGGTTATCCCCAACCCCGCGCGGGTAGGGTCAGGCCGCGCGAGCTCTCTAGGTAAAAGGGGTCCCAGACAACAAGTGGCGACGACTGAAAGTGGCCTTCAACACTTGCGCTGTCGAGCTGCTAGTCGATACAATCTTGTCGGACAGCTTCGCGGGCGCAGGCCCCGCAGATGAGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAGC",
//	"TCGCACGCGGGGTAGTTGCCGTCAGCGCCTCGAGAGGTTGCCTGCCTCTCTATTTTCGCGGTTGGTCTAGTGAGGCTCAATTTGTCCTTCCTATCGGATTATCCCAAACTCCGCGCGGCTAGGGTCAGACCGCACGAGCTCTTTAGGTAAAGGTGGTCCCGCACAACAAGTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGTCGAGCTGCTAGT",
//	"ACGGGCACCGCTGTAAAAAATCGACTCGCGTTTCGTAATCGGACCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTCGGCGCACTTAATCAGACGGCTAAAGCCGTGACCGGAAAATTGCCGGGACGGACTTGTCCAGCCCATA",
//	"ACGGGCACCGCTATAAAAAATCGACTCGCGTTTCGTAATTGGTCCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTTACCGCACATAATCAGGCGGCTAAAGCCGTGAACGGAAAATGGCCGGGACGGAGTTGTCCAGCCCATAAGACGGCGTATCTGAGGCTCGTAC",
//	"TCGTATGCAGGGAAGTTGCTGTCAACGGCTCAAGAGGTTGCCCACACCCCCTTTTTCGCGGTTGGTCTAGCGAGGCTCAATCTCTCCTTCCTATCGGACTATCCCGAACTCCGCACGGCTAGGGTAAGGCCGCACGAGCTCTTTAGGTCAAGGCGGTCCCACACAACAAGTGGCGATGAATGAAAATGACGTTCAGCAGAGGCGCTGTCGAGCTGCTAGTCAATATAATCTTGTCGGTTAGCTTCGCGGGCGCAGGCCCCGCAAATGGGTAGCCCACCCCGAACTTCCTAGCTCAGTTCGAGGGGACACTTTCAGCCGTGCGCCGTACTCTGCAACAGTCACTAGATTCGCGACGACAGGATAACCCGGACTTGGCTGCTGCCGGAGGCAAGACGGATCCTTGCTCGCGGACAACGCGTGGCCCCGATCATGTAGTAAATTGAGTATTTTATCGTTCTCCTTAGCCGGCTAGATTTATCCAAAAGCCAGCTTGCTCGTCCTAAGTTCGCCGCGGTTTATGGTGCGAGACTCATTTCAACGTATCTATTGTTGAAGTAATGGTCCAGGAACGCCCAAGTCGTTTTGTCTCCCATGATGAGAACACTGTGGTATAAGTCTACAGTTTCGTGCCCATTAACTAGAATACCGGTAACCACACTAGGGAACACGTACTTACGTTGTGGACTCTTTCAAGATAGTCAGAAACCTTCATCGACAAGCGCGTTATCCAATAGGATGCTGATTACTACTTCCAACTGTCTTTGGTCCCTCGATTAGAAGTATACATCTAGAGGGCGCACGTCTACAATTCCACTTACGGCGATTCGTACTGAGAACGCAGTTCCATGGAACGCTACACTCTTCCTGATGTATATTAAGACGGTGGAGTGTCCTAAGACTATGCTTGAGCTTTGGTCTAGTGCCAATTAGAATTGGTGCGGGAGTGACTAATAGCGGTAGACATTTCGGTCGCGTCGTGACCAGCTAGCGCCCTACTACAACGTTCTCCAATTGCCCGCGTGCGCGAAGCAGGAATGAGCCTCAGTTACGCCGTCCTATAGGCTGGGCAAGTCCGGCCCGGCCATTTTCCGGTCACGGCTATAGCCGTCTGATTATGTGCGCCAAGCTCCTCAGGTGATAACAGCATAAGAGAGTGGTCCAGTTACGGAACGTGTGTCAATTATATAAAGCGGTGCCCG",
//	"ACGGGCACCGCTGTAAAAAATCGACTCCCGTATCGTAATTGGACCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTCGGCGCACATAATCA",
//	"CGCACGCGGGGTAGTTGCCGTCAGCGCCTCGAGAGGTTGACCACCCCTCTGTTTTCGCGGTTGGTCTAGCGAGGCTCAATTTGTCCTTCCTATCGGATT",
//	"TCGCACGTGGGGTAGTTGCCGTCAGCGCCTCGAGAGGTTGCCCACCCCTCTGTTTTCGCGGTTGGTCTAGCGAGGCTCAATTTGTCCTTCCTATCGGATT",
//	"GTCCTTCCTATCGGATTATCCCGAACTCCGCGCGGCTAGGGTCAGGCCGCACGAGCTCTTTA",
//	"CGGAAAATGGCCGGGACGGAGTTGTCCAGCCCATAAGACGGCGTATCTGAGGCTCGTACC",
//	"GGGGTCCCACACAACAAGTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGTCGAGCTGCTAGT",
//	"GTCACCTGAGGAGCTCGGCGCACATAATCAGACGGCTAAAGCCGTGACCGGAAAATTGCCGGGA",
//	"TCGTAATTGGACCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTCGGCGCACATAATCAG",
//	"ATGAGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAGCCGT",
//	"CACAACAAGTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGTCGAGCTGCTAGTCGATACAA",
//	"TAATCAGACGGCTAAAGCCGTGACCGGAAAATTGCCGGGACGGACTTGTCCAGCCCATAA",
//	"CCACCCCTCTGTTTTCGCGGTTGGTCTAGCGAGGCTCAATTTGTCCTTCCTATCGGATTATCCCGAAC",
//	"CTTCCTATCGGATTATCCCGAACTCCGCGCGGCTAGGGTCAGGCCGCACGAGCTCTTTAGGTAAAGGGGGTCCCACACAACAAGTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGTCGA",
//	"CTTCCTATCGGATTATCCCGAACTCCGCGCGGCTAGGGTCAGGCCGCACGAGCTCTTTAAGTAAAAGGGCTCCCACACAATAAGTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGCCGAGCAGCTAGGCGATACAATCTTGTCAGACAGCTTCGCGGGCGCAGGCCCCGCAAATTAGTAGCCCACCTTAAACTTCACAACTCAGTTCGAGCGAACACGTTTAGGCGTGAGCCATACTCTGCAACAGTCACTAGATTCGCTAAGAGAGGATAACCCGGGCTAGGCCGTTGCCGGAGCCAAGACAGAACCTTGCTCTCGCACAACGCGTGGTCCCGAACATGTATTGGATCGAGTAATTTATCGGTCTCTCTATCCGGCCGGACTTATCCAAATGCCAGCTTGCTCGTCCTAAGCTCGCCGCGGTATATGGTGCGAGACTTGTCTCGACGTATTCATTGTTCAGACCATGGTCCGGGGTCGCCCAAGTCGTTTTATGTCCCATCATGAGAATAAAGTGGTATGAGTTTCCAGTTTCGTGCCAAGTAACCAGAATTCTGGTGGCCACACTAGGGGACACGTAATTAAGTTGTGGTCTCTTCCCAGTTGGTCAGAGACCTTCATCGATCAGCCCGCTATCCAATAGAATGCCGATTGTTACTTCCAACTATCTTTGGTCCCCCGGTTTAAAATATACATCAAGAGGTTGTACGTCTACGATTTCACATACGGGGCTTCGCATTGGGAACGCAGCTCCCTGGATCGGTACGCCCTTCCCAACGTATACGGAGACGGTGGGGTGTCCTAAGACTATTCTTGAGCTTTAGTCTAGAGCCAATTGGAATTGAAGCGAGAGTGACAGATAACAGGGGAGATTTCGGTCGCATTGTGGCCGGCTCGTGCCCTACTACAACGTTCTCCAATTGCCCGGGTGAGCCGAGCGGGTACGAGCCTCAGACACGGCGTGTTATGGGCTGGACAAGTCCGTCCCGGCCATTTTCCGGTCACGCGTTTAGCCGCCTGATTATGTGCGCCGAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGTCCAATTACG",
//	"GGAAAATGGCCGGGACGGAGTTGTCCAGCCCATAAGACGGCGTATCTGAGGCTCGTACCGGCTTCGCTCACGCGGGCAATTGGAGAACGTTGTAGTAGGGCACGAGCCGGCCACAATGCGACCGAAATCTCTCCCGTTATCTGTCACTCTCGCATCA",
//	"GGAAAATGGCCGGGACGGAGTTGTCCAGCCCATAAGACGGCGTATCTGAGGCTCGTACCCGCTTCGCTCACGCGGGCAATTGGAGAACGTTGTAGTAAGGCACGAGCCGGCCACAATGCGACCGAAATCTCTCCCGTTATCTGTCACTCTCGCATCA",
//	"CTGAGGAGCTCGGCGCACATAATCAGACGGCTAAAGCCGTGACCGGAAAATTGCCGGGATGGACTTGTCCAGCCCATAATACGCCGTGTCTGAGGCTCGTACCCGCTCCGCTCACGCGG",
//	"CTGAGGAGCTCGGCGCACATAATCAGACGGCTAAAGCCGTGACCGGAAAATTGCCGGGACGGACTTGTCCAGCCCATA",
//	"CGTAATTGGACCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTCGGCGCACATAATCAGGCGGCTAAACGCGTGACCGGAAAATGGCCGGGACGGACTTGTCCAGCCCATAACACGCCGTGTCTGAGGCTCGTACCCGCTCGGCTCACCCGGGCAATTGGAGAACGTTGTAGTAGGGCACGAGCCGGCCACAATGCGACCGAAATCTCCCCTGTTATCTGTCACTCTCGCTTCAATTCCAATTGGCTCTAGACTAAAGCTCAAGAATAGTCTTAGGACACCCCACCGTCTCCGTATACGTTGGGAAGGGCGTACCGATCCAGGGAGCTGCGTTCCCAATGCGAAGCCCCGTATGTGAAATCGTAGACGTACAACCTCTTGATGTATATTTTAAACCGGGGGACCAAAGATAGTTGGAAGTAACAATCGGCATTCTATTGGATAGCGGGCTGATCGATGAAGGTCTCTGACCAACTGGGAAGAGACCACAACTTAATTACGTGTCCCCTAGTGTGGCCACCAGAATTCTGGTTACTTGGCACGAAACTGGAAACTCATACCACTTTATTCTCATGATGGGACATAAAACGACTTGGGCGACCCCGGACCATGGTCTGAACAATGAATACGTCGAGACAAGTCTCGCACCATATACCGCGGCGAGCTTAGGACGAGCAAGCTGGCATTTGGATAAGTCCGGCCGGATAGAGAGACCGATAAATTACTCGATCCAATACATGTTCGGGACCACGCGTTGTGCGAGAGCAAGGTTCTGTCTTGGCTCCGGCAACGGCCTAGCCCGGGTTATCCTCTCTTAGCGAATCTAGTGACTGTTGCAGAGTATGGCTCACGCCTAAACGTGTTCGCTCGAACTGAGTTGTGAAGTTTAAGGTGGGCTACTAATTTGCGGGGCCTGCGCCCGCGAAGCTGTCTGACAAGATTGTATCGCCTAGCTGCTCGGCAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCACTTATTGTGTGGGAGCCCTTTTACTTAAAGAGCTCGTGCGGCCTGACCCTAGCCGCGCGGAGTTCGGGATAATCCGATAGGAAG",
//	"CGTAATTGGACCACTCTTTTCTGCTGTTGTCACCTGAGGAGCTCGGCGCACATAATCAGACGGCTAAAGCCGTGACCGGAAAATTGC",
//	"AGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAGCCGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATAACCC",
//	"AGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAGCCGTGCTCTGCAACAGTCACTAGATTCGCGATGAGAGAATAACCCGGACTAGGCCGTTGCCAGAGGCAGGACAGAACCTTGCTCTCGCATAACACGTGGCCCCAAACATGTATTGGATCGAGTATTTTATCGGCCTCCCTATTCGGCCGGACTTATCCAAATGCCAGCTTGCTCGTCCGAAGCTCGCCGCGGTTTATGGTGCGAGAC",
//	"GTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGTCGAGCTGCTAGTCGATACAACCTTGTCGGACAGCTTCGCGGGCGCAGGCCCCGCAAATGAGTAGCCCACCCTAAACTTCA",
//	"GTGGCGACGACTGAAAGTGGCCTTCACCAGTTGCGCTGTCGAGCTGCTAGTCGATACAATCTTGTCGGTCAGCTTC",
//	"AATCAGACGGCTAAAGCCGTGACCGGAAAATTGCCGGGACGGACTTGTCCAGCCCATAAAACGCCGTGTCTGAGGCTCGTACCCGCTCCGCTCACGCGGGCAATTGGAGAACGTTGTAGTAGGGCACGAGCCGACCACAATGCGACCGAAATCGCTCCCGTTATCTGTCACTCTCGCGTCAATTCCAATTGGCTCTAGATCAAAGCCCAAGCATAGTCTTAAGACACCCCACCGTCTCCGTATACATTGGGAAGGGCGTACCGGTCCATGGAGCTGTGTTCCCAATGCGGAGCCCCGTATGTGAAATAGTAGACGTACACCCTCTTGATGTATATTTTAAACCGGGGGACCAAAGATAGTTGGGAGTAACAAT",
//	"AATCAGACGGCTAAAGCCGTGACCGGAAAATTGCCGGGACGGACTTGTCCAGCCCATAATACGCCGTGTCTGAGGCTCGTACCCGCTCCGCTCACGCGG",
//	"TGTTTTCGCGGTTGGTCTAGCGAGGCTCAATTTGTCCTTCCTATCGGATTATCCCGAACTCCGCGCGGCTAGGGTCAGGCCGCACGAGCTCT",
//	"TGTTTTCGCGGTTGGTCTAGCGAGGCTCAATTTGTCCTTCCTATCGGATTATCCCGAACCCCGCGCGGCTAAGGTCAGGCCGCACGAGCCCTTTAGGTAAAAGGGGTCCCACACAACAAGTGGCGACGAGTGCAAGTGGCCTTCACCAGTTGCGCTGTCGGGCTGCTTGTCGATACAATCTTGTCGGACAGCTTCGCGGGCGCAGGCCCCGCAAATGAGTAGCCCACCCTAAACTTCA",
//	"CTGATTATGTGCGCCGAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGTCCAATTACGA",
//	"GGCACGAGCCGGCCACAATGCGACCGAAATCTCTCCCGTTATCTGTCACTCTCGCATCAA",
//	"GGACTTGTCCAGCCCATAATACGCCGTGTCTGAGGCTCGTACCCGCTCCGCTCACGCGGGC",
//	"TAAAGAGCTCGTGCGGCCTGACCCTAGCCGCGCGGAGTTCGGGATAATCCGATAGGAAGGAC",
//	"CTTTTAGGCGTGAGCCGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATAACCCGGGCTAGGCCGTTGCCGGAGGCAAGACAGAACCTTGCTCTCGCATAACACGTGGTCCCGAACATGTATTGGATCGAGTATTTTATCGGTTTCCCTATTCGGCCGGACTTATCCAAATGCCAGCTTGCTCGTCCGAAGCTCGCCGCAGTTTATGGTGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCGGGGTCGC",
//	"CTTTTAGGCGTGAGCCGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATAACCCTGGCTAGGCCGTTGCCGGAGGCAAGACAGAACCTTGTTCTCGCACAACACGTGGTCCCGAACACGTATTGGATCGAGTATTTTATCGATCTCCCTATCCGACTGGACTTATCCAAATGCCAGCTTGCTCGTCCGAAGCTCGCCGCGGTTTATGGTGCGAGAC",
//	"GGACTTATCCAAATGCCAGCTTGCTCGTCCGAAGCTCGCCGCGGTTTATGGTGCGAGACTTGTCTCAACGTATT",
//	"CTTGTCGGACAGCTTCGCGGGCGCAGGCCCCGCAAATGAGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAG",
//	"TGGCCTTCACCAGTTGCGCTGTCGAGCTGCTAGTCGATACAATCTTGTCGGTCAGCTTCGCGGGCGCGGACCCCGCAAATGAGTAGCCCACCCTAAATTTCACAACTCAGTTCGTGCAAACACTTTTAGGCGTGAGCTGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATA",
//	"TGGCCTTCACCAGTTGCGCTGTCGAGCTGCTAGTCGATACAATCTTGTCGGTCAGCTTCACGGGCGCGGGCCCCGCGAATGAGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAG",
//	"GTACACCCTCTTGATGTATATTTTAAACCGGGGGACCAAAGATAGTTGGGAGTAACAATCGGCA",
//	"TGATTATGTGCGCCGAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGTCCAATTACGATACGGGAGTCGATTTTTTACAGCGGTGCCCGT",
//	"TGATTATGTGCGCCGAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGTCCAATTACGAAACGCGAGTCGATTCTTTACAGCGGTGCCCG",
//	"GCACGAGCCGGCCACAATGCGACCGAAATCTCTCCCGTTATCTGTCACTCTCGCATCAATTCCAACTGGCTCTAGACCAAAGCTCAAGGATAGTCTTCAGGCACCCCACCGTCCCCGTATAAATTAGGAAAGGCGTACCGATCCACGGAGCTGAGTTCCCAATGCGAAGCCCCGTAAGTGGAATCGTAGACGAACACCCTCTTGATGTATATTTTAAATCGGGGGACCATAGATAGTTGGAAGTAACAATCGGCATTCTATTGGATAACGGCCTGGTCGATGAGAGTCTCTGACCAACCGGGAAGAGACCACAATTTAATTACGTGTCCCCTAGTGCGGTCACCAGAATTCTAGTTACGTGGCACGAAACTGTAGACTCATACCACTTAATTCTCA",
//	"GCACGAGCCGGCCACAATGCGACCGAAATCTCTCCCGTTATCTGTCACTCTCGCATCAACTCCAATTGGCTCTAGACCAAAGCTCAAGGATAGTCTTAAGGCACCCCACCGTCCCCGTATACATTAGGAAAGGCGCACCGATCCACGGAGCTGAGTTCCCAATGCGAAGCCCCGTATGTGGACTCGTAGACGTACACCCTCTTGATGGATATTTTAAATCGGGGGACCAAAGATAGTTGGAAGTAACAATCGGCATTCTATTGGATAACGGGCTGGTCGGTGAGAGTCTCTGACCAGCCGGGAATAGACCACAACTTAATTACGTGTCCCCTAGTGTGGTCACCAGAATTCTAGTTACGTGGCACGAAACTGTAGACTCATACCACTTAATTCTCA",
//	"ACTTGTCCAGCCCATAATACGCCGTGTCTGAGGCTCGTACCCGCTCCGCTCACGCGGGCGATTGGAGAACGTTGTAGTAGGGCACGAGCCGGCCACAATGCGACTGAAATCTCTCCCGTTATGTGTCACGCTCACATCAATTCCAATTGGCTCTAGACCAAAGCTCAAGCATATTCTTAAGACACCCCACCGTCTCCGTATGCATTAGGAAGGGCGTACCGATCCATGGAGCTGTGTTCCCAATGCGAAGCCCCGTATGTGAAATCATAGATGTACACCCTCTTGATGTATATTTTAAACCGGGGGACCAAAGATAGTTGGGAGTAACAAT",
//	"ACTTGTCCAGCCCATAATACGCCGTGTCTGAGGCTCGTACCCGCTCCGCTCACGCGGGCAATTGGAGAATGTTGTCGTAGGGCACGAGCCGGCCACAATGCGACCGAAATCTCTCCCGTTATGTGTCACTCTCGCATCAATTCCAATTGGCTCTAGACCAAAGCTCAAGCATACTCTTAAGACACCCCACCGTCTCCGTATACATTAGAAAGGGCGTACCGATCCATGGAGCTGAGTTCCCAATGCGAAGCCCCGTATGTGAAATCGTAGACGTACGCCCTCTTGATGTATATTTTAAACCGGGGGACCAAAGATAGTTGGGAGTAACAATCGGCA",
//	"AGAGCTCGTGCGGCCTGACCCTAGCCGCGCGGAGTTCGGGATAATCCGATAGGAAGGACAAATTGAGCCTCGCTAGACCAACCGCGAAAACA",
//	"AGAGCTCGTGCGGCCTGACCCTAGCCGCGCGGAGTTCGGGATAATCCGATAGGAAGGACGAATTGAGCCTCGCTAGACAAATCGCGAAAACAGAGGGGTGGGCAACCTCTCGAGGCGCTGACGGCAACTACCCCGCGTGCGA",
//	"GTTTATGGTGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCGGGGTCGCC",
//	"CCAGCTTGCTCGTCCGAAGCTCGCCGCGGTTTATGGTGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCGGGGTCGC",
//	"CCAGCTTGCTCGTCCGAAGCTCGCCGCGGTTTATGGTGCGAGACTTGTCTCAACGTATTAATTGTTCAGACCATGGTCCGGGGTCGCCCAAGTCTTTTTATGTCCCATCATGAGAGTAAAGTGGTATGAGTCTACAGTTTCGTGTCAAGTAACTAGAAATCTGGTGACCACACTAGGGGACACGTAATTAAGTTGTGGTCTCTCCCCAGTTGGCCAGGGACTCTCACCGGCCAGCCCGTTATCCAATAGAATGCCGATTGTTACTCCCAACTATCTTTGGTCCCCC",
//	"AATGAGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAGC",
//	"AACACTTTTAGGCGTGAGCTGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATAACCCGGG",
//	"CCCTCTTGATGTATATTTTAAACCGGGGGACCAAAGATAGTTGGGAGTAACAATCGGCATTCTATTGGATAACGGGCTGGCCG",
//	"GGTCACCAGAATTCTAGTTACGTGGCACGAAACTGTAGACTCATACCACTTAATTCTCAT",
//	"GTTCGGGATAATCCGATAGGAAGGACAAATTGAGCCTCGCTAGACCAACCGCGAAAACAGAGGGGTGG",
//	"TTTATGGTGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCGGGGTCGCCCAAGTCATTTTATGTCCCATCATGAGAGTAGAGTGGTATGAGTCTACAGTTTCGTGCCAAGTAACTAGAGATCTGGTGACCACACTAGGGGACACGTAATTAAGTTGTGGTCTCTTCCCAGTTAGCCAGAGACCCTCATCGGCCAGCCCGTTATCCAATAGAATGCCGATTGTTACTCCCAACTATCTTTGGTCCC",
//	"TTTATGGTGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCGGGGTCGCCTGGGTCGTTTTATGTCCCATGATGAGAGTAAAATGGTATGAGTCTACAGTTTCGTGCCAAGTAACTAGAAATCCGGTGACCACACTAGGGGACGCGTAATTAAGTTGTGGTCTCTTCCCAGTTGGCCAGAGACCCTTATCGGCCAGCCCGTTATCCAATAGAATGCCGATTGTTACTCCCAACTATCTTTGGTCCC",
//	"CGGCCAGCCCGTTATCCAATAGAATGCCGATTGTTACTCCCAACTATCTTTGGTCCCCCGGTTTAAAATATACATCAAGAGGG",
//	"ATGAGTAGCCCACCCTAAACTTCACAACTCAGTTCGAGCGAACACTTTTAGGCGTGAGCTGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATA",
//	"TTAGGCGTGAGCTGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATAACCCGGGTTAGGCCGCTGCCGGAGGCAAGACAGAACCTTGTTCGTGCATGACCCGTGGTCCCGAGCATGTATTGGGTCGAGTATTTTATCGGTCTCCCTATCCGGCCGAACTTATCCAAATGCCTGCTTGCTCGTCCTAAGCTTGCCGCGGTTTATGATGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCTG",
//	"TTAGGCGTGAGCTGTACTCTGCAACAGTCACTAGATTCGCGATGAGAGGATAACCCGGGCTAGGCCGCTGCCGGAAGCAAGACAGAACCTTGCTCGTGCATAACCCGTGGTCCCGAACATGCATTGGGTCGAGTATTTTAGCGGTCTCCCTATCCGGCCGGACTTATCCAAATGCCAGCTTGCTCGTCCTAAGGTCGCCGCGGTTTATGATGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCTG",
//	"GGGGGACCAAAGATAGTTGGGAGTAACAATCGGCATTCTATTGGATAACGGGCTGGCCGAT",
//	"GGGGGACCAAAGATAGTTGGGAGTAACAATCGGCATTCTATTGGATAACGGGCTGGCCGGTGAGAGTCCCTGGCCAACTGGGGAGAGACCACAACTTAATTACGTGTCCCCTAGTGTGGTCACCAGATTTCTAGTTACTTGACACGAAACTGTAGACTCATACCACTTTACTCTCATGATGGGACATAAAAAGACTTGGGCGACCCCGGACCATGGTCTGAACAATTAATACGTTGAGACAAGTCTCGCACCATAAACCGCGGCGAGCTTCGGACGAGCAAGCTGG",
//	"GTCACCAGAATTCTAGTTACGTGGCACGAAACTGTAGACTCATACCACTTAATTCTCATGATGGGACATCAAACGACTTGGGCAACCCAGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCATCATAAACCGC",
//	"GTCACCAGAATTCTAGTTACGTGGCACGAAACTGTAGACTCATACCACTTAATTCTCATAATGGGACATCAAACGACTTGGGCGACCCAGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCATCATAAACCGC",
//	"AATCCGATAGGAAGGACAAATTGAGCCTCGCTAGACCAACCGCGAAAACAGAGGGGTGGGCAACCTCTCGAGGCGCTGACGGCAACTACCCCACGTGCGA",
//	"AATCCGATAGGAAGGACAAATTGAGCCTCGCTAGACCAACCGCGAAAACAGAGGGGTGGTCAACCTCTCGAGGCGCTGACGGCAACTACCCCGCGTGCG",
//	"ATCGGCCAGCCCGTTATCCAATAGAATGCCGATTGTTACTCCCAACTATCTTTGGTCCCCC",
//	"TGCCGATTGTTACTCCCAACTATCTTTGGTCCCCCGGTTTAAAATATACATCAAGAGGGCGTACGTCTACGATTTCACATACGGGGCTTCGCATTGGGAACTCAGCTCCATGGATCGGTACGCCCTTTCTAATGTATACGGAGACGGTGGGGTGTCTTAAGAGTATGCTTGAGCTTTGGTCTAGAGCCAATTGGAATTGATGCGAGAGTGACACATAACGGGAGAGATTTCGGTCGCATTGTGGCCGGCTCGTGCCCTACGACAACATTCTCCAATTGCCCGCGTGAGCGGAGCGGGTACGAGCCTCAGACACGGCGTATTATGGGCTGGACAAGT",
//	"TGCCGATTGTTACTCCCAACTATCTTTGGTCCCCCGGTTTAAAATATACATCAAGAGGGTGTAC",
//	"GCCGCGGTTTATGATGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCTGGGT",
//	"GGGACCAAAGATAGTTGGGAGTAACAATCGGCATTCTATTGGATAACGGGCTGGCCGATAAGGGTCTCTGGCCAACTGGGAAGAGACCACAACTTAATTACGCGTCCCCTAGTGTGGTCACCGGATTTCTAGTTACTTGGCACGAAACTGTAGACTCATACCATTTTACTCTCATCATGGGACATAAAACGACCCAGGCGACCCCGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCACCATAAA",
//	"GGGACCAAAGATAGTTGGGAGTAACAATCGGCATTCTATTGGATAACGGGCTGGCCGATGAGGGTCTCTGGCTAACTGGGAAGAGACCACAACTTAATTACGTGTCCCCTAGTGTGGTCACCAGATCTCTAGTTACTTGGCACGAAACTGTAGACTCATACCACTCTACTCTCATGATGGGACATAAAATGACTTGGGCGACCCCGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCACCATAAA",
//	"AATACGTTGAGACAAGTCTCGCACCATAAACCGCGGCGAGCTTCGGACGAGCAAGCTGGCATTTGGATAAGTCC",
//	"ACCCAGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCATCATAAACCGCGGC",
//	"GCCCGCGTGAGCGGAGCGGGTACGAGCCTCAGACACGGCGTATTATGGGCTGGACAAGTCC",
//	"ATTGTTACTCCCAACTATCTTTGGTCCCCCGGTTTAAAATATACATCAAGAGGGTGTACATCTATGATTTCACATACGGGGCTTCGCATTGGGAACACAGCTCCATGGATCGGTACGCCCTTCCTAATGCATACGGAGACGGTGGGGTGTCTTAAGAATATGCTTGAGCTTTGGTCTAGAGCCAATTGGAATTGATGTGAGCGTGACACATAACGGGAGAGATTTCAGTCGCATTGTGGCCGGCTCGTGCCCTACTACAACGTTCTCCAATCGCCCGCGTGAGCGGAGCGGGTACGAGCCTCAGACACGGCGTATTATGGGCTGGACAAGT",
//	"ATTGTTACTCCCAACTATCTTTGGTCCCCCGGTTTAAAATATACATCAAGAGGGTGTACGTCTACTATTTCACATACGGGGCTCCGCATTGGGAACACAGCTCCATGGACCGGTACGCCCTTCCCAATGTATACGGAGACGGTGGGGTGTCTTAAGACTATGCTTGGGCTTTGATCTAGAGCCAATTGGAATTGACGCGAGAGTGACAGATAACGGGAGCGATTTCGGTCGCATTGTGGTCGGCTCGTGCCCTACTACAACGTTCTCCAATTGCCCGCGTGAGCGGAGCGGGTACGAGCCTCAGACACGGCGTTTTATGGGCTGGACAAGTCCGTCCCGGCAATTTTCCGGTCACGGCTTTAGCCGTCTGATT",
//	"GCGGTTTATGATGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCTGGGTTGCCCAAGTCGTTTGATGTCCCATCATGAGAATTAAGTGGTATGAGTCTACAGTTTCGTGCCACGTAACTAGAATTCTGGTGAC",
//	"GCGGTTTATGATGCGAGACTTGTCTCAACGTATTCATTGTTCAGACCATGGTCCTGGGTCGCCCAAGTCGTTTGATGTCCCATTATGAGAATTAAGTGGTATGAGTCTACAGTTTCGTGCCACGTAACTAGAATTCTGGTGAC",
//	"GGCGACCCCGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCACCATAAAC",
//	"GTCTCGCACCATAAACCGCGGCGAGCTTCGGACGAGCAAGCTGGCATTTGGATAAGTCCAGTCGGATAGGGAGATCGATAAAATACTCGATCCAATACGTGTTCGGGACCACGTGTTGTGCGAGAACAAGGTTCTGTCTTGCCTCCGGCAACGGCCTAGCCAGGGTTATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACGGCTCACGCCTAAAAG",
//	"GTCTCGCACCATAAACCGCGGCGAGCTTCGGACGAGCAAGCTGGCATTTGGATAAGTCCGGCCGAATAGGGAGGCCGATAAAATACTCGATCCAATACATGTTTGGGGCCACGTGTTATGCGAGAGCAAGGTTCTGTCCTGCCTCTGGCAACGGCCTAGTCCGGGTTATTCTCTCATCGCGAATCTAGTGACTGTTGCAGAGCACGGCTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACT",
//	"CAGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCATCATAAACCGCGGCGACCTTAGGACGAGCAAGCTGGCATTTGGATAAGTCCGGCCGGATAGGGAGACCGCTAAAATACTCGACCCAATGCATGTTCGGGACCACGGGTTATGCACGAGCAAGGTTCTGTCTTGCTTCCGGCAGCGGCCTAGCCCGGGTTATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACAGCTCACGCCTAA",
//	"CAGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCATCATAAACCGCGGCAAGCTTAGGACGAGCAAGCAGGCATTTGGATAAGTTCGGCCGGATAGGGAGACCGATAAAATACTCGACCCAATACATGCTCGGGACCACGGGTCATGCACGAACAAGGTTCTGTCTTGCCTCCGGCAGCGGCCTAACCCGGGTTATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACAGCTCACGCCTAA",
//	"CCGCGTGAGCGGAGCGGGTACGAGCCTCAGACACGGCGTATTATGGGCTGGACAAGTCCATCCCGGCAATTTTCCGGTCACGGCTTTAGCCGTCTGATTATGTGCGCCGAGCTCCTCAG",
//	"CCGCGTGAGCGGAGCGGGTACGAGCCTCAGACACGGCGTATTATGGGCTGGACAAGTCCGTCCCGGCAATTTTCCGGTCACGGCTTTAGCCGTCTGATT",
//	"TTATGGGCTGGACAAGTCCGTCCCGGCAATTTTCCGGTCACGGCTTTAGCCGTCTGATTA",
//	"ATGAGAATTAAGTGGTATGAGTCTACAGTTTCGTGCCACGTAACTAGAATTCTGGTGACC",
//	"GCGACCCCGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCACCATAAACCGCGGCGAGCTTCGGACGAGCAAGCTGG",
//	"GCGACCCCGGACCATGGTCTGAACAATGAATACGTTGAGACAAGTCTCGCACCATAAACTGCGGCGAGCTTCGGACGAGCAAGCTGGCATTTGGATAAGTCCGGCCGAATAGGGAAACCGATAAAATACTCGATCCAATACATGTTCGGGACCACGTGTTATGCGAGAGCAAGGTTCTGTCTTGCCTCCGGCAACGGCCTAGCCCGGGTTATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACGGCTCACGCCTAAAAG",
//	"GGGTTATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACGGCTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACT",
//	"ACGGCTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACTCAT",
//	"CCCGGGTTATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACAGCTCACGCCTAAAAGTGTT",
//	"TCCCGGCAATTTTCCGGTCACGGCTTTAGCCGTCTGATTATGTGCGCCGAGCTCCTCAGGTGAC",
//	"TATGGGCTGGACAAGTCCGTCCCGGCAATTTTCCGGTCACGGCTTTAGCCGTCTGATTAAGTGCGCCGAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGTCCGATTACGAAACGCGAGTCGATTTTTTACAGCGGTGCCCGT",
//	"TATGGGCTGGACAAGTCCGTCCCGGCAATTTTCCGGTCACGGCTTTAGCCGTCTGATTATGTGCGCCGAGCTCCTCAG",
//	"TGAGAATTAAGTGGTATGAGTCTACAGTTTCGTGCCACGTAACTAGAATTCTGGTGACCGCACTAGGGGACACGTAATTAAATTGTGGTCTCTTCCCGGTTGGTCAGAGACTCTCATCGACCAGGCCGTTATCCAATAGAATGCCGATTGTTACTTCCAACTATCTATGGTCCCCCGATTTAAAATATACATCAAGAGGGTGTTCGTCTACGATTCCACTTACGGGGCTTCGCATTGGGAACTCAGCTCCGTGGATCGGTACGCCTTTCCTAATTTATACGGGGACGGTGGGGTGCCTGAAGACTATCCTTGAGCTTTGGTCTAGAGCCAGTTGGAATTGATGCGAGAGTGACAGATAACGGGAGAGATTTCGGTCGCATTGTGGCCGGCTCGTGC",
//	"TGAGAATTAAGTGGTATGAGTCTACAGTTTCGTGCCACGTAACTAGAATTCTGGTGACCACACTAGGGGACACGTAATTAAGTTGTGGTCTATTCCCGGCTGGTCAGAGACTCTCACCGACCAGCCCGTTATCCAATAGAATGCCGATTGTTACTTCCAACTATCTTTGGTCCCCCGATTTAAAATATCCATCAAGAGGGTGTACGTCTACGAGTCCACATACGGGGCTTCGCATTGGGAACTCAGCTCCGTGGATCGGTGCGCCTTTCCTAATGTATACGGGGACGGTGGGGTGCCTTAAGACTATCCTTGAGCTTTGGTCTAGAGCCAATTGGAGTTGATGCGAGAGTGACAGATAACGGGAGAGATTTCGGTCGCATTGTGGCCGGCTCGTGC",
//	"GCTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACTCATCTGCGGGGCCTGCGCCCGCGAAGCTGTCCGACAAGATTGTATCGACTAGCAGCTCGACAGCGCAAGTGTTGAAGGCCACTTTCAGTCGTCGCCACTTGTTGTCTGGGACCCCTTTTACCTAGAGAGCTCGCGCGGCCTGACCCTACCCGCGCGGGGTTGGGGATAACCCGATAGGAAGGACAAATTGAGCCTCGCTAGACCAACCGCGAAAACGGAGGGGTGGTCAACCTCTCGAGGCGCTGACGGCATCTACCCCGCGTGCCA",
//	"GCTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACTCATT",
//	"TATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACAGCTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACTCAT",
//	"TATCCTCTCATCGCGAATCTAGTGACTGTTGCAGAGTACAGCTCACGCCTAAAAGTGTTTGCACGAACTGAGTTGTGAAATTTAGGGTGGGCTACTCATTTGCGGGGTCCGCGCCCGCGAAGCTGACCGACAAGATTGTATCGACTAGCAGCTCGACAGCGCAACTGGTGAAGGCCA",
//	"GCAATTTTCCGGTCACGGCTTTAGCCGTCTGATTATGTGCGCCGAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGTCCAATTACG",
//	"GCAATTTTCCGGTCACGGCTTTAGCCGTCTGATTATGTGCGCCGAGCTCCTCAGGTGACTACAGCAGAAAAGAGTGGTCCAATTACGAAACGCGAGTCGATTCTTTACAGCGGTGCCCGT",
//	"TTGATGCGAGAGTGACAGATAACGGGAGAGATTTCGGTCGCATTGTGGCCGGCTCGTGCC",
//	"CTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACTCATTCGCGGGGCCCGCGCCCGTGAAGCTGACCGACAAGATTGTATCGACTAGCAGCTCGACAGCGCAACTGGTGAAGGCCA",
//	"CTCACGCCTAAAAGTGTTCGCTCGAACTGAGTTGTGAAGTTTAGGGTGGGCTACTCATTTGCGGGGCCTGCGCCCGCGAAGCTGTCCGACAAG",
//	"GAAGCTGACCGACAAGATTGTATCGACTAGCAGCTCGACAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCAC",
//	"TGATGCGAGAGTGACAGATAACGGGAGAGATTTCGGTCGCATTGTGGCCGGCTCGTGCCTTACTACAACGTTCTCCAATTGCCCGCGTGAGCGAAGCGGGTACGAGCCTCAGATACGCCGTCTTATGGGCTGGACAACTCCGTCCCGGCCATTTTCC",
//	"TGATGCGAGAGTGACAGATAACGGGAGAGATTTCGGTCGCATTGTGGCCGGCTCGTGCCCTACTACAACGTTCTCCAATTGCCCGCGTGAGCGAAGCCGGTACGAGCCTCAGATACGCCGTCTTATGGGCTGGACAACTCCGTCCCGGCCATTTTCC",
//	"TGAAGTTTAGGGTGGGCTACTCATTTGCGGGGCCTGCGCCCGCGAAGCTGTCCGACAAGATTGTATCGACAAGCAGCCCGACAGCGCAACTGGTGAAGGCCACTTGCACTCGTCGCCACTTGTTGTGTGGGACCCCTTTTACCTAAAGGGCTCGTGCGGCCTGACCTTAGCCGCGCGGGGTTCGGGATAATCCGATAGGAAGGACAAATTGAGCCTCGCTAGACCAACCGCGAAAACA",
//	"TGAAGTTTAGGGTGGGCTACTCATTTGCGGGGCCTGCGCCCGCGAAGCTGTCCGACAAGGTTGTATCGACTAGCAGCTCGACAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCAC",
//	"TTGTATCGACTAGCAGCTCGACAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCACTTGTTGTG",
//	"GGTACGAGCCTCAGATACGCCGTCTTATGGGCTGGACAACTCCGTCCCGGCCATTTTCCG",
//	"ACTAGCAGCTCGACAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCACTTGTTGTGTGGGACCCC",
//	"ACTAGCAGCTCGACAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCACTTGTTGTGCGGGACCACCTTTACCTAAAGAGCTCGTGCGGTCTGACCCTAGCCGCGCGGAGTTTGGGATAATCCGATAGGAAGGACAAATTGAGCCTCACTAGACCAACCGCGAAAATAGAGAGGCAGGCAACCTCTCGAGGCGCTGACGGCAACTACCCCGCGTGCGA",
//	"GTACGAGCCTCAGATACGCCGTCTTATGGGCTGGACAACTCCGTCCCGGCCATTTTCCGGTCACGGCTTTAGCCGCCTGATTATGTGCGCCAAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGTCCAATTACGAAACGCGAGTCGATTTTTTATAGCGGTGCCCGT",
//	"GTACGAGCCTCAGATACGCCGTCTTATGGGCTGGACAACTCCGTCCCGGCCATTTTCCGTTCACGGCTTTAGCCGCCTGATTATGTGCGGTAAGCTCCTCAGGTGACAACAGCAGAAAAGAGTGGACCAATTACGAAACGCGAGTCGATTTTTTATAGCGGTGCCCGT",
//	"TCGACAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCACTTGTTGTGTGGGACCCCTTTTACCTAAAGAGCTCGTGCGGCCTAACCCTATCCGCGCGGGGTTCGGGATAATCCGATAGGAAGGACAAATTGAGCCTCGCTAGACCCACCACGAAAAAAGGGGGGTGGTCAACCTCTCGAGGCGCTGACGGCAACTACCCCACGTGCGA",
//	"TCGACAGCGCAACTGGTGAAGGCCACTTTCAGTCGTCGCCACTTGTTGTGTGGGACCCCCTTTACCTAAAGAGCTCGTGCGGCCTGACCCTAGCCGCGCGGAGTTCGGGATAATCCGATAGGAAG"
//};
//
