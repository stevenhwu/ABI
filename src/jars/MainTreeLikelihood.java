package jars;

import java.io.OutputStream;
import java.io.PrintStream;

import likelihood.LikelihoodCalculation;
import core.DataImporter;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
public class MainTreeLikelihood {

	/**
	 * @param args
	 */

	public static void main(String[] args) {

		if (args.length == 0) {
			System.out
					.println("Usage: java -jar testTreeLikelihood.jar alignmentFile treeFile");
		}

		String trueAlignmentFile = args[0];
		String truePhylogenyFile = args[1];

		System.setErr(new PrintStream(new OutputStream() {
			@Override
			public void write(int b) {
			}
		}));

		DataImporter dataImporter = new DataImporter();
		SimpleAlignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);

		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);

		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment);

		System.out.println(li.getTreeLikelihood());



	}

}
