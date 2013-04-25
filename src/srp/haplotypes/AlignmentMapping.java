package srp.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequence;
import dr.math.MathUtils;

public class AlignmentMapping {

	private static final int GAP = '-';
	
	private ArrayList<Integer>[] mapToSrp; // each [] = position, each ArrayList = map to which read
	private HashSet<Character>[] setsOfAvailableChar;

	private ArrayList<Character>[] listOfAvailableChar;
	
	private HashMap<String, Integer> seqNameToSeqID; // map sequence_name >xxx to int

	private ArrayList<ShortRead> shortReads;

	private int length;
	private Integer srpCount;

	
	
	public AlignmentMapping(Alignment srpAlignment){
		
		init( srpAlignment.getSiteCount() );
		for (int i = 0; i < srpAlignment.getSequenceCount(); i++) {
			Sequence s = srpAlignment.getSequence(i);
			addSequence(s);
		}
		for (int i = 0; i < length; i++) {
			listOfAvailableChar[i] = new ArrayList<Character>(setsOfAvailableChar[i]);
		}
		setsOfAvailableChar=null;
	}


	private void init(int l){
			length = l;
			mapToSrp = new ArrayList[length];
			setsOfAvailableChar = new HashSet[length];
			listOfAvailableChar = new ArrayList[length];
			
			for (int i = 0; i < this.length; i++) {
				mapToSrp[i] = new ArrayList<Integer>(); 
				setsOfAvailableChar[i] = new HashSet<Character>();
			}
	
			seqNameToSeqID = new HashMap<String, Integer>();
			shortReads = new ArrayList<ShortRead>();
	//		fullSrp = new HashMap<>();
			srpCount = 0;
			
			
		}


	private void addSequence(Sequence s) {

		ShortRead srp = new ShortRead(s);
		if (srp.getIsValid()){
			seqNameToSeqID.put(srp.getName(), srpCount);
			shortReads.add(srpCount, srp);
			
			for (int j = srp.getStart(); j < srp.getEnd(); j++) {
				mapToSrp[j].add(srpCount);
				setsOfAvailableChar[j].add(srp.getFullSrpCharAt(j));
			}
			srpCount++;
		}
		
	}

	@Override
	public String toString() {
		String s = "";
		for (int i = 0; i < mapToSrp.length; i++) {
			s += mapToSrp[i].toString();
			s += "\n";

		}
		return s;

	}

	public ArrayList<Integer> getMapToSrp(int pos) {
		return mapToSrp[pos];
	}
	public Integer mapNameToID(String name){
		return seqNameToSeqID.get(name);
	}

	public int getLength() {
		return length;
	}

	public int getSrpCount() {
		return srpCount;
	}

	public ShortRead getShortRead(int i) {
		return shortReads.get(i);
	}

	public char getShortReadCharAt(int index, int c) {
		return shortReads.get(index).getFullSrpCharAt(c);
	}

	public String getSrpFull(int i) {

		return shortReads.get(i).getFullSrp();
	}

	public String getSrpFragment(int i) {

		return shortReads.get(i).getFragmentSrp();
	}

	public int getSrpStart(int i) {

		return shortReads.get(i).getStart();
	}

	public int getSrpEnd(int i) {

		return shortReads.get(i).getEnd();
	}

	public int getSrpLength(int i) {
		return shortReads.get(i).getLength();
	}

	public String getSrpName(int i) {
		return shortReads.get(i).getName();
	}


	public int[] nextBase() {

		int pos = MathUtils.nextInt(length);
		int newChar = nextBaseAt(pos);
	
		return new int[]{pos, newChar};
	}
	public int nextBaseAt(int pos){
		int newChar = GAP;
		int size = mapToSrp[pos].size();
		if (size != 0) {
			int srpIndex = mapToSrp[pos].get(MathUtils.nextInt(size));
			newChar = getShortReadCharAt(srpIndex, pos);
		}
		
		return newChar;
	}


	public int[] nextBaseUniform() {
		int newChar = GAP;
		int pos = MathUtils.nextInt(length);
		int size = listOfAvailableChar[pos].size();
		if (size != 0) {
			newChar = listOfAvailableChar[pos].get(MathUtils.nextInt(size));
		}
		
		
		return new int[]{pos, newChar};
	}
}
