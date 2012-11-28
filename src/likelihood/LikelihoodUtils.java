package likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.StringTokenizer;

import javax.print.attribute.IntegerSyntax;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.CharSetUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathUtils;

import com.google.common.base.CharMatcher;
import com.google.common.base.Charsets;
import com.google.common.base.Strings;
import com.google.common.collect.Iterables;
import com.google.common.primitives.Booleans;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.Chars;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Primitives;

public class LikelihoodUtils {

	private static final int ZERO = 0;
	private static final int ONE = 1;


	public static int hamDist(String s1, String s2)  {
//
//		if (s1.length() != s2.length()) {
//			System.err.println("Different length\t"+  s1.length() +"\t"+ s2.length()) ;
//			System.exit(-1);
//		}
		 
		int count = 0;
		int a;
		for (int i = 0; i < s1.length(); i++) {

			a = (s1.charAt(i) - s2.charAt(i));
			count+= (a==ZERO)? 0:1;
		}
		return count;
	}



	public static int hamDist(char[] c1, char[] c2) {


		int count = 0;

		for (int i = 0; i < c1.length; i++) {
//			System.out.println(c1[i] +"\t"+ c2[i] +"\t"+ (c1[i]-c2[i]));
			if ( (c1[i]-c2[i])!=0 ){
//			if (s1.charAt(i) != s2.charAt(i) ){
				count++;
			}
		}
		return count;
	}





	static final Byte bb = new Byte((byte) 0);
	
	
	public static int hamDist(char[] srCharArray, char[] hapCharArray, int start) {
		
		int count = 0;
		int a;
		
		for (int i = 0; i < srCharArray.length; i++) {
			a = (srCharArray[i] - hapCharArray[i+start]);
			
			count += (a==ZERO)? 0:1;

		}

		return count;
	}

//static int cc = 500;
//static int c2 = 0;//
	static int counts[] = new int[1000];
	static char[] srCharArray;
	static char[] hapCharArray;
	static int a, noComb;
	static char cc;
	public static int[] hamDistAll(char[] srCharArray2, char[] hapCharArray2) {
		srCharArray = srCharArray2;
		hapCharArray = hapCharArray2;
		
		int noComb = hapCharArray.length - srCharArray.length + 1; 		
//				int noComb = hLength - srLength + 1;
//System.out.println(noComb);
//		Arrays.fill(counts, 0);
//		int x = noComb*100;
//		int counts[] = new int[550];
//		if (noComb>counts.length) {
//			Arrays.fill(counts, 1);
//			c2++;
//			return counts;
//		counts = new int[noComb];
//			Arrays.fill(counts, 0);
//		int[] counts = Arrays.copyOfRange(counts2, 0, noComb);
//		int[] counts = Arrays.copyOf(counts2, noComb);
//		Arrays.fill(counts, 0, noComb, 0);
//		int[] counts = ArrayUtils.clone(counts2);
//		Ints.
		
//			int[] x = ArrayUtils.EMPTY_INT_ARRAY;
//			Ints.ensureCapacity(array, minLength, padding)
//			System.out.println(x.length);
//			System.out.println(noComb);
//		}
//		else {
//			c2++;
//		int a;
//		char cc;
			for (int i = 0; i < noComb; i++) {
				counts[i]=0;
			}
//		}

		
		for (int j = 0; j < srCharArray.length; j++) {
			cc = srCharArray[j];
			for (int i = 0; i < noComb; i++) {
				int a = (cc - hapCharArray[j+i]);
				counts[i] += (a==ZERO)? ZERO:ONE;
			}
		}

		return counts;
	}
}
