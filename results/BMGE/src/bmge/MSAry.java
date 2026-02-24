/*
  ########################################################################################################

  BMGE: Block Mapping and Gathering with Entropy
    
  Copyright (C) 2010-2024  Institut Pasteur
  
  This file is part of BMGE.

  BMGE is  free software:  you can  redistribute it  and/or modify it under  the terms of the  GNU General 
  Public License as  published by the  Free Software Foundation,  either version 3 of the License,  or (at 
  your option) any later version.

  BMGE is  distributed in  the hope  that it will  be useful,  but WITHOUT ANY WARRANTY;  without even the 
  implied warranty  of MERCHANTABILITY or  FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU  General Public 
  License for more details.

  You should have received a copy of the  GNU General Public License along with this program.  If not, see
  <http://www.gnu.org/licenses/>.
  
  Contact:
   Alexis Criscuolo                                                            alexis.criscuolo@pasteur.fr
   Genome Informatics & Phylogenetics (GIPhy)                                             giphy.pasteur.fr
   Bioinformatics and Biostatistics Hub                                 research.pasteur.fr/team/hub-giphy
   Institut Pasteur, Paris, FRANCE           research.pasteur.fr/team/bioinformatics-and-biostatistics-hub

   Simonetta Gribaldo                                                        simonetta.gribaldo@pasteur.fr
   Evolutionary Biology of the Microbial Cell                         https://research.pasteur.fr/en/b/6jQ
   Institut Pasteur, Paris, FRANCE                                             https://research.pasteur.fr

  ########################################################################################################
*/

package bmge;

import java.util.stream.IntStream;
import java.util.Arrays;
import java.util.BitSet;

public class MSAry extends MSA {

    final static char[] ALPHABET = {'-','?','N','R','X','Y'};
    final static int ALPHABET_SIZE = 6;
    final static byte GAP = (byte) 0;  // gap
    final static byte UNK = (byte) 1;  // unknown
    final static byte N   = (byte) 2;  // unknown character state
    final static byte R   = (byte) 3;  // A or G
    final static byte X   = (byte) 4;  // unknown character state
    final static byte Y   = (byte) 5;  // C or T

    final static char[] RESIDUES = {'R','Y'};
    final static int RESIDUE_NUMBER = 2;
    final static byte RR  = (byte) 0;  // A or G
    final static byte RY  = (byte) 1;  // C or T

    final private byte[][] msa;
    final private int[] rid;

    public MSAry(final int size, final int length) {
	super(size, length);
	msa = new byte[size][length];
	rid = IntStream.range(0, size).toArray();
    }

    private final static byte charCode(final char cs) {
	final byte c = (byte) Arrays.binarySearch(ALPHABET, Character.toUpperCase(cs));
	return (c < 0) ? UNK : c;
    }
    
    public final void setSequence(final int row, final String ryseq) {
	int p = ryseq.length();
	while ( --p >= 0 ) setCharCode(row, p, charCode(ryseq.charAt(p)));
    }

    public final void setCharCode(final int row, final int col, final byte code) {
	msa[row][col] = code;
    }

    public final byte getCharCode(final int row, final int col) {
	return msa[row][col];
    }

    public final String getCharState(final int row, final int col) {
	return String.valueOf(ALPHABET[getCharCode(row,col)]);
    }

    public final String getSequence(final int row) {
	final StringBuilder sb = new StringBuilder("x".repeat(length()));
	int p = -1; for (byte c: msa[row]) sb.setCharAt(++p, ALPHABET[c]);
	return sb.toString();
    }

    public final double[] getFreq(final int col) {
	int r = RESIDUE_NUMBER;
	final double[] freq = new double[r];
	int total = size();
	for (int i: rid)
	    switch ( getCharCode(i, col) ) {
	    case N:  freq[RR] += 0.5; freq[RY] += 0.5; continue;
	    case R:  freq[RR]++;                       continue;
	    case X:  freq[RR] += 0.5; freq[RY] += 0.5; continue;
	    case Y:                   freq[RY] += 0.5; continue;
	    default: --total;                          continue;
	    }
	r = RESIDUE_NUMBER;
	if ( total > 0 ) while ( --r >= 0 ) freq[r] /= total;
	return freq;
    }
	
    public final String getMajorityCharState(final int col) {
	byte maj = UNK;
	int max = ALPHABET_SIZE; 
	final int[] freq = new int[max];
	for (int i: rid) 
	    switch ( getCharCode(i, col) ) {
	    case N:  max = ( max < (++freq[N])) ? freq[maj=N] : max; continue;
	    case R:  max = ( max < (++freq[R])) ? freq[maj=R] : max; continue;
	    case X:  max = ( max < (++freq[X])) ? freq[maj=X] : max; continue;
	    case Y:  max = ( max < (++freq[Y])) ? freq[maj=Y] : max; continue;
	    }
	return String.valueOf(ALPHABET[maj]);
    }

    public final double gapFreq(final int col) {
	double g = 0; for (int i: rid) g += ( getCharCode(i, col) == GAP ) ? 1 : 0;
	return g / size();
    }
	
    public final MSA trim(final BitSet mask) {
	final MSAry tmsa = new MSAry(size(), mask.cardinality());
	int p;
	for (int i: rid) {
	    tmsa.setLabel(i, getLabel(i));
	    //b = p = -1; while ( (b=mask.nextSetBit(++b)) != -1 ) tmsa.setCharCode(i, ++p, getCharCode(i, b));
	    p = -1; for (int b: mask.stream().toArray()) tmsa.setCharCode(i, ++p, getCharCode(i, b));
	}
	return tmsa;
    }

    public final MSAaa toMSAaa() {
	return toMSAco().toMSAaa();
    }

    public final MSAco toMSAco() {
	int l = length() / 3;
	final MSAco co = new MSAco(size(), l);
	for (int i: rid) {
	    co.setLabel(i, getLabel(i));
	    co.setSequence(i, getSequence(i).substring(0, 3*l));
	}
	return co;
    }
    
    public final MSAnt toMSAnt() {
	final MSAnt nt = new MSAnt(size(), length());
	int p;
	for (int i: rid) {
	    nt.setLabel(i, getLabel(i));
	    p = length(); while ( --p >= 0 ) nt.setCharCode(i, p, getCharCode(i, p));
	}
	return nt;
    }

    public final MSAry toMSAry() {
	final MSAry ry = new MSAry(size(), length());
	for (int i: rid) {
	    ry.setLabel(i, getLabel(i));
	    int p = -1; for (byte c: msa[i]) ry.setCharCode(i, ++p, c);
	}
	return ry;
    }
	
}

