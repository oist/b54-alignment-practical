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

public class MSAnt extends MSA {

    final static char[] ALPHABET = {'-','?','A','B','C','D','G','H','K','M','N','R','S','T','U','V','W','X','Y'};
    final static int ALPHABET_SIZE = 19;
    final static byte GAP = (byte) 0;  // gap
    final static byte UNK = (byte) 1;  // unknown
    final static byte A   = (byte) 2;  // Adenosine
    final static byte B   = (byte) 3;  // C or G or T 
    final static byte C   = (byte) 4;  // Cytosine
    final static byte D   = (byte) 5;  // A or G or T
    final static byte G   = (byte) 6;  // Guanine
    final static byte H   = (byte) 7;  // A or C or T
    final static byte K   = (byte) 8;  // G or T
    final static byte M   = (byte) 9;  // A or C
    final static byte N   = (byte) 10; // unknown character state
    final static byte R   = (byte) 11; // A or G
    final static byte S   = (byte) 12; // C or G
    final static byte T   = (byte) 13; // Thymine / U  Uracil
    final static byte U   = (byte) 14; // Uracil
    final static byte V   = (byte) 15; // A or C or G
    final static byte W   = (byte) 16; // A or T 
    final static byte X   = (byte) 17; // unknown character state
    final static byte Y   = (byte) 18; // C or T

    final static char[] RESIDUES = {'A','C','G','T'};
    final static int RESIDUE_NUMBER = 4;
    final static byte RA  = (byte) 0;  // Adenosine
    final static byte RC  = (byte) 1;  // Cytosine
    final static byte RG  = (byte) 2;  // Guanine
    final static byte RT  = (byte) 3;  // Thymine / U Uracil

    private final static double THD = 1.0 / 3.0;
    private final static double FTH = 0.25;

    private final byte[][] msa;
    private final int[] rid;

    public MSAnt(final int size, final int length) {
	super(size, length);
	msa = new byte[size][length];
	rid = IntStream.range(0, size).toArray();
    }

    private final static byte charCode(final char cs) {
	final byte c = (byte) Arrays.binarySearch(ALPHABET, Character.toUpperCase(cs));
	return (c < 0) ? UNK : c;
    }
    
    public final void setSequence(final int row, final String seq) {
	int p = seq.length();
	while ( --p >= 0 ) setCharCode(row, p, charCode(seq.charAt(p)));
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
	    case A:  freq[RA]++;                                                         continue;
	    case B:                   freq[RC] += THD; freq[RG] += THD; freq[RT] += THD; continue;
	    case C:                   freq[RC]++;                                        continue;
	    case D:  freq[RA] += THD;                  freq[RG] += THD; freq[RT] += THD; continue;
	    case G:                                    freq[RG]++;                       continue;
	    case H:  freq[RA] += THD; freq[RC] += THD;                  freq[RT] += THD; continue;
	    case K:                                    freq[RG] += 0.5; freq[RT] += 0.5; continue;
	    case M:  freq[RA] += 0.5; freq[RC] += 0.5;                                   continue;
	    case N:  freq[RA] += FTH; freq[RC] += FTH; freq[RG] += FTH; freq[RT] += FTH; continue;
	    case R:  freq[RA] += 0.5;                  freq[RG] += 0.5;                  continue;
	    case S:                   freq[RC] += 0.5; freq[RG] += 0.5;                  continue;
	    case T:                                                     freq[RT]++;      continue;
	    case U:                                                     freq[RT]++;      continue;
	    case V:  freq[RA] += THD; freq[RC] += THD; freq[RG] += THD;                  continue;
	    case W:  freq[RA] += 0.5;                                   freq[RT] += 0.5; continue;
	    case X:  freq[RA] += FTH; freq[RC] += FTH; freq[RG] += FTH; freq[RT] += FTH; continue;
	    case Y:                   freq[RC] += 0.5;                  freq[RT] += 0.5; continue;
	    default: --total;                                                            continue;
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
	    case A:  max = ( max < (++freq[A])) ? freq[maj=A] : max; continue;
	    case B:  max = ( max < (++freq[B])) ? freq[maj=B] : max; continue;
	    case C:  max = ( max < (++freq[C])) ? freq[maj=C] : max; continue;
	    case D:  max = ( max < (++freq[D])) ? freq[maj=D] : max; continue;
	    case G:  max = ( max < (++freq[G])) ? freq[maj=G] : max; continue;
	    case H:  max = ( max < (++freq[H])) ? freq[maj=H] : max; continue;
	    case K:  max = ( max < (++freq[K])) ? freq[maj=K] : max; continue;
	    case M:  max = ( max < (++freq[M])) ? freq[maj=M] : max; continue;
	    case N:  max = ( max < (++freq[N])) ? freq[maj=N] : max; continue;
	    case R:  max = ( max < (++freq[R])) ? freq[maj=R] : max; continue;
	    case S:  max = ( max < (++freq[S])) ? freq[maj=S] : max; continue;
	    case T:  max = ( max < (++freq[T])) ? freq[maj=T] : max; continue;
	    case U:  max = ( max < (++freq[T])) ? freq[maj=T] : max; continue;
	    case V:  max = ( max < (++freq[V])) ? freq[maj=V] : max; continue;
	    case W:  max = ( max < (++freq[W])) ? freq[maj=W] : max; continue;
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
	final MSAnt tmsa = new MSAnt(size(), mask.cardinality());
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
	for (int i: rid) {
	    nt.setLabel(i, getLabel(i));
	    int p = -1; for (byte c: msa[i]) nt.setCharCode(i, ++p, c);
	}
	return nt;
    }
	
    public final MSAry toMSAry() {
	final MSAry ry = new MSAry(size(), length());
	int p;
	for (int i: rid) {
	    ry.setLabel(i, getLabel(i));
	    p = length();
	    while ( --p >= 0 ) 
		switch ( getCharCode(i, p) ) {
		case GAP: ry.setCharCode(i, p, MSAry.GAP); continue;
		case UNK: ry.setCharCode(i, p, MSAry.UNK); continue;
		case A:   ry.setCharCode(i, p, MSAry.R);   continue;
		case C:   ry.setCharCode(i, p, MSAry.Y);   continue;
		case G:   ry.setCharCode(i, p, MSAry.R);   continue;
		case N:   ry.setCharCode(i, p, MSAry.N);   continue;
		case R:   ry.setCharCode(i, p, MSAry.R);   continue;
		case T:   ry.setCharCode(i, p, MSAry.Y);   continue;
		case U:   ry.setCharCode(i, p, MSAry.Y);   continue;
		case X:   ry.setCharCode(i, p, MSAry.X);   continue;
		case Y:   ry.setCharCode(i, p, MSAry.Y);   continue;
		default:  ry.setCharCode(i, p, MSAry.X);   continue;
		}
	}
	return ry;
    }
	
}

