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

public class MSAaa extends MSA {

    final static char[] ALPHABET = {'-','?','A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'};
    final static int  ALPHABET_SIZE = 25;
    final static byte GAP = (byte) 0;  // gap
    final static byte UNK = (byte) 1;  // unknown
    final static byte A   = (byte) 2;  // Alanine
    final static byte B   = (byte) 3;  // D or N
    final static byte C   = (byte) 4;  // Cysteine
    final static byte D   = (byte) 5;  // Aspartic acid
    final static byte E   = (byte) 6;  // Glutamic acid
    final static byte F   = (byte) 7;  // Phenylalanine 
    final static byte G   = (byte) 8;  // Glycine
    final static byte H   = (byte) 9;  // Histidine
    final static byte I   = (byte) 10; // Isoleucine
    final static byte K   = (byte) 11; // Lysine
    final static byte L   = (byte) 12; // Leucine
    final static byte M   = (byte) 13; // Methionine
    final static byte N   = (byte) 14; // Asparagine
    final static byte P   = (byte) 15; // Proline
    final static byte Q   = (byte) 16; // Glutamine
    final static byte R   = (byte) 17; //  Arginine
    final static byte S   = (byte) 18; // Serine
    final static byte T   = (byte) 19; // Threonine
    final static byte V   = (byte) 20; // Valine
    final static byte W   = (byte) 21; // Tryptophan
    final static byte X   = (byte) 22; // unknown character state
    final static byte Y   = (byte) 23; // Tyrosine
    final static byte Z   = (byte) 24; // Q or E

    final static char[] RESIDUES = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'}; //## WARNING: should be in the same order as the blosum matrices
    final static int  RESIDUE_NUMBER = 20;
    final static byte RA  = (byte) 0;  // Alanine
    final static byte RR  = (byte) 1;  // Arginine
    final static byte RN  = (byte) 2;  // Asparagine
    final static byte RD  = (byte) 3;  // Aspartic acid
    final static byte RC  = (byte) 4;  // Cysteine
    final static byte RQ  = (byte) 5;  // Glutamine
    final static byte RE  = (byte) 6;  // Glutamic acid
    final static byte RG  = (byte) 7;  // Glycine
    final static byte RH  = (byte) 8;  // Histidine
    final static byte RI  = (byte) 9;  // Isoleucine
    final static byte RL  = (byte) 10; // Leucine
    final static byte RK  = (byte) 11; // Lysine
    final static byte RM  = (byte) 12; // Methionine
    final static byte RF  = (byte) 13; // Phenylalanine 
    final static byte RP  = (byte) 14; // Proline
    final static byte RS  = (byte) 15; // Serine
    final static byte RT  = (byte) 16; // Threonine
    final static byte RW  = (byte) 17; // Tryptophan
    final static byte RY  = (byte) 18; // Tyrosine
    final static byte RV  = (byte) 19; // Valine

    private final byte[][] msa;
    private final int[] rid;

    public MSAaa(final int size, final int length) {
	super(size, length);
	msa = new byte[size][length];
	rid = IntStream.range(0, size).toArray();
    }

    private static final byte charCode(final char cs) {
	final byte c = (byte) Arrays.binarySearch(ALPHABET, Character.toUpperCase(cs));
	return (c < 0) ? UNK : c;
    }
    
    public final void setSequence(final int row, final String aaseq) {
	int p = aaseq.length();
	while ( --p >= 0 ) setCharCode(row, p, charCode(aaseq.charAt(p)));
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

    public final double[] getFreq(final int col) { //## WARNING: should be in the same order as the blosum matrices
	int r = RESIDUE_NUMBER;
	final double[] freq = new double[r];
	int total = size();
	for (int i: rid)
	    switch ( getCharCode(i, col) ) {
	    case A:  freq[RA]++;                         continue;
	    case B:  freq[RD] += 0.5;  freq[RN] += 0.5;  continue;
	    case C:  freq[RC]++;                         continue;
	    case D:  freq[RD]++;                         continue;
	    case E:  freq[RE]++;                         continue;
	    case F:  freq[RF]++;                         continue;
	    case G:  freq[RG]++;                         continue;
	    case H:  freq[RH]++;                         continue;
	    case I:  freq[RI]++;                         continue;
	    case K:  freq[RK]++;                         continue;
	    case L:  freq[RL]++;                         continue;
	    case M:  freq[RM]++;                         continue;
	    case N:  freq[RN]++;                         continue;
	    case P:  freq[RP]++;                         continue;
	    case Q:  freq[RQ]++;                         continue;
	    case R:  freq[RR]++;                         continue;
	    case S:  freq[RS]++;                         continue;
	    case T:  freq[RT]++;                         continue;
	    case V:  freq[RV]++;                         continue;
	    case W:  freq[RW]++;                         continue;
	    case X:  r = RESIDUE_NUMBER;
		     while ( --r >= 0 ) freq[r] += 0.05; continue;
	    case Y:  freq[RY]++;                         continue;
	    case Z:  freq[RE] += 0.5;  freq[RQ] += 0.5;  continue;
	    default: --total;                            continue;
	    }
	r = RESIDUE_NUMBER;
	if ( total > 0 ) while ( --r >= 0 ) freq[r] /= total;
	return freq;
    }

    public final String getMajorityCharState(final int col) {
	final int[] freq = new int[ALPHABET_SIZE];
	byte maj = UNK; int max = 0; 
	for (int i: rid) 
	    switch ( getCharCode(i, col) ) {
	    case A:  max = ( max < (++freq[A])) ? freq[maj=A] : max; continue;
	    case B:  max = ( max < (++freq[B])) ? freq[maj=B] : max; continue;
	    case C:  max = ( max < (++freq[C])) ? freq[maj=C] : max; continue;
	    case D:  max = ( max < (++freq[D])) ? freq[maj=D] : max; continue;
	    case E:  max = ( max < (++freq[E])) ? freq[maj=E] : max; continue;
	    case F:  max = ( max < (++freq[F])) ? freq[maj=F] : max; continue;
	    case G:  max = ( max < (++freq[G])) ? freq[maj=G] : max; continue;
	    case H:  max = ( max < (++freq[H])) ? freq[maj=H] : max; continue;
	    case I:  max = ( max < (++freq[I])) ? freq[maj=I] : max; continue;
	    case K:  max = ( max < (++freq[K])) ? freq[maj=K] : max; continue;
	    case L:  max = ( max < (++freq[L])) ? freq[maj=L] : max; continue;
	    case M:  max = ( max < (++freq[M])) ? freq[maj=M] : max; continue;
	    case N:  max = ( max < (++freq[N])) ? freq[maj=N] : max; continue;
	    case P:  max = ( max < (++freq[P])) ? freq[maj=P] : max; continue;
	    case Q:  max = ( max < (++freq[Q])) ? freq[maj=Q] : max; continue;
	    case R:  max = ( max < (++freq[R])) ? freq[maj=R] : max; continue;
	    case S:  max = ( max < (++freq[S])) ? freq[maj=S] : max; continue;
	    case T:  max = ( max < (++freq[T])) ? freq[maj=T] : max; continue;
	    case V:  max = ( max < (++freq[V])) ? freq[maj=V] : max; continue;
	    case W:  max = ( max < (++freq[W])) ? freq[maj=W] : max; continue;
	    case X:  max = ( max < (++freq[X])) ? freq[maj=X] : max; continue;
	    case Y:  max = ( max < (++freq[Y])) ? freq[maj=Y] : max; continue;
	    case Z:  max = ( max < (++freq[Z])) ? freq[maj=Z] : max; continue;
	    }
	return String.valueOf(ALPHABET[maj]);
    }

    public final double gapFreq(final int col) {
	double g = 0; for (int i: rid) g += ( getCharCode(i, col) == GAP ) ? 1 : 0;
	return g / size();
    }
	
    public final MSA trim(final BitSet mask) {
	final MSAaa tmsa = new MSAaa(size(), mask.cardinality());
	int p;
	for (int i: rid) {
	    tmsa.setLabel(i, getLabel(i));
	    //b = p = -1; while ( (b=mask.nextSetBit(++b)) != -1 ) tmsa.setCharCode(i, ++p, getCharCode(i, b));
	    p = -1; for (int b: mask.stream().toArray()) tmsa.setCharCode(i, ++p, getCharCode(i, b));
	}
	return tmsa;
    }

    public final MSAaa toMSAaa() {
	final MSAaa aa = new MSAaa(size(), length());
	int p;
	for (int i: rid) {
	    aa.setLabel(i, getLabel(i));
	    p = -1; for (byte c: msa[i]) aa.setCharCode(i, ++p, c);
	}
	return aa;
    }

    public final MSAco toMSAco() {
	final MSAco co = new MSAco(size(), length());
	int p;
	for (int i: rid) {
	    co.setLabel(i, getLabel(i));
	    p = length();
	    while ( --p >= 0 ) 
		switch ( getCharCode(i, p) ) {
		case GAP: co.setCharCode(i, p, MSAco.GAP); continue;
		case UNK: co.setCharCode(i, p, MSAco.UNK); continue;
		case A:   co.setCharCode(i, p, MSAco.GCX); continue;
		case R:   co.setCharCode(i, p, MSAco.MGX); continue;
		case N:   co.setCharCode(i, p, MSAco.AAY); continue;
		case D:   co.setCharCode(i, p, MSAco.GAY); continue;
		case C:   co.setCharCode(i, p, MSAco.TGY); continue;
		case Q:   co.setCharCode(i, p, MSAco.CAR); continue;
		case E:   co.setCharCode(i, p, MSAco.GAR); continue;
		case G:   co.setCharCode(i, p, MSAco.GGX); continue;
		case H:   co.setCharCode(i, p, MSAco.CAY); continue;
		case I:   co.setCharCode(i, p, MSAco.ATH); continue;
		case L:   co.setCharCode(i, p, MSAco.YTX); continue;
		case K:   co.setCharCode(i, p, MSAco.AAR); continue;
		case M:   co.setCharCode(i, p, MSAco.ATG); continue;
		case F:   co.setCharCode(i, p, MSAco.TTY); continue;
		case P:   co.setCharCode(i, p, MSAco.CCX); continue;
		case S:   co.setCharCode(i, p, MSAco.WSX); continue;
		case T:   co.setCharCode(i, p, MSAco.ACX); continue;
		case W:   co.setCharCode(i, p, MSAco.TGG); continue;
		case Y:   co.setCharCode(i, p, MSAco.TAY); continue;
		case V:   co.setCharCode(i, p, MSAco.GTX); continue;
		case X:   co.setCharCode(i, p, MSAco.XXX); continue;
		default:  co.setCharCode(i, p, MSAco.UNK); continue;
		}
	}
	return co;
    }	    
	    
    public final MSAnt toMSAnt() {
	return toMSAco().toMSAnt();
    }

    public final MSAry toMSAry() {
	final MSAry ry = new MSAry(size(), 3*length());
	int b, p;
	for (int i: rid) {
	    ry.setLabel(i, getLabel(i));
	    b = p = -1;
	    while ( ++p < length() ) 
		switch ( getCharCode(i, p) ) {
		case GAP: ry.setCharCode(i, ++b, MSAry.GAP); ry.setCharCode(i, ++b, MSAry.GAP); ry.setCharCode(i, ++b, MSAry.GAP); continue; // ---
		case UNK: ry.setCharCode(i, ++b, MSAry.UNK); ry.setCharCode(i, ++b, MSAry.UNK); ry.setCharCode(i, ++b, MSAry.UNK); continue; // ???
		case A:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.X);   continue; // RYX
		case R:   ry.setCharCode(i, ++b, MSAry.X);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.X);   continue; // XRX
		case N:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   continue; // RRY
		case D:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   continue; // RRY
		case C:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   continue; // YRY
		case Q:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   continue; // YRR
		case E:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   continue; // RRR
		case G:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.X);   continue; // RRX
		case H:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   continue; // YRY
		case I:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.X);   continue; // RYX
		case L:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.X);   continue; // YYX
		case K:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   continue; // RRR
		case M:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.R);   continue; // RYR
		case F:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.Y);   continue; // YYY
		case P:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.X);   continue; // YYX
		case S:   ry.setCharCode(i, ++b, MSAry.X);   ry.setCharCode(i, ++b, MSAry.X);   ry.setCharCode(i, ++b, MSAry.X);   continue; // XXX
		case T:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.X);   continue; // RYX
		case W:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.R);   continue; // YRR
		case Y:   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   continue; // YRY
		case V:   ry.setCharCode(i, ++b, MSAry.R);   ry.setCharCode(i, ++b, MSAry.Y);   ry.setCharCode(i, ++b, MSAry.X);   continue; // RYX
		case X:   ry.setCharCode(i, ++b, MSAry.X);   ry.setCharCode(i, ++b, MSAry.X);   ry.setCharCode(i, ++b, MSAry.X);   continue; // XXX
		default:  ry.setCharCode(i, ++b, MSAry.UNK); ry.setCharCode(i, ++b, MSAry.UNK); ry.setCharCode(i, ++b, MSAry.UNK); continue; // ???
		}
	}
	return ry;
    }	    
	    

}

