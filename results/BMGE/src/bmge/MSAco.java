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

public class MSAco extends MSA {

    final static String[] ALPHABET = {"---","???","AAA","AAC","AAG","AAR","AAT","AAU","AAY","ACA","ACB","ACC","ACD","ACG","ACH","ACK","ACM","ACR","ACS","ACT","ACU","ACV","ACW","ACX","ACY","AGA","AGC","AGG","AGR","AGT","AGU","AGY","ATA","ATC","ATG","ATH","ATT","AUA","AUC","AUG","AUH","AUU","CAA","CAC","CAG","CAR","CAT","CAU","CAY","CCA","CCB","CCC","CCD","CCG","CCH","CCK","CCM","CCR","CCS","CCT","CCU","CCV","CCW","CCX","CCY","CGA","CGB","CGC","CGD","CGG","CGH","CGK","CGM","CGR","CGS","CGT","CGU","CGV","CGW","CGX","CGY","CTA","CTB","CTC","CTD","CTG","CTH","CTK","CTM","CTR","CTS","CTT","CTV","CTW","CTX","CTY","CUA","CUB","CUC","CUD","CUG","CUH","CUK","CUM","CUR","CUS","CUU","CUV","CUW","CUX","CUY","GAA","GAC","GAG","GAR","GAT","GAU","GAY","GCA","GCB","GCC","GCD","GCG","GCH","GCK","GCM","GCR","GCS","GCT","GCU","GCV","GCW","GCX","GCY","GGA","GGB","GGC","GGD","GGG","GGH","GGK","GGM","GGR","GGS","GGT","GGU","GGV","GGW","GGX","GGY","GTA","GTB","GTC","GTD","GTG","GTH","GTK","GTM","GTR","GTS","GTT","GTV","GTW","GTX","GTY","GUA","GUB","GUC","GUD","GUG","GUH","GUK","GUM","GUR","GUS","GUU","GUV","GUW","GUX","GUY","MGX","TAC","TAT","TAY","TCA","TCB","TCC","TCD","TCG","TCH","TCK","TCM","TCR","TCS","TCT","TCV","TCW","TCX","TCY","TGC","TGG","TGT","TGY","TTA","TTC","TTG","TTR","TTT","TTY","UAC","UAU","UAY","UCA","UCB","UCC","UCD","UCG","UCH","UCK","UCM","UCR","UCS","UCU","UCV","UCW","UCX","UCY","UGC","UGG","UGU","UGY","UUA","UUC","UUG","UUR","UUU","UUY","WSB","WSD","WSH","WSK","WSM","WSR","WSS","WSV","WSW","WSX","WSY","XXX","YTB","YTD","YTH","YTK","YTM","YTR","YTS","YTV","YTW","YTX","YTY","YUB","YUD","YUH","YUK","YUM","YUR","YUS","YUV","YUW","YUX","YUY"};
    final static int ALPHABET_SIZE = 271;
    final static short GAP = (short) 0;   // -
    final static short UNK = (short) 1;   // ?
    final static short AAA = (short) 2;   // K
    final static short AAC = (short) 3;   // N
    final static short AAG = (short) 4;   // K
    final static short AAR = (short) 5;   // K
    final static short AAT = (short) 6;   // N
    final static short AAU = (short) 7;   // N
    final static short AAY = (short) 8;   // N
    final static short ACA = (short) 9;   // T
    final static short ACB = (short) 10;  // T
    final static short ACC = (short) 11;  // T
    final static short ACD = (short) 12;  // T
    final static short ACG = (short) 13;  // T
    final static short ACH = (short) 14;  // T
    final static short ACK = (short) 15;  // T
    final static short ACM = (short) 16;  // T
    final static short ACR = (short) 17;  // T
    final static short ACS = (short) 18;  // T
    final static short ACT = (short) 19;  // T
    final static short ACU = (short) 20;  // T
    final static short ACV = (short) 21;  // T
    final static short ACW = (short) 22;  // T
    final static short ACX = (short) 23;  // T
    final static short ACY = (short) 24;  // T
    final static short AGA = (short) 25;  // R
    final static short AGC = (short) 26;  // S
    final static short AGG = (short) 27;  // R
    final static short AGR = (short) 28;  // R
    final static short AGT = (short) 29;  // S
    final static short AGU = (short) 30;  // S
    final static short AGY = (short) 31;  // S
    final static short ATA = (short) 32;  // I
    final static short ATC = (short) 33;  // I
    final static short ATG = (short) 34;  // M
    final static short ATH = (short) 35;  // I
    final static short ATT = (short) 36;  // I
    final static short AUA = (short) 37;  // I
    final static short AUC = (short) 38;  // I
    final static short AUG = (short) 39;  // M
    final static short AUH = (short) 40;  // I
    final static short AUU = (short) 41;  // I
    final static short CAA = (short) 42;  // Q
    final static short CAC = (short) 43;  // H
    final static short CAG = (short) 44;  // Q
    final static short CAR = (short) 45;  // Q
    final static short CAT = (short) 46;  // H
    final static short CAU = (short) 47;  // H
    final static short CAY = (short) 48;  // H
    final static short CCA = (short) 49;  // P
    final static short CCB = (short) 50;  // P
    final static short CCC = (short) 51;  // P
    final static short CCD = (short) 52;  // P
    final static short CCG = (short) 53;  // P
    final static short CCH = (short) 54;  // P
    final static short CCK = (short) 55;  // P
    final static short CCM = (short) 56;  // P
    final static short CCR = (short) 57;  // P
    final static short CCS = (short) 58;  // P
    final static short CCT = (short) 59;  // P
    final static short CCU = (short) 60;  // P
    final static short CCV = (short) 61;  // P
    final static short CCW = (short) 62;  // P
    final static short CCX = (short) 63;  // P
    final static short CCY = (short) 64;  // P
    final static short CGA = (short) 65;  // R
    final static short CGB = (short) 66;  // R
    final static short CGC = (short) 67;  // R
    final static short CGD = (short) 68;  // R
    final static short CGG = (short) 69;  // R
    final static short CGH = (short) 70;  // R
    final static short CGK = (short) 71;  // R
    final static short CGM = (short) 72;  // R
    final static short CGR = (short) 73;  // R
    final static short CGS = (short) 74;  // R
    final static short CGT = (short) 75;  // R
    final static short CGU = (short) 76;  // R
    final static short CGV = (short) 77;  // R
    final static short CGW = (short) 78;  // R
    final static short CGX = (short) 79;  // R
    final static short CGY = (short) 80;  // R
    final static short CTA = (short) 81;  // L
    final static short CTB = (short) 82;  // L
    final static short CTC = (short) 83;  // L
    final static short CTD = (short) 84;  // L
    final static short CTG = (short) 85;  // L
    final static short CTH = (short) 86;  // L
    final static short CTK = (short) 87;  // L
    final static short CTM = (short) 88;  // L
    final static short CTR = (short) 89;  // L
    final static short CTS = (short) 90;  // L
    final static short CTT = (short) 91;  // L
    final static short CTV = (short) 92;  // L
    final static short CTW = (short) 93;  // L
    final static short CTX = (short) 94;  // L
    final static short CTY = (short) 95;  // L
    final static short CUA = (short) 96;  // L
    final static short CUB = (short) 97;  // L
    final static short CUC = (short) 98;  // L
    final static short CUD = (short) 99;  // L
    final static short CUG = (short) 100; // L
    final static short CUH = (short) 101; // L
    final static short CUK = (short) 102; // L
    final static short CUM = (short) 103; // L
    final static short CUR = (short) 104; // L
    final static short CUS = (short) 105; // L
    final static short CUU = (short) 106; // L
    final static short CUV = (short) 107; // L
    final static short CUW = (short) 108; // L
    final static short CUX = (short) 109; // L
    final static short CUY = (short) 110; // L
    final static short GAA = (short) 111; // E
    final static short GAC = (short) 112; // D
    final static short GAG = (short) 113; // E
    final static short GAR = (short) 114; // E
    final static short GAT = (short) 115; // D
    final static short GAU = (short) 116; // D
    final static short GAY = (short) 117; // D
    final static short GCA = (short) 118; // A
    final static short GCB = (short) 119; // A
    final static short GCC = (short) 120; // A
    final static short GCD = (short) 121; // A
    final static short GCG = (short) 122; // A
    final static short GCH = (short) 123; // A
    final static short GCK = (short) 124; // A
    final static short GCM = (short) 125; // A
    final static short GCR = (short) 126; // A
    final static short GCS = (short) 127; // A
    final static short GCT = (short) 128; // A
    final static short GCU = (short) 129; // A
    final static short GCV = (short) 130; // A
    final static short GCW = (short) 131; // A
    final static short GCX = (short) 132; // A
    final static short GCY = (short) 133; // A
    final static short GGA = (short) 134; // G
    final static short GGB = (short) 135; // G
    final static short GGC = (short) 136; // G
    final static short GGD = (short) 137; // G
    final static short GGG = (short) 138; // G
    final static short GGH = (short) 139; // G
    final static short GGK = (short) 140; // G
    final static short GGM = (short) 141; // G
    final static short GGR = (short) 142; // G
    final static short GGS = (short) 143; // G
    final static short GGT = (short) 144; // G
    final static short GGU = (short) 145; // G
    final static short GGV = (short) 146; // G
    final static short GGW = (short) 147; // G
    final static short GGX = (short) 148; // G
    final static short GGY = (short) 149; // G
    final static short GTA = (short) 150; // V
    final static short GTB = (short) 151; // V
    final static short GTC = (short) 152; // V
    final static short GTD = (short) 153; // V
    final static short GTG = (short) 154; // V
    final static short GTH = (short) 155; // V
    final static short GTK = (short) 156; // V
    final static short GTM = (short) 157; // V
    final static short GTR = (short) 158; // V
    final static short GTS = (short) 159; // V
    final static short GTT = (short) 160; // V
    final static short GTV = (short) 161; // V
    final static short GTW = (short) 162; // V
    final static short GTX = (short) 163; // V
    final static short GTY = (short) 164; // V
    final static short GUA = (short) 165; // V
    final static short GUB = (short) 166; // V
    final static short GUC = (short) 167; // V
    final static short GUD = (short) 168; // V
    final static short GUG = (short) 169; // V
    final static short GUH = (short) 170; // V
    final static short GUK = (short) 171; // V
    final static short GUM = (short) 172; // V
    final static short GUR = (short) 173; // V
    final static short GUS = (short) 174; // V
    final static short GUU = (short) 175; // V
    final static short GUV = (short) 176; // V
    final static short GUW = (short) 177; // V
    final static short GUX = (short) 178; // V
    final static short GUY = (short) 179; // V
    final static short MGX = (short) 180; // R
    final static short TAC = (short) 181; // Y
    final static short TAT = (short) 182; // Y
    final static short TAY = (short) 183; // Y
    final static short TCA = (short) 184; // S
    final static short TCB = (short) 185; // S
    final static short TCC = (short) 186; // S
    final static short TCD = (short) 187; // S
    final static short TCG = (short) 188; // S
    final static short TCH = (short) 189; // S
    final static short TCK = (short) 190; // S
    final static short TCM = (short) 191; // S
    final static short TCR = (short) 192; // S
    final static short TCS = (short) 193; // S
    final static short TCT = (short) 194; // S
    final static short TCV = (short) 195; // S
    final static short TCW = (short) 196; // S
    final static short TCX = (short) 197; // S
    final static short TCY = (short) 198; // S
    final static short TGC = (short) 199; // C
    final static short TGG = (short) 200; // W
    final static short TGT = (short) 201; // C
    final static short TGY = (short) 202; // C
    final static short TTA = (short) 203; // L
    final static short TTC = (short) 204; // F
    final static short TTG = (short) 205; // L
    final static short TTR = (short) 206; // L
    final static short TTT = (short) 207; // F
    final static short TTY = (short) 208; // F
    final static short UAC = (short) 209; // Y
    final static short UAU = (short) 210; // Y
    final static short UAY = (short) 211; // Y
    final static short UCA = (short) 212; // S
    final static short UCB = (short) 213; // S
    final static short UCC = (short) 214; // S
    final static short UCD = (short) 215; // S
    final static short UCG = (short) 216; // S
    final static short UCH = (short) 217; // S
    final static short UCK = (short) 218; // S
    final static short UCM = (short) 219; // S
    final static short UCR = (short) 220; // S
    final static short UCS = (short) 221; // S
    final static short UCU = (short) 222; // S
    final static short UCV = (short) 223; // S
    final static short UCW = (short) 224; // S
    final static short UCX = (short) 225; // S
    final static short UCY = (short) 226; // S
    final static short UGC = (short) 227; // C
    final static short UGG = (short) 228; // W
    final static short UGU = (short) 229; // C
    final static short UGY = (short) 230; // C
    final static short UUA = (short) 231; // L
    final static short UUC = (short) 232; // F
    final static short UUG = (short) 233; // L
    final static short UUR = (short) 234; // L
    final static short UUU = (short) 235; // F
    final static short UUY = (short) 236; // F
    final static short WSB = (short) 237; // S
    final static short WSD = (short) 238; // S
    final static short WSH = (short) 239; // S
    final static short WSK = (short) 240; // S
    final static short WSM = (short) 241; // S
    final static short WSR = (short) 242; // S
    final static short WSS = (short) 243; // S
    final static short WSV = (short) 244; // S
    final static short WSW = (short) 245; // S
    final static short WSX = (short) 246; // S
    final static short WSY = (short) 247; // S
    final static short XXX = (short) 248; // X
    final static short YTB = (short) 249; // L
    final static short YTD = (short) 250; // L
    final static short YTH = (short) 251; // L
    final static short YTK = (short) 252; // L
    final static short YTM = (short) 253; // L
    final static short YTR = (short) 254; // L
    final static short YTS = (short) 255; // L
    final static short YTV = (short) 256; // L
    final static short YTW = (short) 257; // L
    final static short YTX = (short) 258; // L
    final static short YTY = (short) 259; // L
    final static short YUB = (short) 260; // L
    final static short YUD = (short) 261; // L
    final static short YUH = (short) 262; // L
    final static short YUK = (short) 263; // L
    final static short YUM = (short) 264; // L
    final static short YUR = (short) 265; // L
    final static short YUS = (short) 266; // L
    final static short YUV = (short) 267; // L
    final static short YUW = (short) 268; // L
    final static short YUX = (short) 269; // L
    final static short YUY = (short) 270; // L
    private final static byte[] CO2AA = {MSAaa.GAP,MSAaa.UNK,MSAaa.K,MSAaa.N,MSAaa.K,MSAaa.K,MSAaa.N,MSAaa.N,MSAaa.N,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.T,MSAaa.R,MSAaa.S,MSAaa.R,MSAaa.R,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.I,MSAaa.I,MSAaa.M,MSAaa.I,MSAaa.I,MSAaa.I,MSAaa.I,MSAaa.M,MSAaa.I,MSAaa.I,MSAaa.Q,MSAaa.H,MSAaa.Q,MSAaa.Q,MSAaa.H,MSAaa.H,MSAaa.H,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.P,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.R,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.E,MSAaa.D,MSAaa.E,MSAaa.E,MSAaa.D,MSAaa.D,MSAaa.D,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.A,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.G,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.V,MSAaa.R,MSAaa.Y,MSAaa.Y,MSAaa.Y,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.C,MSAaa.W,MSAaa.C,MSAaa.C,MSAaa.L,MSAaa.F,MSAaa.L,MSAaa.L,MSAaa.F,MSAaa.F,MSAaa.Y,MSAaa.Y,MSAaa.Y,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.C,MSAaa.W,MSAaa.C,MSAaa.C,MSAaa.L,MSAaa.F,MSAaa.L,MSAaa.L,MSAaa.F,MSAaa.F,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.S,MSAaa.X,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L,MSAaa.L};

    private final short[][] msa;
    private final byte[][] msaa;
    private final int[] rid;

    public MSAco(final int size, final int length) { //## WARNING: length = no. codons
	super(size, length);
	msa = new short[size][length];
	msaa = new byte[size][length];
	rid = IntStream.range(0, size).toArray();
    }

    private final static short charCode(final String cod) {
	final short c = (short) Arrays.binarySearch(ALPHABET, cod.toUpperCase());
	return (c < 0) ? UNK : c;
    }

    public final void setSequence(final int row, final String coseq) {
	int p = coseq.length() / 3;
	while ( --p >= 0 ) setCharCode(row, p, charCode(coseq.substring(3*p,3*(p+1))));
    }

    public final void setCharCode(final int row, final int col, final short code) {
	msa[row][col] = code;
	msaa[row][col] = CO2AA[code];
    }

    public final short getCharCode(final int row, final int col) {
	return msa[row][col];
    }

    public final String getCharState(final int row, final int col) {
	return ALPHABET[getCharCode(row,col)];
    }

    public final String getSequence(final int row) {
	StringBuilder sb = new StringBuilder(3*length());
	for (short c: msa[row]) sb = sb.append(ALPHABET[c]);
	return sb.toString();
    }

    public final double[] getFreq(final int col) { //## WARNING: returns aa freq
	int r = MSAaa.RESIDUE_NUMBER;
	final double[] freq = new double[r];
	int total = size();
	for (int i: rid) 
	    switch ( msaa[i][col] ) {
	    case MSAaa.A:  freq[MSAaa.RA]++;                               continue;
	    case MSAaa.B:  freq[MSAaa.RD] += 0.5;  freq[MSAaa.RN] += 0.5;  continue;
	    case MSAaa.C:  freq[MSAaa.RC]++;                               continue;
	    case MSAaa.D:  freq[MSAaa.RD]++;                               continue;
	    case MSAaa.E:  freq[MSAaa.RE]++;                               continue;
	    case MSAaa.F:  freq[MSAaa.RF]++;                               continue;
	    case MSAaa.G:  freq[MSAaa.RG]++;                               continue;
	    case MSAaa.H:  freq[MSAaa.RH]++;                               continue;
	    case MSAaa.I:  freq[MSAaa.RI]++;                               continue;
	    case MSAaa.K:  freq[MSAaa.RK]++;                               continue;
	    case MSAaa.L:  freq[MSAaa.RL]++;                               continue;
	    case MSAaa.M:  freq[MSAaa.RM]++;                               continue;
	    case MSAaa.N:  freq[MSAaa.RN]++;                               continue;
	    case MSAaa.P:  freq[MSAaa.RP]++;                               continue;
	    case MSAaa.Q:  freq[MSAaa.RQ]++;                               continue;
	    case MSAaa.R:  freq[MSAaa.RR]++;                               continue;
	    case MSAaa.S:  freq[MSAaa.RS]++;                               continue;
	    case MSAaa.T:  freq[MSAaa.RT]++;                               continue;
	    case MSAaa.V:  freq[MSAaa.RV]++;                               continue;
	    case MSAaa.W:  freq[MSAaa.RW]++;                               continue;
	    case MSAaa.X:  r = MSAaa.RESIDUE_NUMBER;
		     while ( --r >= 0 ) freq[r] += 0.05;                   continue;
	    case MSAaa.Y:  freq[MSAaa.RY]++;                               continue;
	    case MSAaa.Z:  freq[MSAaa.RE] += 0.5;  freq[MSAaa.RQ] += 0.5;  continue;
	    default: --total;                                              continue;
	    }
	r = MSAaa.RESIDUE_NUMBER;
	if ( total > 0 ) while ( --r >= 0 ) freq[r] /= total;
	return freq;
    }
	
    public final String getMajorityCharState(final int col) {
	final int[] freq = new int[ALPHABET_SIZE];
	short r, maj = UNK;
	for (int i: rid) if ( (r=getCharCode(i, col)) > UNK ) ++freq[r];
	r = ALPHABET_SIZE; int max = 0; 
	while ( --r >= 0 ) if ( max < freq[r] ) max = freq[maj=r];
	return ALPHABET[maj];
    }

    public final double gapFreq(final int col) {
	double g = 0; for (int i: rid) g += ( getCharCode(i, col) == GAP ) ? 1 : 0;
	return g / size();
    }
	
    public final MSA trim(final BitSet mask) {
	final MSAco tmsa = new MSAco(size(), mask.cardinality());
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
	    p = length(); while ( --p >= 0 ) aa.setCharCode(i, p, msaa[i][p]);
	}
	return aa;
    }

    public final MSAco toMSAco() {
	final MSAco co = new MSAco(size(), length());
	int p;
	for (int i: rid) {
	    co.setLabel(i, getLabel(i));
	    p = -1; for (short c: msa[i]) co.setCharCode(i, ++p, c);
	}
	return co;
    }
	
    public final MSAnt toMSAnt() {
	final MSAnt nt = new MSAnt(size(), 3*length());
	for (int i: rid) {
	    nt.setLabel(i, getLabel(i));
	    nt.setSequence(i, getSequence(i));
	}
	return nt;
    }

    public final MSAry toMSAry() {
	return toMSAnt().toMSAry();
    }
}

