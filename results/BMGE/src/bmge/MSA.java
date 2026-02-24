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

import java.util.BitSet;

public abstract class MSA {

    private final int size;
    private final int length;
    private final String[] lbl;
    
    public MSA(final int size, final int length) {
	this.size = size;
	this.length = length;
	lbl = new String[size];
    }

    public final int size() {
	return size;
    }

    public final int length() {
	return length;
    }

    public final void setLabel(final int row, final String label) {
	lbl[row] = label;
    }

    public final String getLabel(final int row) {
	return lbl[row];
    }

    public abstract void setSequence(final int row, final String seq);
    public abstract String getCharState(final int row, final int col);
    public abstract String getMajorityCharState(final int col);
    public abstract double[] getFreq(final int col);
    public abstract double gapFreq(final int col);
    public abstract String getSequence(final int row);
    public abstract MSAaa toMSAaa();
    public abstract MSAco toMSAco();
    public abstract MSAnt toMSAnt();
    public abstract MSAry toMSAry();
    public abstract MSA trim(final BitSet mask);
}

