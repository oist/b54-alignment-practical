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

public class Utils {

    private Utils() {}

    // returns the square of a specified double value
    public static final double square (final double x) {
	return x*x;
    }
    
    // returns true when the specified array contains only zero
    public static final boolean isZero (final double[] a) {
	for (double x: a) if ( x != 0 ) return false;
	return true;
    }

    // returns the max value of the specified array
    public static final double max (final double[] a) {
	double max = Double.NEGATIVE_INFINITY;
	for (double x: a) max = ( max < x ) ? x : max;
	return max;
    }

    // returns log_b(a)
    public static final double logb (final double a, final double b) {
	return Math.log(a) / Math.log(b);
    }
	
    // formula (1) in Criscuolo and Gribaldo (2010) with S = I
    public static final double entropy (final double[] freq) {
	double h = max(freq);
	if ( h == 0 || h == 1 ) return 0;
	h = 0; for (double f: freq) h -= ((f > 0) ? f * logb(f, freq.length) : 0);
	return h;
    }

    // formula (1) in Criscuolo and Gribaldo (2010) with specified S
    public static final double entropy (final double[] freq, final SMatrix simm) {
	double h = max(freq);
	if ( h == 0 || h == 1 ) return 0;
	SMatrix mat = simm.mmult(new SMatrix(freq));
	mat = mat.mult(1.0 / mat.trace());
	double[] eig = mat.eigen().getRealEigenvalues();
	h = 0; for (double l: eig) h -= ((l > 0) ? l * logb(l, eig.length) : 0);
	return h;
    }
}
