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
import java.util.Locale;

public class SMatrix {

    private final int size;
    private final double[][] mat;
    private final int[] id;

    // square matrix of size 0
    public SMatrix() {
	this(0);
    }
    
    // square matrix of specified size 
    public SMatrix(final int size) {
	this.size = size;
	mat = new double[size][size];
	id = IntStream.range(0, size).toArray();
    }
    
    // square matrix built from a specified 2d array
    public SMatrix(final double[][] a) {
	size = a.length;
	mat = new double[size][size];
	id = IntStream.range(0, size).toArray();
	for (int i: id) mat[i] = Arrays.copyOf(a[i], size);
    }
    
    // square diagonal matrix built from a specified array
    public SMatrix(final double[] v) {
	size = v.length;
	mat = new double[size][size];
	id = IntStream.range(0, size).toArray();
	for (int i: id) mat[i][i] = v[i];
    }

    // square matrix of specified size with a specified value
    public SMatrix(final int size, final double val) {
	this.size = size;
	mat = new double[size][size];
	id = IntStream.range(0, size).toArray();
	for (int i: id) Arrays.fill(mat[i], val);
    }

    // square identity matrix
    public static SMatrix identity(final int size) {
	final double[] v = new double[size];
	Arrays.fill(v, 1.0);
	return new SMatrix(v);
    }
    
    public final int size() {
	return size;
    }

    public final void set(final int row, final int col, final double val) {
	mat[row][col] = val;
    }
    
    public double get(final int row, final int col) {
	return mat[row][col];
    }

    // square matrix multiplication
    public final SMatrix mmult(final SMatrix m) {
	final SMatrix mp = new SMatrix(size);
	double mij;
	for (int i: id) 
	    for (int j: id) {
		mij = 0;
		for (int k: id) mij += mat[i][k] * m.get(k,j);
		mp.set(i, j, mij);
	    }
	return mp;
    }
	  
    public final SMatrix mult(final double s) {
	final SMatrix mp = new SMatrix(size);
	for (int i: id) 
	    for (int j: id) mp.set(i, j, s * mat[i][j]);
	return mp;
    }
	
    public final double trace() {
	double tr = 0; for (int i: id) tr += mat[i][i];
	return tr;
    }

    public final Eigen eigen() {
	return new Eigen(mat);
    }

    public final void print(final String format) {
	for (int i: id) {
	    for (int j: id) System.out.print(String.format(Locale.US, " %"+format+"f", mat[i][j]));
	    System.out.println("");
	}
    }
}
