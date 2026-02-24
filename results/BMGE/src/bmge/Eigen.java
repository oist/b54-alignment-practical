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

import java.util.Arrays;
import java.util.stream.IntStream;

/*
  NOTE:
  This class is a complete rewriting of the class EigenvalueDecomposition from the package JAMA.
   + JAMA:                        https://math.nist.gov/javanumerics/jama/
   + EigenvalueDecomposition:     https://math.nist.gov/javanumerics/jama/doc/Jama/EigenvalueDecomposition.html
   + EigenvalueDecomposition src: https://github.com/fiji/Jama/blob/master/src/main/java/Jama/EigenvalueDecomposition.java
*/

public class Eigen {
    
    //## constants
    private static final double EPS = Math.pow(2.0,-52.0);
    
    //## data
    private static int n;                //## NOTE: size
    private static int n_1;              //## NOTE: size-1
    private final  int[] id0n;           //## NOTE: int array from 0 to n-1
    private final  int[] id0n_1;         //## NOTE: int array from 0 to n-2
    private final  int[] id1n;           //## NOTE: int array from 1 to n-1
    private static boolean issymmetric;  //## NOTE: true if symmetric
    private final  double[] d, e;        //## NOTE: eigenvalues
    private final  double[][] V;         //## NOTE: eigenvectors
    /*private final  double[][] inV;       //## NOTE: inverse eigenvectors*/

    //## stuffs
    private static int i, i_1, i1, j, k, l, m, m1, m_1, nn, iter;
    private static double c, c2, c3, f, g, h, p, q, r, s, s2, t, w, x, y, z;
    private static double dl1, el1, exshift, hh, norm, scale, tst1, ra, sa, vr, vi;
    private static boolean notlast;
    private static int[] id0i, id0i1, idmhigh1;
    private static double[][] H;
    private static double[] ort;
    
    public Eigen (final double[][] mat) {
	n_1 = n = mat.length;
	--n_1;
	
	id0n   = IntStream.range(0, n).toArray();
	id0n_1 = IntStream.range(0, n_1).toArray();
	id1n   = IntStream.range(1, n).toArray();
	
	V   = new double[n][n];
	d   = new double[n];
	e   = new double[n];
	
	issymmetric = true;
	i = n;
	while ( issymmetric && (j=--i) >= 0 )
	    while ( issymmetric && --j >= 0 )
		issymmetric = ( mat[i][j] == mat[j][i] );
	
	if ( issymmetric ) {
	    for (int i: id0n) System.arraycopy(mat[i], 0, V[i], 0, n);
	    tred2();  //## tridiagonalizing
	    tql2();   //## diagonalizing
	}
	else {
	    H = new double[n][n];
	    ort = new double[n];
	    for (int i: id0n) System.arraycopy(mat[i], 0, H[i], 0, n);
	    orthes(); //## reducing to Hessenberg form.
	    hqr2();   //## reducing Hessenberg to real Schur form.
	}
	/*
	  inV = new double[n][n];
	  inv();        //## inverting V
	*/
    }

    //## returns the eigenvector matrix
    public final double[][] getEigenvectors ()    { return Arrays.copyOf(V, n);   }
    //## returns the real parts of the eigenvalues
    public final double[] getRealEigenvalues ()   { return Arrays.copyOf(d, n);   }
    //## returns the imaginary parts of the eigenvalues
    public final double[] getImagEigenvalues ()   { return Arrays.copyOf(e, n);   }
    //## returns the inverse eigenvector matrix
    /*public final double[][] getInvEigenvectors () { return Arrays.copyOf(inV, n); }*/
        
    //## returns x*x
    private static final double square(final double x) { return x*x; }
    
    //## returns sqrt(a^2 + b^2) without under/overflow.
    private static final double hypot(final double a, final double b) { return ( Math.abs(a) > Math.abs(b) ) ? Math.abs(a) * Math.sqrt(1+square(b/a)) : ( b != 0 ) ? Math.abs(b) * Math.sqrt(1+square(a/b)) : 0; }
    
    //## Symmetric Householder reduction to tridiagonal form.
    private final void tred2 () {
	for (int j: id0n) d[j] = V[n_1][j];
	i = n;
	while ( --i > 0 ) {
	    // Scale to avoid under/overflow.
	    i_1 = i-1;
	    id0i = IntStream.range(0, i_1).toArray(); //## NOTE: [0,i[
	    scale = h = 0;
	    for (int k: id0i) scale += Math.abs(d[k]);
	    if ( scale == 0.0 ) {
		e[i] = d[i_1];
		for (int j: id0i) {
		    d[j] = V[i_1][j];
		    V[i][j] = V[j][i] = 0.0;
		}
	    }
	    else {
		// Generate Householder vector.
		for (int k: id0i) {
		    d[k] /= scale;
		    h += square(d[k]);
		}
		f = d[i_1];
		g = ( f > 0 ) ? -Math.sqrt(h) : Math.sqrt(h);
		e[i] = scale * g;
		h -= f * g;
		d[i_1] = f - g;
		for (int j: id0i) e[j] = 0.0;
		// Apply similarity transformation to remaining columns.
		for (int j: id0i) {
		    V[j][i] = f = d[j];
		    g = e[j] + V[j][j] * f;
		    for (int k: IntStream.range(j+1, i_1).toArray()) {
			g += V[k][j] * d[k];
			e[k] += V[k][j] * f;
		    }
		    e[j] = g;
		}
		f = 0;
		for (int j: id0i) {
		    e[j] /= h;
		    f += e[j] * d[j];
		}
		hh = f / (2*h);
		for (int j: id0i) e[j] -= hh * d[j];
		for (int j: id0i) {
		    f = d[j];
		    g = e[j];
		    for (int k: IntStream.range(j, i_1).toArray()) V[k][j] -= (f * e[k] + g * d[k]);
		    d[j] = V[i_1][j];
		    V[i][j] = 0;
		}
	    }
	    d[i] = h;
	}
	// Accumulate transformations
	for (int i: id0n_1) {
	    V[n_1][i] = V[i][i];
	    V[i][i] = 1;
	    i1 = i+1;
	    id0i1 = IntStream.range(0, i1).toArray(); //## NOTE [0,i]
	    if ( (h=d[i1]) != 0.0 ) {
		for (int k: id0i1) d[k] = V[k][i1] / h;
		for (int j: id0i1) {
		    g = 0;
		    for (int k: id0i1) g += V[k][i1] * V[k][j];
		    for (int k: id0i1) V[k][j] -= g * d[k];
		}
	    }
	    for (int k: id0i1) V[k][i1] = 0;
	}
	for (int j: id0n) {
	    d[j] = V[n_1][j];
	    V[n_1][j] = 0;
	}
	V[n_1][n_1] = 1;
	e[0] = 0.0;
    }
    
    //## Symmetric tridiagonal QL algorithm.
    private final void tql2 () {
	for (int i: id1n) e[i-1] = e[i];
	e[n_1] = 0.0;
	f = tst1 = 0;
	for (int l: id0n) {
	    // Find small subdiagonal element
	    tst1 = Math.max(tst1, Math.abs(d[l])+Math.abs(e[l]));
	    m = l; while ( m < n && Math.abs(e[m]) > EPS*tst1 ) ++m;
	    // If m == l, d[l] is an eigenvalue,   otherwise, iterate.
	    if ( m > l) {
		iter = 0;
		do {
		    ++iter;
		    // Compute implicit shift
		    g = d[l];
		    p = (d[l+1] - g) / (2.0 * e[l]);
		    r = ( p < 0 ) ? -hypot(p, 1.0) : hypot(p, 1.0);
		    d[l] = e[l] / (p + r);
		    d[l+1] = e[l] * (p + r);
		    h = g - d[l];
		    dl1 = d[(i=(l+1))];
		    while ( ++i < n ) d[i] -= h;
		    f = f + h;
		    // Implicit QL transformation
		    p = d[m];
		    c = 1.0;
		    c2 = c3 = c;
		    el1 = e[l+1];
		    s = s2 = 0.0;
		    i = m;
		    while ( --i >= l ) {
			i1 = i+1;
			c3 = c2;
			c2 = c;
			s2 = s;
			g = c * e[i];
			h = c * p;
			r = hypot(p, e[i]);
			e[i1] = s * r;
			s = e[i] / r;
			c = p / r;
			p = c * d[i] - s * g;
			d[i1] = h + s * (c * g + s * d[i]);
			// Accumulate transformation
			for (int k: id0n) {
			    h = V[k][i1];
			    V[k][i1] = s * V[k][i] + c * h;
			    V[k][i] = c * V[k][i] - s * h;
			}
		    }
		    p = -s * s2 * c3 * el1 * e[l] / dl1;
		    e[l] = s * p;
		    d[l] = c * p;
		} while ( Math.abs(e[l]) > EPS*tst1 );
	    }
	    d[l] = d[l] + f;
	    e[l] = 0;
	}
	// Sort eigenvalues and corresponding vectors
	for (int i: id0n_1) {
	    p = d[k=i];
	    for (int j: IntStream.range(i+1, n).toArray()) //## NOTE: [i+1,n[
		if ( d[j] < p ) {
		    k = j;
		    p = d[j];
		}
	    if ( k != i ) {
		d[k] = d[i];
		d[i] = p;
		for (int j: id0n) {
		    p = V[j][i];
		    V[j][i] = V[j][k];
		    V[j][k] = p;
		}
	    }
	}
    }
    
    //## Nonsymmetric reduction to Hessenberg form.
    private final void orthes () {
	final int low = 0,
	    low1 = low+1,
	    high = n_1,
	    high1 = high+1;

	for (int m: IntStream.range(low1, high).toArray()) { //## NOTE: [low1,high[
	    idmhigh1 = IntStream.range(m, high1).toArray();  //## NOTE: [m,high]
	    m_1 = m-1;
	    // Scale column.
	    scale = 0.0; for (int i: idmhigh1) scale += Math.abs(H[i][m_1]);
	    if ( scale != 0.0 ) {
		// Compute Householder transformation.
		h = 0.0; i = high1; while ( --i >= m ) h += square((ort[i]=H[i][m_1]/scale));
		g = ( ort[m] > 0 ) ? -Math.sqrt(h) : Math.sqrt(h);
		h -= ort[m] * g;
		ort[m] -= g;
		// Apply Householder similarity transformation:  H = (I-u*u'/h)*H*(I-u*u')/h)
		for (int j: IntStream.range(m, n).toArray()) {     //## NOTE: [m,n[
		    f = 0.0; i = high1; while ( --i >= m ) f += ort[i]*H[i][j];
		    f /= h;  for (int i: idmhigh1) H[i][j] -= f*ort[i];
		}
		for (int i: IntStream.range(0, high1).toArray()) {  //## NOTE: [0,high]
		    f = 0.0; j = high1; while ( --j >= m ) f += ort[j]*H[i][j];
		    f /= h;  for (int j: idmhigh1) H[i][j] -= f*ort[j];
		}
		ort[m] = scale * ort[m];
		H[m][m_1] = scale*g;
	    }
	}
	// Accumulate transformations (Algol's ortran).
	for (int i: id0n) { Arrays.fill(V[i], 0); V[i][i] = 1.0; }
	m = high;
	while ( --m >= low1 ) 
	    if ( H[m][(m_1=(m-1))] != 0.0 ) {
		idmhigh1 = IntStream.range(m, high1).toArray();  //## NOTE: [m,high]
		for (int i: idmhigh1) if ( i != m ) ort[i] = H[i][m_1];
		for (int j: idmhigh1) {
		    g = 0;                    for (int i: idmhigh1) g += ort[i] * V[i][j];
		    g = (g/ort[m])/H[m][m_1]; for (int i: idmhigh1) V[i][j] += g * ort[i];
		}
	    }
    }
    
    //## Complex scalar division.
    private transient double cdivr, cdivi;
    private void cdiv(double xr, double xi, double yr, double yi) {
	double rr,dd;
	if ( Math.abs(yr) > Math.abs(yi) ) {
	    rr = yi/yr;
	    dd = yr + rr*yi;
	    cdivr = (xr + rr*xi)/dd;
	    cdivi = (xi - rr*xr)/dd;
	}
	else {
	    rr = yr/yi;
	    dd = yi + rr*yr;
	    cdivr = (rr*xr + xi)/dd;
	    cdivi = (rr*xi - xr)/dd;
	}
    }
    
    //## Nonsymmetric reduction from Hessenberg to real Schur form.
    private final void hqr2 () {
	// Initialize
	nn = this.n;
	n = nn-1;    // ########## WARNING!
	n_1 = n-1;   // ########## WARNING!
	final int low = 0,
	    high = nn-1;
	exshift = p = q = r = s = z = 0;
	// Store roots isolated by balanc and compute matrix norm
	norm = 0;
	for (int i: id0n) {
	    if ( i < low || i > high ) {
		d[i] = H[i][i];
		e[i] = 0.0;
	    }
	    j = ( i-1 > 0 ) ? i-2 : -1; while ( ++j < nn ) norm += Math.abs(H[i][j]);
	}
	// Outer loop over eigenvalue index
	iter = 0;
	while ( n >= low ) {
	    // Look for single small sub-diagonal element
	    l = n;
	    while ( l > low ) {
		s = ( (s=(Math.abs(H[l-1][l-1])+Math.abs(H[l][l]))) == 0 ) ? norm : s;
		if ( Math.abs(H[l][l-1]) < EPS*s ) break;
		--l;
	    }
	    // Check for convergence
	    if ( l == n ) {           // One root found
		H[n][n] += exshift;
		d[n] = H[n][n];
		e[n] = iter = 0;
		--n;    // ########## WARNING!
		--n_1;  // ########## WARNING!
	    }
	    else
		if ( l == n_1 )  {   // Two roots found
		    w = H[n][n_1] * H[n_1][n];
		    p = (H[n_1][n_1] - H[n][n]) / 2.0;
		    q = p*p + w;
		    z = Math.sqrt(Math.abs(q));
		    H[n][n] += exshift;
		    H[n_1][n_1] += exshift;
		    x = H[n][n];
		    if ( q >= 0 ) {   // Real pair
			z = ( p >= 0 ) ? p + z : p - z;
			d[n_1] = x + z;
			d[n] = ( z != 0.0 ) ? x - w/z : d[n_1];
			e[n_1] = e[n] = 0;
			x = H[n][n_1];
			s = Math.abs(x) + Math.abs(z);
			p = x/s;
			q = z/s;
			r = Math.sqrt(p*p + q*q);
			p /= r;
			q /= r;
			// Row modification
			j = n-2;
			while ( ++j < nn ) {
			    z = H[n_1][j];
			    H[n_1][j] = q*z + p*H[n][j];
			    H[n][j] = q*H[n][j] - p*z;
			}
			// Column modification
			i = -1;
			while ( ++i <= n ) {
			    z = H[i][n_1];
			    H[i][n_1] = q*z + p*H[i][n];
			    H[i][n] = q*H[i][n] - p*z;
			}
			// Accumulate transformations
			i = low-1;
			while ( ++i <= high ) {
			    z = V[i][n_1];
			    V[i][n_1] = q*z + p*V[i][n];
			    V[i][n] = q*V[i][n] - p*z;
			} 
		}
		else {   // Complex pair
		    d[n_1] = x+p;
		    d[n] = x+p;
		    e[n_1] = z;
		    e[n] = -z;
		}
		n -= 2;   //############### WARNING!!
		n_1 -= 2; //############### WARNING!!
		iter = 0;
	    }
	    else {                   // No convergence yet
		// Form shift
		x = H[n][n];
		y = ( l < n ) ? H[n_1][n_1] : 0;
		w = ( l < n ) ? H[n][n_1]*H[n_1][n] : 0;
		if ( iter == 10 ) { // Wilkinson's original ad hoc shift
		    exshift += x;
		    i = low-1; while ( ++i <= n ) H[i][i] -= x;
		    s = Math.abs(H[n][n_1]) + Math.abs(H[n_1][n-2]);
		    x = y = 0.75*s;
		    w = -0.4375*s*s;
		}
		if ( iter == 30 && (s=(square((y-x)/2.0) + w)) > 0 ) { // MATLAB's new ad hoc shift
		    s = ( y < x ) ? -Math.sqrt(s) : Math.sqrt(s);
		    s = x - w / ((y-x)/2.0 + s);
		    i = low-1; while ( ++i <= n ) H[i][i] -= s;
		    exshift += s;
		    x = y = w = 0.964;
		}
		++iter;
		// Look for two consecutive small sub-diagonal elements
		m = n-2;
		while ( m >= l ) {
		    m1 = m+1;
		    z = H[m][m];
		    r = x - z;
		    s = y - z;
		    p = (r*s-w)/H[m1][m] + H[m][m1];
		    q = H[m1][m1] - z - r - s;
		    r = H[m+2][m1];
		    s = Math.abs(p) + Math.abs(q) + Math.abs(r);
		    p /= s;
		    q /= s;
		    r /= s;
		    if ( m == l ) break;
		    if ( Math.abs(H[m][m-1])*(Math.abs(q)+Math.abs(r)) < EPS*Math.abs(p)*(Math.abs(H[m-1][m-1])+Math.abs(z)+Math.abs(H[m1][m1])) ) break;
		    --m;
		}
		i = m+1;
		while ( ++i <= n ) {
		    H[i][i-2] = 0.0;
		    if ( i > m+2 ) H[i][i-3] = 0.0;
		}
		// Double QR step involving rows l:n and columns m:n
		k = m-1;
		while ( ++k < n ) {
		    notlast = ( k != n_1 );
		    if ( k != m ) {
			p = H[k][k-1];
			q = H[k+1][k-1];
			r = ( notlast ) ? H[k+2][k-1] : 0;
			if ( (x=(Math.abs(p)+Math.abs(q)+Math.abs(r))) == 0.0 ) continue;
			p /= x;
			q /= x;
			r /= x;
		    }
		    s = ( p < 0 ) ? -Math.sqrt(p*p + q*q + r*r): Math.sqrt(p*p + q*q + r*r);
		    if ( s != 0 ) {
			if ( k != m ) H[k][k-1] = -s*x; else if (l != m) H[k][k-1] = -H[k][k-1];
			p += s;
			x = p/s;
			y = q/s;
			z = r/s;
			q /= p;
			r /= p;
			// Row modification
			j = k-1;
			while ( ++j < nn ) {
			    p = H[k][j] + q*H[k+1][j];
			    if ( notlast ) {
				p += r * H[k+2][j];
				H[k+2][j] -= p*z;
			    }
			    H[k][j] -= p*x;
			    H[k+1][j] -= p*y;
			}
			// Column modification
			i = -1;
			while ( ++i <= ((n<k+3)?n:k+3) ) {
			    p = x*H[i][k] + y*H[i][k+1];
			    if ( notlast ) {
				p += z*H[i][k+2];
				H[i][k+2] -= p*r;
			    }
			    H[i][k] -= p;
			    H[i][k+1] -= p*q;
			}
			// Accumulate transformations
			i = low-1;
			while ( ++i <= high ) {
			    p = x*V[i][k] + y*V[i][k+1];
			    if ( notlast ) {
				p += z*V[i][k+2];
				V[i][k+2] -= p*r;
			    }
			    V[i][k] -= p;
			    V[i][k+1] -= p*q;
			}
		    }  // (s != 0)
		}  // k loop
	    }  // check convergence
	}  // while (n >= low)

	// Backsubstitute to find vectors of upper triangular form
	if ( norm == 0.0 ) return;
	n = nn;              //############### WARNING!!
	n_1 = n-1;           //############### WARNING!!
	while ( --n >= 0 ) { //############### WARNING!!
	    --n_1;           //############### WARNING!!
	    p = d[n];
	    q = e[n];
	    // Real vector
	    if ( q == 0 ) {
		l = n;
		H[n][n] = 1;
		i = n;
		while ( --i >= 0 ) {
		    w = H[i][i] - p;
		    r = 0; for (int j: IntStream.range(l, n+1).toArray()) r += H[i][j] * H[j][n];
		    if ( e[i] < 0 ) {
			z = w;
			s = r;
		    }
		    else {
			l = i;
			if ( e[i] == 0 ) H[i][n] = ( w != 0 ) ? -r/w : -r/(EPS*norm);
			else {  // Solve real equations
			    i1 = i+1;
			    x = H[i][i1];
			    y = H[i1][i];
			    q = square(d[i] - p) + square(e[i]);
			    t = (x*s - z*r)/q;
			    H[i][n] = t;
			    H[i1][n] = ( Math.abs(x) > Math.abs(z) ) ? -(r + w*t)/x : -(s + y*t)/z;
			}
			// Overflow control
			t = Math.abs(H[i][n]);
			if ( (EPS * t) * t > 1 ) for (int j: IntStream.range(i, n+1).toArray()) H[j][n] /= t; 
		    }
		}
	    }
	    else
		if ( q < 0 ) {   // Complex vector
		    l = n_1;
		    // Last vector component imaginary so matrix is triangular
		    if ( Math.abs(H[n][n_1]) > Math.abs(H[n_1][n]) ) {
			H[n_1][n_1] = q / H[n][n_1];
			H[n_1][n] = -(H[n][n] - p) / H[n][n_1]; }
		    else {
			cdiv(0.0, -H[n_1][n], H[n_1][n_1]-p, q);  //################################ TO BE REIMPLEMENTED
			H[n_1][n_1] = cdivr;
			H[n_1][n] = cdivi;
		    }
		    H[n][n_1] = 0;
		    H[n][n] = 1;
		    i = n_1;
		    while ( --i >= 0 ) {
			ra = sa = 0;
			for (int j: IntStream.range(l, n+1).toArray()) { //##NOTE: [l,n]
			    ra += H[i][j]*H[j][n_1];
			    sa += H[i][j]*H[j][n];
			}
			w = H[i][i] - p;
			if ( e[i] < 0 ) { z = w; r = ra; s = sa; }
			else {
			    l = i;
			    if ( e[i] == 0 ) {
				cdiv(-ra, -sa, w, q);  //################################ TO BE REIMPLEMENTED
				H[i][n_1] = cdivr;
				H[i][n] = cdivi;
			    }
			    else {  // Solve complex equations
				i1 = i+1;
				x = H[i][i1];
				y = H[i1][i];
				vr = square(d[i] - p) + square(e[i]) - q*q;
				vi = (d[i] - p) * 2.0 * q;
				if ( vr == 0.0 && vi == 0.0 ) vr = EPS * norm * (Math.abs(w) + Math.abs(q) + Math.abs(x) + Math.abs(y) + Math.abs(z));
				cdiv(x*r-z*ra+q*sa, x*s-z*sa-q*ra, vr, vi); //################################ TO BE REIMPLEMENTED
				H[i][n_1] = cdivr;
				H[i][n] = cdivi;
				if ( Math.abs(x) > (Math.abs(z) + Math.abs(q)) ) {
				    H[i1][n_1] = (-ra - w*H[i][n-1] + q*H[i][n]) / x;
				    H[i1][n] = (-sa - w*H[i][n] - q*H[i][n_1]) / x;
				}
				else {
				    cdiv(-r-y*H[i][n-1], -s-y*H[i][n], z, q); //################################ TO BE REIMPLEMENTED
				    H[i1][n_1] = cdivr;
				    H[i1][n] = cdivi;
				}
			    }
			    // Overflow control
			    t = Math.max(Math.abs(H[i][n_1]), Math.abs(H[i][n]));
			    if ( (EPS * t) * t > 1 )
				for (int j: IntStream.range(i, n+1).toArray()) { //##NOTE: [i,n]
				    H[j][n_1] /= t;
				    H[j][n] /= t;
				}
			}
		    }
		}
	}
	// Vectors of isolated roots
	for (int i: id0n)
	    if ( i < low || i > high )
		for (int j: IntStream.range(i, nn+1).toArray()) V[i][j] = H[i][j];
	// Back transformation to get eigenvectors of original matrix
	j = nn;
	while ( --j >= low ) {
	    for (int i: IntStream.range(low, nn).toArray()) { //## NOTE: [low,high] because high = nn-1
		z = 0; for (int k: IntStream.range(low, ((j<high)?j+1:nn)).toArray()) z += V[i][k] * H[k][j];
		V[i][j] = z;
	    }
	}
	this.n = nn;
	this.n_1 = --nn;
    }

    //## Inversing eigenvector matrix V
    /*
    private final void inv() {
	//## Initializing
	nn = 2*this.n;
	H = new double[n][nn];
	for (int i: id0n ) {
	    H[i][n+i] = 1;
	    for (int j: id0n) H[i][j] = V[i][j];
	}
	//## Gauss-Jordan
	for (int i: id1n) {
	    i_1 = i-1;
	    if ( H[i_1][0] < H[i][0] ) {
		j = nn;
		while ( --j >= 0 ) {
		    r = H[i][j];
		    H[i][j] = H[i_1][j];
		    H[i_1][j] = r;
		}
	    }
	}
	for (int i: id0n)
	    for (int j: id0n)
		if ( i != j ) {
		    r = H[j][i]/H[i][i];
		    k = nn; while ( --k >= 0 ) H[j][k] -= r*H[i][k];
		}
	for (int i: id0n) {
	    r = H[i][i];
	    for (int j: id0n) inV[i][j] = H[i][n+j] / r;
	}
    }
    */
   
    public static void main(String[] args) {
	int size = 4;
	double[][] mp = new double[size][size];
	mp[0][0] = 0.9513; mp[0][1] = 0.0040; mp[0][2] = 0.0430; mp[0][3] = 0.0017;
	mp[1][0] = 0.0040; mp[1][1] = 0.9092; mp[1][2] = 0.0003; mp[1][3] = 0.0865;
	mp[2][0] = 0.0995; mp[2][1] = 0.0008; mp[2][2] = 0.8989; mp[2][3] = 0.0008;
	mp[3][0] = 0.0030; mp[3][1] = 0.1036; mp[3][2] = 0.0004; mp[3][3] = 0.8940;
	Eigen e = new Eigen(mp);
	double[] rev =  e.getRealEigenvalues ();
	double[] iev =  e.getImagEigenvalues ();
	double[][] ev = e.getEigenvectors();
	/*double[][] nev = e.getInvEigenvectors();*/
	System.out.println(Arrays.toString(rev));
	System.out.println(Arrays.toString(iev));
	int row = -1; while ( ++row < size ) System.out.println(Arrays.toString(ev[row]));
	/*row = -1; while ( ++row < size ) System.out.println(Arrays.toString(nev[row]));*/
    }

    
}
