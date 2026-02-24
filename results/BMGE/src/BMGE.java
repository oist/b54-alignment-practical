/*
  ########################################################################################################

  BMGE: Block Mapping and Gathering with Entropy
    
  Copyright (C) 2010-2024  Institut Pasteur
  
  This program  is free software:  you can  redistribute it  and/or modify it  under the terms  of the GNU
  General Public License as published by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,  but WITHOUT ANY WARRANTY;  without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
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

import bmge.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Locale;
import java.util.stream.IntStream;

public class BMGE {

    //### constants   ###############################################################
    static final String VERSION  = "2.0                                                      Copyright (C) 2010-2024 Institut Pasteur";
    static final String NOTHING  = "N.o./.T.h.I.n.G";
    static final byte   NONE     = (byte) 0;
    static final byte   AA       = (byte) 10;
    static final byte   CO       = (byte) 11;
    static final byte   NT       = (byte) 12;
    static final byte   RY       = (byte) 13;
    static final byte   FASTA    = (byte) 20;
    static final byte   HTML     = (byte) 21;
    static final byte   NEXUS    = (byte) 22;
    static final byte   PHYLIP   = (byte) 23;
    static final double EPS      = 1e-5;

    //### options   #################################################################
    static File              infile;   // -i
    static byte format; 
    static ArrayList<String> outforms; // -o
    static ArrayList<File>   outfiles; 
    static String            stype;    // -t
    static byte              datatype; 
    static String            model;    // -m
    static SMatrix           sm;
    static int               wsize;    // -w
    static String            hmm;      // -e
    static double            hmin;
    static double            hmax;
    static double            gmax;     // -g
    static int               minwidth; // -b
    static boolean           verbose;  // -v

    //### io   ######################################################################
    static String filename;
    static BufferedReader in;
    static BufferedWriter out;

    //### data   ####################################################################
    static int n;       // no. sequences
    static int l;       // no. characters
    static MSA msa;     // multiple sequence alignment
    static int[] cid;   // column ids
    static int[] rid;   // row ids
    static double[] g;  // character gap proportion
    static double[] e;  // character entropy
    static double[] h;  // smoothed entropy
    static BitSet mask; // character mask
    static MSA mmsa;    // masked multiple sequence alignment
    static MSA cmsa;    // converted multiple sequence alignment

    //### stuffs   ##################################################################
    static byte of, ot;
    static int b, o, p, t;
    static int hw, c_ws, c_hw, start, end, lb_1, lb, lb1;
    static boolean conserved, ok, oc;
    static double x, up, dn, grate;
    static String line, seq, blank, cs, mcs;
    static StringBuilder sb;
    static ArrayList<String> als1, als2;
    static ArrayList<ArrayList<Integer>> blocks;
    static BitSet bs;
   
    public static void main(String[] args) throws IOException {
	
	//#############################################################################################################
	//#############################################################################################################
	//### reading options                                                                                       ###
	//#############################################################################################################
	//#############################################################################################################
	if ( args.length < 4 ) { usage(); System.exit(1); }
	
	infile  = new File(NOTHING);        // -i
	outforms = new ArrayList<String>(); // -o
	outfiles = new ArrayList<File>();   // -o
	stype = NOTHING;                    // -t
	model = NOTHING;                    // -m
	wsize = 3;                          // -w
	hmm = "0:0.5";                      // -e
	gmax = 0.5;                         // -g
	minwidth = 3;                       // -b
	verbose = false;                    // -v
	o = -1;
	while ( ++o < args.length ) {
	    if ( args[o].equals("-i") )     {   infile = new File(args[++o]);           continue; }
	    if ( args[o].equals("-t") )     {    stype = args[++o];                     continue; }
	    if ( args[o].equals("-m") )     {    model = args[++o];                     continue; }
	    if ( args[o].equals("-w") ) try {    wsize = Integer.parseInt(args[++o]);   continue; } catch ( NumberFormatException e ) { System.err.println("incorrect value (option -w): " + args[o]); System.exit(1); }
	    if ( args[o].equals("-e") )     {      hmm = args[++o];                     continue; }
	    if ( args[o].equals("-g") ) try {     gmax = Double.parseDouble(args[++o]); continue; } catch ( NumberFormatException e ) { System.err.println("incorrect value (option -g): " + args[o]); System.exit(1); }
	    if ( args[o].equals("-b") ) try { minwidth = Integer.parseInt(args[++o]);   continue; } catch ( NumberFormatException e ) { System.err.println("incorrect value (option -b): " + args[o]); System.exit(1); }
	    if ( args[o].equals("-v") )     {  verbose = true;                          continue; }
	    if ( args[o].startsWith("-o") ) { outforms.add(args[o].substring(1)); outfiles.add(new File(args[++o])); continue; }
	    if ( args[o].equals("-h") )     {  usage();                           System.exit(0); }
	}

	//### infile   ###################################################################
	filename = infile.toString();
	if ( filename.equals(NOTHING) )     { System.err.println("no input file (option -i)");                        System.exit(1); }
	if ( ! infile.exists() )            { System.err.println("infile does not exist (option -i): " + filename);   System.exit(1); }

	//### outfile(s)   ###############################################################
	if ( outfiles.size() == 0 )         { System.err.println("no output file (option -o)");                       System.exit(1); }
	of = NONE; for (String form: outforms) if ( form.indexOf("oh") != -1 ) of = HTML;

	//### data type   ################################################################
	if ( stype.equals(NOTHING) )        { System.err.println("no sequence coding (option -t)");                   System.exit(1); }
	switch ( stype.toUpperCase() )  {
	case "AA":
	case "AMINO": datatype = AA;   break;
	case "CO":
	case "CODON": datatype = CO;   break;
	case "DNA":
	case "NT":
	case "RNA":   datatype = NT;   break;
	default:      datatype = NONE; break;
	}
	if ( datatype == NONE )             { System.err.println("incorrect sequence coding (option -t): " + stype);  System.exit(1); }

	//### model   ####################################################################
	if ( model.equals(NOTHING) )
	    switch ( datatype ) {
	    case AA:
	    case CO: model = "BLOSUM30";    break;
	    case NT: model = "DNAPAM180:2"; break;
	    }
	if ( ! (model.startsWith("BLOSUM") ||
		model.startsWith("DNAPAM") ||
		model.startsWith("ID")) )   { System.err.println("incorrect model (option -m): " + model);            System.exit(1); }
	if ( model.startsWith("BLOSUM") ) {
	    if ( datatype == NT )           { System.err.println("incompatible model (option -m): " + model);         System.exit(1); }
	    switch ( model ) {
	    case "BLOSUM30":  sm = Model.BLOSUM30;  break;
	    case "BLOSUM35":  sm = Model.BLOSUM35;  break;
	    case "BLOSUM40":  sm = Model.BLOSUM40;  break;
	    case "BLOSUM45":  sm = Model.BLOSUM45;  break;
	    case "BLOSUM50":  sm = Model.BLOSUM50;  break;
	    case "BLOSUM55":  sm = Model.BLOSUM55;  break;
	    case "BLOSUM60":  sm = Model.BLOSUM60;  break;
	    case "BLOSUM62":  sm = Model.BLOSUM62;  break;
	    case "BLOSUM65":  sm = Model.BLOSUM65;  break;
	    case "BLOSUM70":  sm = Model.BLOSUM70;  break;
	    case "BLOSUM75":  sm = Model.BLOSUM75;  break;
	    case "BLOSUM80":  sm = Model.BLOSUM80;  break;
	    case "BLOSUM85":  sm = Model.BLOSUM85;  break;
	    case "BLOSUM90":  sm = Model.BLOSUM90;  break;
	    case "BLOSUM95":  sm = Model.BLOSUM95;  break;
	    case "BLOSUM100": sm = Model.BLOSUM100; break;
	    default:                          System.err.println("model non available (option -m): " + model);        System.exit(1);
	    }
	}
	if ( model.startsWith("DNAPAM") ) {
	    if ( datatype != NT )           { System.err.println("incompatible model (option -m): " + model);         System.exit(1); }
	    line = model.replaceAll("DNAPAM", "");
	    t = 1;
	    try {
		if ( (o=line.indexOf(":")) < 0 ) p = Integer.parseInt(line);
		else {
		    p = Integer.parseInt(line.substring(0,o));
		    t = Integer.parseInt(line.substring(++o));
		}
	    } catch ( NumberFormatException e ) { System.err.println("incorrect model (option -m): " + model);        System.exit(1); }
	    sm = Model.pam(p, t);
	}
	if ( model.startsWith("ID") )
	    switch ( datatype ) {
	    case AA:
	    case CO: sm = Model.id(20); break;
	    case NT: sm = Model.id(4);  break;
	    }

	//### window size   ##############################################################
	if ( wsize % 2 == 0 )		    { System.err.println("window size should be odd (option -w): " + wsize);  System.exit(1); }

	//### entropy threshold(s)   #####################################################
	hmin = 0;
	try {
	    if ( (o=hmm.indexOf(":")) < 0 ) hmax = Double.parseDouble(hmm);
	    else {
		hmin = Double.parseDouble(hmm.substring(0,o));
		hmax = Double.parseDouble(hmm.substring(++o));
	    }
	} catch ( NumberFormatException e ) { System.err.println("incorrect value (option -h): " + hmm);              System.exit(1); }
	

	//#############################################################################################################
	//#############################################################################################################
	//### reading infile                                                                                        ###
	//#############################################################################################################
	//#############################################################################################################
	filename = infile.toString();
	in = Files.newBufferedReader(Path.of(filename));

	//### determining format   #######################################################
	line = "";
	while ( line.length() == 0 ) try { line = in.readLine().trim(); } catch ( NullPointerException e ) { break; }
	format = ( line.length() == 0 ) ? NONE : ( line.startsWith(">") ) ? FASTA : PHYLIP;

	als1 = new ArrayList<String>(); // labels
	als2 = new ArrayList<String>(); // sequences
	l = 0;                          // MSA length
	switch ( format ) {
	    
	//### reading FASTA file   #######################################################
	case FASTA:
	    als1.add(line.substring(1).trim());
	    sb = new StringBuilder("");
	    while ( true ) {
		try { line = in.readLine().trim(); }  catch ( NullPointerException e ) { break; }
		if ( line.startsWith(">") ) {
		    seq = sb.toString().replaceAll(" ", "");
		    if ( l > 0 && seq.length() != l ) { System.err.println("problem with input file: different sequence length"); System.exit(1); }
		    l = seq.length();
		    als2.add(seq);
		    sb = new StringBuilder("");
		    als1.add(line.substring(1).trim());
		    continue;
		}
		sb = sb.append(line);
	    }
	    seq = sb.toString().replaceAll(" ", "");
	    if ( l > 0 && seq.length() != l ) { System.err.println("problem with input file: different sequence length"); System.exit(1); }
	    l = seq.length();
	    als2.add(seq);
	    break;
	    
	//### reading PHYLIP sequential file   ###########################################
	case PHYLIP:
	    while ( true ) {
		try { line = in.readLine().trim(); }  catch ( NullPointerException e ) { break; }
		if ( line.length() == 0 ) continue;
		b = line.indexOf(" ");
		als1.add(line.substring(0, b));
		seq = line.substring(b).replaceAll(" ", "");
		if ( l > 0 && seq.length() != l ) { System.err.println("problem with input file: different sequence length"); System.exit(1); }
		l = seq.length();
		als2.add(seq.toString());
	    }
	    break;
	    
	default:
	    System.err.println("problem with input file (option -i): " + filename);
	    System.exit(1);
	}
	in.close();
	
	    
	//#############################################################################################################
	//#############################################################################################################
	//### storing MSA                                                                                           ###
	//#############################################################################################################
	//#############################################################################################################
	n = als1.size();   // MSA size
	line = ((format == FASTA) ? "FASTA input:   " : "PHYLIP input:  " );
	switch ( datatype ) {
	case AA:    msa = new MSAaa(n, l); System.out.println(line + n + " sequences, " + l + " amino acid characters  (" + infile.toString() + ")"); break;
	case NT:    msa = new MSAnt(n, l); System.out.println(line + n + " sequences, " + l + " nucleotide characters  (" + infile.toString() + ")"); break;
	case CO:
	    if ( l % 3 != 0 ) { System.err.println("not codon sequences as their length is not a multiple of 3"); System.exit(1); }
	    l /= 3; msa = new MSAco(n, l); System.out.println(line + n + " sequences, " + l + " codon characters  ("      + infile.toString() + ")"); break;
	}
	cid = IntStream.range(0, l).toArray();
	rid = IntStream.range(0, n).toArray();
	for (int i: rid) {
	    msa.setLabel(i, als1.get(i));
	    msa.setSequence(i, als2.get(i));
	}
	als1.clear();
	als2.clear();

	
	//#############################################################################################################
	//#############################################################################################################
	//### estimating entropy values                                                                             ###
	//#############################################################################################################
	//#############################################################################################################
	g = new double[l];  // character gap proportion
	e = new double[l];  // character entropy
	h = new double[l];  // smoothed entropy
	if ( of != HTML && hmin <= 0 && hmax >= 1 ) {
	    Arrays.fill(g, 1);
	    Arrays.fill(e, 1);
	    Arrays.fill(h, 1);
	}
	else {
	    hw = (wsize - 1) / 2;
	    up = dn = 0;
	    c_ws = -wsize;  // = c - wsize
	    c_hw = -hw - 1; // = c - (wsize - 1) / 2
	    for (int c: cid) {
		g[c] = msa.gapFreq(c);
		e[c] = h[c] = Utils.entropy(msa.getFreq(c), sm);
		x = 1.0 - g[c];  // computing h, i.e. 
		up += x * e[c];  //  formula (3) in
		dn += x;         //  the BMGE paper
		
		++c_ws;
		++c_hw;
		if ( c_ws >= 0 ) {
		    x = ( dn < EPS ) ? 0 : up / dn;                    //## x = wgt avg of e in [c-wsize , c]
		    h[c_hw] = ( x < EPS ) ? 0 : ( x > 1.0 ) ? 1.0 : x; //## smoothed entropy at position c-hw
		    //h[c_hw] = Math.min(e[c_hw], h[c_hw]);
		    x = 1.0 - g[c_ws];  // updating h
		    up -= x * e[c_ws];  // updating h
		    dn -= x;            // updating h
		}
	    }
	}
	

	//#############################################################################################################
	//#############################################################################################################
	//### determining blocks to gather                                                                          ###
	//#############################################################################################################
	//#############################################################################################################
	mask = new BitSet(l); // mask.get(c) == keep c

	//### initial masking   ##########################################################
	for (int c: cid) mask.set(c, h[c] <= hmax);

	//### masking by block merging  ##################################################
	if ( of != HTML && hmin <= 0 && hmax >= 1 ) {
	    blocks = new ArrayList<ArrayList<Integer>>();
	    ok = false;
	    while ( ! ok  ) {
		//## assessing conserved/variable blocks
		blocks.clear();
		b = start = 0; conserved = true; // first block should be conserved
		for (int c: cid) 
		    if ( (conserved && ! mask.get(c)) || (mask.get(c) && ! conserved) ) {
			blocks.add(new ArrayList<Integer>());
			blocks.get(b).add(start);           // field 0: start position
			blocks.get(b).add(start=c);         // field 1: end position
			conserved = ! conserved;
			++b;
		    }
		blocks.add(new ArrayList<Integer>());
		blocks.get(b).add(start);
		blocks.get(b).add(l); 
		if ( ! conserved ) { // last block should be conserved
		    blocks.add(new ArrayList<Integer>());
		    blocks.get(++b).add(l);
		    blocks.get(b).add(l);
		}	
		
		//## attempting to merge variable blocks with neighbor conserved blocks
		ok = true;
		b = blocks.size();
		while ( --b >= 0 ) {
		    // conserved (C) block (always at even positions)
		    if ( b % 2 == 0 ) { 
			mask.set(blocks.get(b).get(0), blocks.get(b).get(1));
			continue; 
		    }
		    // variable (V) block
		    lb_1 = blocks.get(b-1).get(1) - blocks.get(b-1).get(0); // l(C) before V
		    lb   = blocks.get(b).get(1)   - blocks.get(b).get(0);   // l(V) 
		    lb1  = blocks.get(b+1).get(1) - blocks.get(b+1).get(0); // l(C) after V
		    if ( lb < lb_1 && lb < lb1 ) {                          // merging only if l(V) lower than both l(C)
			start = blocks.get(b-1).get(0);
			end = blocks.get(b+1).get(1);
			grate = up = dn = 0;
			p = end;
			while ( --p >= start ) {
			    grate += g[p];
			    up += (1.0 - g[p]) * h[p];
			    dn += (1.0 - g[p]);
			}
			if ( grate/(end-start) < 0.33 &&                      // gap rate < 33%   and
			     up/dn <= hmax ) {                                // wgt entropy < threshold ==> merging(CVC)
			    mask.set(blocks.get(b).get(0), blocks.get(b).get(1));
			    ok = false;
			    continue;
			}
		    }
		    //System.out.print(blocks.get(b).toString());
		}
	    }
	}
	
	//### masking gapped regions  ####################################################
	if ( gmax < 1.0 ) for (int c: cid) if ( g[c] > gmax ) mask.clear(c);
	
	//### masking too conserved regions  #############################################
	if ( hmin > 0.0 ) for (int c: cid) if ( e[c] <= hmin ) mask.clear(c);
	    
	//### masking too short blocks  ##################################################
	if ( minwidth > 0 ) {
	    end = 0;
	    while ( (start=mask.nextSetBit(end)) != -1 ) {
		end = mask.nextClearBit(start);
		if ( end - start < minwidth ) mask.clear(start, end);
	    }
	}

	//### printing details  ##########################################################
	if ( verbose ) {
	    System.out.println("#pos.\tentr.\tsmooth\tgap\tselected");
	    for (int c: cid) System.out.println(String.format(Locale.US, "%d\t%.3f\t%.3f\t%.3f\t", (c+1), e[c], h[c], g[c]) + ((mask.get(c))?"+":""));
	}
		

	//#############################################################################################################
	//#############################################################################################################
	//### masking MSA                                                                                           ###
	//#############################################################################################################
	//#############################################################################################################
	mmsa = msa.trim(mask);

	
	//#############################################################################################################
	//#############################################################################################################
	//### writing output file(s)                                                                                ###
	//#############################################################################################################
	//#############################################################################################################
	lb = 0; for (int i: rid) lb = ( lb < (t=mmsa.getLabel(i).length()) ) ? t : lb;
	blank = " ".repeat(2*lb);

	o = -1;
	for (String form: outforms) {
	    //### converting  ############################################################
	    ot = (  form.indexOf("aa") != -1 ) ? AA
		: ( form.indexOf("co") != -1 ) ? CO
		: ( form.indexOf("nt") != -1 ) ? NT
		: ( form.indexOf("ry") != -1 ) ? RY
		: datatype;

	    switch ( ot ) {
	    case AA: cmsa = mmsa.toMSAaa(); break;
	    case CO: cmsa = mmsa.toMSAco(); break;
	    case NT: cmsa = mmsa.toMSAnt(); break;
	    case RY: cmsa = mmsa.toMSAry(); break;
	    }
	    oc = ( ot == CO && (form.indexOf("1") != -1 || form.indexOf("2") != -1 || form.indexOf("3") != -1) );
	    if ( oc ) {
		p = 3 * cmsa.length();
		bs = new BitSet(p);
		if ( form.indexOf("1") != -1 ) { b = -3; while ( (b+=3) < p ) bs.set(b); }
		if ( form.indexOf("2") != -1 ) { b = -2; while ( (b+=3) < p ) bs.set(b); }
		if ( form.indexOf("3") != -1 ) { b = -1; while ( (b+=3) < p ) bs.set(b); }
		ot = NT;
		cmsa = cmsa.toMSAnt().trim(bs);
	    }

	    //### writing  ###############################################################
	    of = (  form.indexOf("of") != -1 ) ? FASTA
		: ( form.indexOf("oh") != -1 ) ? HTML
		: ( form.indexOf("ox") != -1 ) ? NEXUS
		: ( form.indexOf("op") != -1 ) ? PHYLIP
		: FASTA;
	    ot = ( of == HTML ) ? datatype : ot;
	    filename = outfiles.get(++o).toString();
	    out = Files.newBufferedWriter(Path.of(filename));

	    line = ((of == FASTA) ? "FASTA output:  " : (of == HTML) ? "HTML output:   " : (of == NEXUS) ? "NEXUS output:  " : "PHYLIP output: ");
	    lb1 = ( of == HTML ) ? msa.length() : cmsa.length();		
	    switch ( ot ) {
	    case AA: System.out.println(line + cmsa.size() + " sequences, " + lb1 + " amino acid characters  (" + filename + ")"); break;
	    case CO: System.out.println(line + cmsa.size() + " sequences, " + lb1 + " codon characters  ("      + filename + ")"); break;
	    case NT: System.out.println(line + cmsa.size() + " sequences, " + lb1 + " nucleotide characters  (" + filename + ")"); break;
	    case RY: System.out.println(line + cmsa.size() + " sequences, " + lb1 + " RY characters  ("         + filename + ")"); break;
	    }
	    
	    switch ( of ) {
		
	    //### writing FASTA file   ###################################################
	    case FASTA:
		for (int i: rid) {
		    out.write(">" + cmsa.getLabel(i)); out.newLine();
		    out.write(cmsa.getSequence(i));    out.newLine();
		}
		break;
		
	    //### writing NEXUS file   ###################################################
	    case NEXUS:
		for (int i: rid) cmsa.setLabel(i, cmsa.getLabel(i).replace(' ', '_').replace(':', '_').replace(',', '_').replace(';', '_').replace('(', '_').replace(')', '_'));
		out.write("#NEXUS"); out.newLine(); out.newLine();
		out.write("begin data;"); out.newLine();
		lb1 = ( oc ) ? cmsa.length() : 3*cmsa.length();
		out.write("   dimensions ntax=" + cmsa.size() +" nchar=" + lb1 + ";"); out.newLine();
		if ( ot == AA ) out.write("   format datatype=PROTEIN;"); else out.write("   format datatype=NUCLEOTIDE;");
		out.newLine();
		out.write("   matrix"); out.newLine();
		for (int i: rid) { out.write(("      " + cmsa.getLabel(i) + blank).substring(0, lb+6) + " " + cmsa.getSequence(i)); out.newLine(); }
		out.write("   ;"); out.newLine();
		out.write("end;"); out.newLine();
		break;

	    //### writing PHYLIP file   ##################################################
	    case PHYLIP:
		for (int i: rid) cmsa.setLabel(i, cmsa.getLabel(i).replace(' ', '_').replace(':', '_').replace(',', '_').replace(';', '_').replace('(', '_').replace(')', '_'));
		lb1 = ( oc ) ? cmsa.length() : 3*cmsa.length();
		out.write(" " + cmsa.size() + " " + lb1); out.newLine();
		for (int i: rid) { out.write((cmsa.getLabel(i) + blank).substring(0, lb) + " " + cmsa.getSequence(i)); out.newLine(); }
		break;

	    //### writing HTML file   ####################################################
	    case HTML:
		out.write("<html>");                                             out.newLine();
		out.write("<head>");                                             out.newLine();
		out.write("<style>");                                            out.newLine();
		out.write("body   {background-color: #FFFFFF; color: #D3D4D4}"); out.newLine();
		out.write("span.a {color: #000000}");                            out.newLine();
		out.write("span.b {color: #000000; font-weight: bold}");         out.newLine();
		out.write("span.c {background-color: #000000; color: #FFFFFF}"); out.newLine();
		out.write("span.d {background-color: #F0F0F0; color: #000000}"); out.newLine();
		out.write("</style>");                                           out.newLine();
		out.write("</head>");                                            out.newLine();
		out.write("<body>");                                             out.newLine();
		out.write("<pre>");                                              out.newLine();
		out.newLine();                                                   out.newLine();
		//##### histograms #####		
		t = 11;
		while ( --t >= 0 ) {
		    switch (t) {
		    case 10: out.write(blank.substring(0 , lb-3) + "<span class=\"a\">1.0-</span>"); break;
		    case 5:  out.write(blank.substring(0 , lb-3) + "<span class=\"a\">0.5-</span>"); break;
		    case 0:  out.write(blank.substring(0 , lb-3) + "<span class=\"a\">0.0-</span>"); break;
		    default: out.write(blank.substring(0 , lb+1));                                   break;
		    }
		    for (int c: cid) 
			if ( h[c] < t/10.0 ) out.write( ((g[c] <= t/10.0) ? ((datatype == CO)?"   ":" ") : ((datatype == CO)?"===":"=")) );
			else                 out.write("<span class=\"b\">" + ((h[c] > (t+0.5)/10.0) ? ((datatype == CO)?":::":":") : ((datatype == CO)?"...":".")) + "</span>");
		    out.newLine();
		}
		//##### bars #####
		out.newLine(); 
		sb = new StringBuilder("");
		lb1 = (datatype == CO) ? 3*msa.length() : msa.length(); 
		p = -1; while ( ++p <= lb1 ) { if ( p % 5 == 0 ) sb = sb.append(((p%10 == 0)?'|':'-')); else sb = sb.append('='); }
		seq = blank.substring(0, lb+1) + "<span class=\"d\">" + sb.substring(1) + "</span>";
		out.write(seq); out.newLine();
		sb = new StringBuilder(" ".repeat(lb1));
		p = 0; while ( (p+=10) < lb1 ) { line = String.valueOf(p); if ( p+line.length()-1 < lb1 ) sb = sb.replace(p-1, p+line.length()-1, line); }
		line = blank.substring(0, lb+1) + "<span class=\"d\">" + sb.toString() + "</span>";
		out.write(line); out.newLine(); 
		out.write(seq + "   "); out.newLine(); 
		//##### alignment #####
		als1 = new ArrayList<String>(msa.length()); for (int c: cid) als1.add(msa.getMajorityCharState(c));
		for (int i: rid) {
		    out.write("<span class=\"a\">" + (msa.getLabel(i) + blank).substring(0 , lb+1) + "</span>"); 
		    for (int c: cid) {
			cs = msa.getCharState(i, c);
			if ( ! mask.get(c) ) out.write(cs); else out.write( ((cs.equals(als1.get(c))) ? "<span class=\"c\">" : "<span class=\"a\">") + cs + "</span>" );
		    }
		    out.newLine();
		}
		out.newLine();
		//##### characters #####
		out.write("<span class=\"a\">"); out.newLine();
		out.write("Characters : " + mask.cardinality() + " selected  " + (msa.length()-mask.cardinality()) + " removed"); out.newLine();
		out.write("  selected: ");
		end = 0;
		while ( (start=mask.nextSetBit(end)) != -1 ) {
		    end = mask.nextClearBit(start);
		    if ( ++start == end ) out.write(" " + (start)); else out.write(" " + (start) + "-" + end);
		}
		out.newLine();
		out.write("</span>"); out.newLine(); 
		out.write("</pre>");  out.newLine(); 
		out.write("</body>"); out.newLine(); 
		out.write("</html>"); out.newLine(); 
		break;
	    }
	    out.close();
	}
    }

    public static void usage() {
	System.out.println("");
	System.out.println(" BMGE v" + VERSION);
	System.out.println("");
	System.out.println(" Block Mapping and Gathering with Entropy");
	System.out.println(" Criscuolo and Gribaldo (2010) doi:10.1186/1471-2148-10-210");
	System.out.println(" https://research.pasteur.fr/software/bmge-block-mapping-and-gathering-with-entropy");
	System.out.println("");
	System.out.println(" USAGE:  BMGE  -i <infile>  -t <datatype>  -o <outfile>  [options]");
	System.out.println("");
	System.out.println(" OPTIONS:");
	System.out.println("   -i <file>           multiple sequence alignment file in FASTA or PHYLIP sequential format (mandatory)");
	System.out.println("   -t <AA|CO|NT>       input data type; AA: amino acid, CO: codon, NT: nucleotide (mandatory)");
	System.out.println("   -o[<string>] <file> output file  name in  FASTA (-o, -of),  NEXUS (-ox),  PHYLIP sequential (-op), or");
	System.out.println("                       HTML (-oh) format;  character state conversion  can be performed by adding suffix");
	System.out.println("                       aa (amino acid), co (codon), nt (nucleotide) or ry (RY coding);  codon position p");
	System.out.println("                       selected by adding suffix 1 (p=1), 2 (p=2) and/or 3 (p=3)");
	System.out.println("   -m BLOSUM<int>      [AA, CO] name of the BLOSUM matrix;  n = 30, 35, 40, ..., 60, 62, 65, ..., 90, 95");
	System.out.println("                       (default: BLOSUM30)");
	System.out.println("   -m DNAPAM<int:real> [NT] name of the DNA PAMn matrix (n > 0) with transition/transversion ratio t > 0");
	System.out.println("                       (default: DNAPAM180:2)");
	System.out.println("   -m DNAPAM<int>      [NT] name of the DNA PAMn matrix with transition/transversion ratio t = 1");
	System.out.println("   -m ID               [AA, CO, NT] identity matrix to compute Shannon instead of von Neumann entropy");
	System.out.println("   -w <int>            sliding window (odd) size for smoothing entropy values (default: 3)");
	System.out.println("   -e <real>           maximum entropy value threshold (default: 0.5)");
	System.out.println("   -e <real:real>      minimum and maximum entropy value thresholds, respectively (default: 0:0.5)");
	System.out.println("   -g <real>           maximum gap rate allowed per character (default: 0.5)");
	System.out.println("   -b <int>            minimum width of the gathered blocks (default: 3)");
	System.out.println("   -v                  verbose mode");
	System.out.println("   -h                  prints this help and exits");
	System.out.println("");
	System.out.println(" EXAMPLES:");
	System.out.println("   BMGE  -i dna.fasta -t NT -m DNAPAM100:3 -b 1   -o   out.fna  -ory out.ry.fasta");
	System.out.println("   BMGE  -i msa.faa   -t AA -m BLOSUM30    -g 0.1 -oh  out.html -opco out.phy");
	System.out.println("   BMGE  -i codon.phy -t CO -m BLOSUM50    -e 0.4 -oaa out.faa  -onco12 out.nex");
	System.out.println("");
    }
    
}

