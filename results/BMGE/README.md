[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/)
[![JAVA](https://img.shields.io/badge/Java-11-be0032?logo=java)](http://www.oracle.com/technetwork/java/javase/downloads/index.html)
[![publication](https://img.shields.io/badge/doi-10.1186/1471.2148.10.210-76B900)](https://doi.org/10.1186/1471-2148-10-210)


# BMGE

_BMGE_ (Block Mapping and Gathering with Entropy) is a command line program written in [Java](https://docs.oracle.com/en/java/) to select regions in a multiple sequence alignment that are suited for phylogenetic inference. For more details, see the associated publication (Criscuolo and Gribaldo 2010).

Since _BMGE_ v2.0, this GitLab repository replaces the [previous ftp repository](ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/).

_BMGE_ v2.0 has major differences compared to the previous version (1.12):

+ complete reimplemention of the source code to obtain signicantly faster running times (especially when compiled with [GraalVM](https://www.graalvm.org/)),

+ no longer need of the matrix package [JAMA](https://math.nist.gov/javanumerics/jama/),

+ fixed bugs,

+ simplified/intuitive options (e.g. `-h` is replaced with `-e` to set entropy value thresholds),

+ modified default values (i.e.. `-b`, `-g`, `-m`) to prevent a too strong degradation of the overall phylogenetic signal in standard usage (as highlighted by e.g. Tan et al. 2015),

+ new verbose option (`-v`) to print each character entropy value in tab-delimited format,

+ stationary-based character trimming no longer supported,

+ gap-based filtering no longer supported for sequence (only for characters),

+ reformatting of the NCBI-formatted FASTA headers no longer supported,

+ lighter HTML output files


## Installation

Clone this repository with the following command line:
```bash
git clone https://gitlab.pasteur.fr/GIPhy/BMGE.git
```


## Compilation and execution

The source code of _BMGE_ is inside the _src_ directory. It requires **Java 11** (or higher) and can be compiled and executed in two different ways.

#### Building an executable jar file

On computers with [Oracle JDK](http://www.oracle.com/technetwork/java/javase/downloads/index.html) (**11** or higher) installed, a Java executable jar file can be created. In a command-line window, go to the _src_ directory and type:
```bash
javac BMGE.java
echo Main-Class: BMGE > MANIFEST.MF
jar -cmvf MANIFEST.MF BMGE.jar BMGE.class bmge/*.class
rm MANIFEST.MF BMGE.class bmge/*.class
```
This will create the executable jar file `BMGE.jar` that can be run with the following command line model:
```bash
java -jar BMGE.jar [options]
```

#### Building a native code binary

On computers with [GraalVM](https://www.graalvm.org/downloads/) installed, a native executable can be built. In a command-line window, go to the _src_ directory, and type:
```bash
javac BMGE.java
native-image BMGE BMGE
rm BMGE.build_artifacts.txt BMGE.class bmge/*.class
```
This will create the native executable `BMGE` that can be run with the following command line model:
```bash
./BMGE [options]
```

## Usage

Run _BMGE_ without option to read the following documentation:
```
 Block Mapping and Gathering with Entropy
 Criscuolo and Gribaldo (2010) doi:10.1186/1471-2148-10-210
 https://research.pasteur.fr/software/bmge-block-mapping-and-gathering-with-entropy

 USAGE:  BMGE  -i <infile>  -t <datatype>  -o <outfile>  [options]

 OPTIONS:
   -i <file>           multiple sequence alignment file in FASTA or PHYLIP sequential format (mandatory)
   -t <AA|CO|NT>       input data type; AA: amino acid, CO: codon, NT: nucleotide (mandatory)
   -o[<string>] <file> output file  name in  FASTA (-o, -of),  NEXUS (-ox),  PHYLIP sequential (-op), or
                       HTML (-oh) format;  character state conversion  can be performed by adding suffix
                       aa (amino acid), co (codon), nt (nucleotide) or ry (RY coding);  codon position p
                       selected by adding suffix 1 (p=1), 2 (p=2) and/or 3 (p=3)
   -m BLOSUM<int>      [AA, CO] name of the BLOSUM matrix;  n = 30, 35, 40, ..., 60, 62, 65, ..., 90, 95
                       (default: BLOSUM30)
   -m DNAPAM<int:real> [NT] name of the DNA PAMn matrix (n > 0) with transition/transversion ratio t > 0
                       (default: DNAPAM180:2)
   -m DNAPAM<int>      [NT] name of the DNA PAMn matrix with transition/transversion ratio t = 1
   -m ID               [AA, CO, NT] identity matrix to compute Shannon instead of von Neumann entropy
   -w <int>            sliding window (odd) size for smoothing entropy values (default: 3)
   -e <real>           maximum entropy value threshold (default: 0.5)
   -e <real:real>      minimum and maximum entropy value thresholds, respectively (default: 0:0.5)
   -g <real>           maximum gap rate allowed per character (default: 0.5)
   -b <int>            minimum width of the gathered blocks (default: 3)
   -v                  verbose mode
   -h                  prints this help and exits

 EXAMPLES:
   BMGE  -i dna.fasta -t NT -m DNAPAM100:3 -b 1   -o   out.fna  -ory out.ry.fasta
   BMGE  -i msa.faa   -t AA -m BLOSUM30    -g 0.1 -oh  out.html -opco out.phy
   BMGE  -i codon.phy -t CO -m BLOSUM50    -e 0.4 -oaa out.faa  -onco12 out.nex
```


## Notes

+ For any detail about the character filtering method implemented by BMGE, see the associated paper (Criscuolo and Gribaldo 2010).

+ By default, _BMGE_ is expected to gather character blocks that are well-suited for phylogenetic inference without altering the overall phylogenetic signal. However, to modify the character filtering stringency, it is highly recommended to set the option `-m` with dedicated similarity matrices, instead of modifying the entropy threshold (i.e. option `-e`; see e.g. Steenwyk et al. 2020). When dealing with nucleotide sequences (`-t NT`), the stringency will increase (i.e. less characters are gathered) when lowering the DNAPAM matrix value, e.g. from DNAPAM180 (default) to DNAPAM50. For amino acid or codon sequences, the stringency will increase when setting large BLOSUM matrix values, e.g. from BLOSUM30 (default) to BLOSUM90.

+ An alternative to increase the stringency of _BMGE_ (at the cost of some reduction of the phylogenetic signal) is to gather only large blocks of conserved characters by setting the option `-b` with large values (e.g. 10).

+ _BMGE_ can also be used to convert input files in FASTA or PHYLIP (sequential) format into FASTA (option `of`), NEXUS (`-ox`) or PHYLIP sequential (`-op`) format. HTML output files (`-oh`) can be also written to have a look at the multiple sequence alignment after character filtering. In complement, translation, back-translation or codon position selection can be easily performed by adding dedicated suffixes to the option `-o`. Of important note, to only perform such transformations, run _BMGE_ with options `-e 1 -g 1` to unset character filtering steps.


## References

Criscuolo A, Gribaldo S (2010) _BMGE (Block Mapping and Gathering with Entropy): a new software for selection of phylogenetic informative regions from multiple sequence alignments_. **BMC Evolutionary Biology**, 10:210. [doi:10.1186/1471-2148-10-210](https://doi.org/10.1186/1471-2148-10-210)

Steenwyk JL, Buida TJ III, Li Y, Shen X-X, Rokas A (2020) _ClipKIT: A multiple sequence alignment trimming software for accurate phylogenomic inference_. **PLoS Biology**, 18(12):e3001007. [doi:10.1371/journal.pbio.3001007](https://doi.org/10.1371/journal.pbio.3001007)

Tan G, Muffato M, Ledergerber C, Herrero J, Goldman N, Gil M, Dessimoz C (2015) _Current Methods for Automated Filtering of Multiple Sequence Alignments Frequently Worsen Single-Gene Phylogenetic Inference_. **Systematic Biology**, 64(5):778-791. [doi:10.1093/sysbio/syv033](https://doi.org/10.1093/sysbio/syv033)


