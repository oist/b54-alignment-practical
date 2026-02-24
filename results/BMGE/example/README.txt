Each of the two files prmA.faa and prmA.fna contains a multiple sequence alignment in FASTA format. 
The file prmA.faa was built from cyanobacterial prmA amino acid sequences, and the file prmA.fna 
corresponds to their back-translated nucleotide sequences. The other files prmA.* were obtained 
with the following command lines:

java -jar BMGE.jar -i prmA.faa -t AA -m BLOSUM95      -o prmA.blosum95.faa -oh prmA.blosum95.html
java -jar BMGE.jar -i prmA.fna -t NT -m DNAPAM150:2   -o prmA.pam150.fna   -oh prmA.pam150.html
java -jar BMGE.jar -i prmA.faa -t AA                -oco prmA.blosum30.fna -oh prmA.blosum30.html
