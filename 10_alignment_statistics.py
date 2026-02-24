#!/usr/bin/env python3
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq

def alignment_stats(aln_path, label=''):
    """Print summary statistics for an alignment."""
    if 'clustal' in aln_path:
        aln = AlignIO.read(aln_path, 'clustal')
    else:
        aln = AlignIO.read(aln_path, 'fasta')
    ncol = aln.get_alignment_length()
    nseq = len(aln)
    total_gaps = sum(str(r.seq).count('-') for r in aln)
    total_chars = nseq * ncol
    gap_pct = 100 * total_gaps / total_chars

    # Count fully conserved columns
    conserved = 0
    for j in range(ncol):
        col = set(str(aln[i].seq)[j] for i in range(nseq))
        col.discard('-')
        if len(col) == 1 and '-' not in {str(aln[i].seq)[j] for i in range(nseq)}:
            conserved += 1

    return {'label': label, 'ncol': ncol, 'gaps': total_gaps,
            'gap_pct': gap_pct, 'conserved': conserved}

if __name__ == "__main__":
    alignments = ["pax6_family.clustalw.fasta", "pax6_family.linsi.fasta",  "pax6_family.muscle.fasta",  "pax6_family.prank.fasta.best.fas"]
    for f in alignments:
        print(alignment_stats(f, label=f))
