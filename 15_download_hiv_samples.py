#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import time

Entrez.email = "EMAILADDRESS"  # since 2010 it's mandatory to give an email address to use NCBI

dentist_sequences = {
    "FLD1": "M90848.1",
    "FLD2": "M90849.1",
    "FLD4": "M90850.1",
    "FLD5": "M90851.1",
    "FLD7": "M90852.1",
    "FLD8": "M90853.1",
    "FLPA1": "M90855.1",
    "FLPA5": "M90857.1",
    "FLPA7": "M90859.1",
    "FLPB2": "M90863.1",
    "FLPB29": "M90871.1",
    "FLPB31": "M90872.1",
    "FLPC12": "M90877.1",
    "FLPC14": "M90878.1",
    "FLPC20": "M90879.1",
    "FLPD1": "M90882.1",
    "FLPD12": "M90886.1",
    "FLPD3": "M90883.1",
    "FLPE4": "M90889.1",
    "FLPE5": "M90890.1",
    "FLPE6": "M90891.1",
    "FLPF1": "M90895.1",
    "FLPF3": "M90896.1",
    "FLPF5": "M90897.1",
    "FLPG1": "M90902.1",
    "FLPG3": "M90903.1",
    "FLPG4": "M90904.1",
    "FLPH1D": "M90907.1",
    "FLPH2D": "M90908.1",
    "FLPH4D": "M90909.1",
    "LC08": "M90965.1",
    "LC35": "M90964.1",
    "LC32": "M90961.1",
    "LC15": "M90944.1",
    "LC06": "M90936.1",
    "FLQ5R3E": "M92130.1",
    "FLQ5R3A": "M92127.1",
    "FLQ5R3C": "M92128.1",
    "FLQ5R3D": "M92129.1",
    "LC24": "M90953.1",
}


sequences_family = []
for name, acc in dentist_sequences.items():
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        rec = SeqIO.read(handle, "fasta")
        handle.close()
        print(rec.id)
        rec.id = f"{name}_{rec.id}"
        #rec.description = name
        sequences_family.append(rec)
        print(f"Success {name}: {len(rec.seq)} nt")
    except Exception as e:
        print(f"Error {name}: {e}")
    time.sleep(1)

SeqIO.write(sequences_family, "hiv_samples.fasta", "fasta")
