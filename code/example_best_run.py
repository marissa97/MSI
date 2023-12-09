from msi import *
import numpy as np
##### example ####
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

outDir= "/home/marissa/Documents/sysmex/MSI/results/example_best/"
## create dummy test fasta file
# shortage 5 "A"
best_bat26="cccttaacctttttcaggtaaaaaaaaaaaaaaaaaaaagggttaaaaatgttgaatggttaaaaaatgtt"
chunks1 = [best_bat26[i:i+20] for i in range(0, len(best_bat26), 20)]
chunks2 = [best_bat26[i:i+20] for i in range(2, len(best_bat26), 20)]
chunks3 = [best_bat26[i:i+20] for i in range(4, len(best_bat26), 20)]
chunks4 = [best_bat26[i:i+20] for i in range(6, len(best_bat26), 20)]

chunks= np.concatenate((chunks1, chunks2,chunks3,chunks4), axis=None)


sequences=[]


for i in range(0,len(chunks)):
    record = SeqRecord(
        Seq(chunks[i]),
        id=str(i),
        name="bat",
        description="test"
    )
    sequences.append(record)  # add code here

SeqIO.write(sequences, outDir+"example.fasta", "fasta")

## ref for bat26 fasta
bat26_temp="cccttaacctttttcaggtaaaaaaaaaaaaaaaaaaaaaaaaaaagggttaaaaatgttgaatggttaaaaaatgtt"
record = SeqRecord(
        Seq(bat26_temp),
        id=str(i),
        name="referece",
        description="reference"
    )
SeqIO.write(record, outDir+"bat26.fasta", "fasta")

### parameter
input_fastq_file = outDir+"example.fasta"
bat26_temp="cccttaacctttttcaggtaaaaaaaaaaaaaaaaaaaaaaaaaaagggttaaaaatgttgaatggttaaaaaatgtt"
bat26_start= 0
bat26_end=78
ad_start= bat26_start+20
ad_end= bat26_start+ 46
read_length=20


## mapping 
bamfile= mapFastatoBam("example", "bat26.fasta", outDir)

## depth and cov calc for adenosin 
dep, cov= bamReads(bamfile, ad_start, ad_end, len(bat26_temp),read_length, outDir)

plotBar( ad_start, ad_end, dep, "Read Depths", "Adenosine Shortage Metrics (Depths)", "dep",outDir)
plotBar( ad_start, ad_end, cov, "Coverage", "Adenosine Shortage Metrics (Coverage)", "cov",outDir)
##
