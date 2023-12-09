from Bio import SeqIO
import subprocess
import pysam
import numpy as np
import matplotlib.pyplot as plt 
## 1. read fastq (QC)

## 2. map fastq to reference genome -> bam/sam file  or map to bat26_temp
## illumina limitation: 
## - reads are expected to normally have the same size, but in reallity, this can shift
## - read length of NGS platforms is short (50 to 300), each position may have different coverage  
## - a read might not cover all the adenine region (depends on the length as well), 

##TODO: not specific for bat26, (this is for reads < 70bp)
def mapFastatoBam(name, refname, outDir):
    # Index the reference genome
    #gatk-launch CreateSequenceDictionary -R [genomeRefFasta] -> for vcf
    subprocess.run(["samtools", "faidx", outDir+refname])
    subprocess.run(["bwa", "index", outDir+refname])
    # Map reads to the reference genome
    subprocess.run(["bwa", "aln",outDir+refname, outDir+name+".fasta"], stdout=open(outDir+name+".sai", "w"))
    subprocess.run(["bwa", "samse",outDir+refname, outDir+name+".sai", outDir+name+".fasta"], stdout=open(outDir+"mapped_reads.sam", "w"))
    # Convert SAM to BAM
    subprocess.run(["samtools", "view", "-bS", outDir+"mapped_reads.sam", "-o", outDir+"mapped_reads.bam"])
    # Sort the BAM file
    subprocess.run(["samtools", "sort", outDir+"mapped_reads.bam", "-o", outDir+name+".bam"])
    # Index the sorted BAM file
    subprocess.run(["samtools", "index", outDir+name+".bam"])
    return outDir+name+".bam"

##3. select region
## if map to human genome -> extracting the sam file to bat26 region 
#(chr: bam_file.get_reference_name(read.reference_id), start: read.reference_start, end)




##4. Calculate depth and coverage to prove MSI 
## depth can be calculated with samtools: samtools depth [options] [in1.sam|in1.bam|in1.cram [in2.sam|in2.bam|in2.cram] [...]] 
## limitation: if not 'a' (mismatch)?
##TODO: usage of cigar string, if reads may be splitted (?)
## if not: 
def bamReads(bamfile, ad_st, ad_end, reflength,read_length, outDir):
    len_ad= ad_end-ad_st
    dep_arr=np.zeros(len_ad)
    # Open the BAM file
    with pysam.AlignmentFile(bamfile, 'rb') as bf:
        # Iterate over reads in the BAM file
        for read in bf.fetch():
            read_name = read.query_name
            sequence = read.query_sequence
            positions= read.positions
            # illumina limitation: if reads < adenosine length (e.g. aaaa, this can be mapped to any position in the adenosin repeats)
            dep_arr= dep_arr+ calcReads(len_ad,ad_st, ad_end, positions, sequence)
    cov_arr= calcCoverage(dep_arr, read_length, reflength) 
    return dep_arr,cov_arr
                

#calculate depth
def calcReads(len_ad,ad_st, ad_end, positions, sequence):
    dep_arr=np.zeros(len_ad)
    positions= np.array(positions)
    positions= positions[np.where((positions >= ad_st) & (positions <= ad_end))]
    sequence=  np.array(list(sequence))
    sequence= sequence[np.where((positions >= ad_st) & (positions <= ad_end))]
    if len(positions)!=0:
        for i in range(0, len(positions)):
            if sequence[i]=="A":
                dep_arr[positions[i]-ad_st]=dep_arr[positions[i]-ad_st]+1
    return dep_arr


## calculating coverage: per position: (number of reads) x (length of each read) / (length of the reference)
## e.g. if in the position coverage <10% -> shortened
def calcCoverage(dep_arr, read_length, reflength):
    cov_arr=np.zeros(len(dep_arr))
    for i in range(0,len(dep_arr)):
        cov_arr[i]= dep_arr[i] * read_length / reflength
    return cov_arr



## 5. Decision making, definition of adenosin region that is shortened
## Metrics 
#dep_arr, cov_arr=bamReads(bamfile,seq, ad_st, ad_end, reflength,read_length)

#- how does the distribution of the depth looked like?
#- mean depth of all region? mean depth of the adenosin region?
# e.g. depth: decision: region hast to have mind 5 reads
## e.g. coverage: decision:if in the position coverage <10% -> shortened
def plotBar( ad_st, ad_end, yval, ylab, title,name,outDir):
    xval= list(range(ad_st, ad_end))
    yval = list(yval)
    fig = plt.figure(figsize = (10, 5))
    # creating the bar plot
    plt.bar(xval, yval, color ='maroon', 
        width = 0.4)
    plt.xlabel("Adenosin Positions")
    plt.ylabel(ylab)
    plt.title(title)
    #plt.show()
    plt.savefig(outDir+name)



### illumina limitations: noise




