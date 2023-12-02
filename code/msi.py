from Bio import SeqIO
import subprocess
import pysam

input_fastq_file = "your_illumina_data.fastq"
bat26_temp="cccttaacctttttcaggtaaaaaaaaaaaaaaaaaaaaaaaaaaagggttaaaaatgttgaatggttaaaaaatgtt"
bat26_start= 100
bat26_end=178
ad_start= bat26_start+19
ad_end= bat26_start+ 46
## 1. read fastq

## 2. map fastq to reference genome -> bam/sam file  or map to bat26_temp
## illumina limitation: 
## - reads are expected to normally have the same size, but in reallity, this can shift
## - read length of NGS platforms is short (50 to 300), each position may have different coverage  
## - a read might not cover all the adenine region (depends on the length as well), 

# Index the reference genome
#subprocess.run(["samtools", "faidx", "reference_genome.fasta"])

# Map reads to the reference genome
#subprocess.run(["bwa", "mem", "reference_genome.fasta", "reads.fastq"], stdout=open("mapped_reads.sam", "w"))

# Convert SAM to BAM
#subprocess.run(["samtools", "view", "-bS", "mapped_reads.sam", "-o", "mapped_reads.bam"])

# Sort the BAM file
#subprocess.run(["samtools", "sort", "mapped_reads.bam", "-o", "sorted_mapped_reads.bam"])

# Index the sorted BAM file
#subprocess.run(["samtools", "index", "sorted_mapped_reads.bam"])

##3. select region
## if map to human genome -> extracting the sam file to bat26 region 
#(chr: bam_file.get_reference_name(read.reference_id), start: read.reference_start, end)

bat26_bam = 'bat26.bam'

##4. Calculate depth to prove MSI 
## depth can be calculated with samtools: samtools depth [options] [in1.sam|in1.bam|in1.cram [in2.sam|in2.bam|in2.cram] [...]] 
## limitation: if not 'a' (mismatch)?

## if not: 
def bamReads(bamfile,seq, ad_st, ad_end, reflength,read_length):
    len_ad= ad_end-ad_st
    dep_arr=np.zeros(len_ad)
    # Open the BAM file
    with pysam.AlignmentFile(bamfile, 'rb') as bamfile:
        # Iterate over reads in the BAM file
        for read in bamfile.fetch():
            # Get the read name and sequence
            read_name = read.query_name
            sequence = read.query_sequence
            start_position = read.reference_start
            end_position = read.reference_end

            ## if shortening in adenine -> msi detection ?
            ## using covarage: 
            # illumina limitation: if reads < adenosine length (e.g. aaaa, this can be mapped to any position in the adenosin repeats)
            
            dep_arr= dep_arr+ calcReads(len_ad,ad_st, ad_end, start_position,end_position, sequence)
    cov_arr= calcCoverage(dep_arr, read_length, reflength)

    return dep_arr,cov_arr
                

    


def calcReads(len_ad,ad_st, ad_end, start_position,end_position, sequence):
    dep_arr=np.zeros(len_ad)
    if ad_st<=start_position<=ad_end and end_position<=ad_end:
        ad= sequence[ad_st:ad_end]
        for i in range(0,len(ad)):
            if ad[i]=="a":
                dep_arr[i]=cov_arr[i]+1
    return dep_arr


## Another approach is calculating coverage: per position: (number of reads) x (length of each read) / (length of the reference)
## e.g. if in the position coverage <10% -> shortened

def calcCoverage(dep_arr, read_length, reflength):
    cov_arr=np.zeros(len(dep_arr))
    for i in range(0,len(dep_arr)):
        cov_arr[i]= dep_arr[i] * read_length/reflength


## 5. Decision making, definition of adenosin region that is shortened
## Metrics 
#dep_arr, cov_arr=bamReads(bamfile,seq, ad_st, ad_end, reflength,read_length)

#- how does the distribution of the depth looked like?
#- mean depth of all region? mean depth of the adenosin region?
# e.g. depth: decision: region hast to have mind 5 reads



## e.g. coverage: decision:if in the position coverage <10% -> shortened

### illumina limitations: noise