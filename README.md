# Project Name: PipelineProject_Japani_Doan
# Project Purpose:
   The purpose of this project is to to develop a Python wrapper to automate the execution of specfic software tools we have discussed in COMP 383 to analyze HCMV with donor 1 and 3. This project will be done completely the command line.
   HCMV is linked here: https://www.ncbi.nlm.nih.gov/pubmed/29158406
# Step 1
   First, "!pip install Bio" needs to be downloaded first. This will download BioPython, which will be needed in coding other parts of this pipeline 
   Next, we need to import os, import Bio, and from Bio import Entrez. Import os will run python as script. Import Bio will be used for addtional assistance in  biological computation. Finally, Bio import Entrez will assist in online search into the NCBI site.
   To start the code, the user will need to obtain the AWS links by using the links.
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363 
Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374 
Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
```Python
#Problem #1
import os 
import Bio 
from Bio import Entrez
import logging

#Make directory and change directory
#os.system("mkdir PipelineProject_Japani_Doan")
#os.chdir("PipelineProject_Japani_Doan")
```

# Step 2 
   Once this is done, using "wget" on each link will retreive and download the file URL into the command line.
```Python
#Use wget to download URLs
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044")
os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045")
```

The "fastq-dump -I --split-files" is then used with each file obtained about using "wget" to split the file in half. This will create file names ending in 1 and 2.
```Python
os.system ("fastq-dump -I --split-files SRR5660030")
os.system ("fastq-dump -I --split-files SRR5660033")
os.system ("fastq-dump -I --split-files SRR5660044")
os.system ("fastq-dump -I --split-files SRR5660045")
```

To reads map to the HCMV genome the Bio import Entrez will be used retrive the HCMV file. This will record and then close the file. It will be stored in the command line.
```Python
#Reads map to the HCMV genome
Entrez.email = "japani.doan2@gmail.com"
handle = Entrez.efetch(db="nucleotide", id=["NC_006273.2"], rettype="fasta")
records = handle.read() 
handle.close()
```
Then, the index of the DNA sequence will be opened with python code to prepare it for BowTie 
```Python
#Building index from set of DNA sequences
with open('Genome.fasta', 'w') as file:
  file.write(records)
```

   BowTie is a tool meant to work through The Burrows-Wheeler matrix. Bowtie works by using "seed" substrings to see any specfic matches within the genome. It is then aligned and ensure that there are no gaps. In this case we will be using BowTie2. "os.system("bowtie2-build fasta HCMV")" will be used to download BowTie2 to build the HCMV fasta. 
  Sequence Alignment/Map will be create in which it will so multiple sequence alignments within a genome. For more information, use this link:  http://samtools.github.io/hts-specs/SAMv1.pdf
 These are the files needed. Using Bowtie make .sam files. Use -1 and -2 as a way to paired-end reads.
```Python
os.system("bowtie2-build Genome.fasta HCMV_Ind")

os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660030_1.fastq  -2 SRR5660030_2.fastq -S 30Genomemap.sam")
os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660033_1.fastq  -2 SRR5660033_2.fastq -S 33Genomemap.sam")
os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660044_1.fastq  -2 SRR5660044_2.fastq -S 44Genomemap.sam")
os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660045_1.fastq  -2 SRR5660045_2.fastq -S 45Genomemap.sam")
```
Python code will be needed to count the reads for transcriptomes for each file. 
```Python
#Will count the reads for transcriptomes
def count(filename):
  with open(filename, 'r') as f:
    count = 0
    for line in f:
      if not line.startwith('@'):
        count += 1
    return count
 ``` 
 # Step 3
 Bowtie2 is then needed to assemble all four transcriptomes together. For this to be possible, 1 assembly via SPAdes will be done. For more information about SPAdes it will be linked here: http://cab.spbu.ru/software/spades/
 Write a SPAdes command to make a log file. 
 ```Python 
 import os

samples = {
    'SRR5660030': ('SRR5660030_1.fastq', 'SRR5660030_2.fastq'),
    'SRR5660033': ('SRR5660033_1.fastq', 'SRR5660033_2.fastq'),
    'SRR5660044': ('SRR5660044_1.fastq', 'SRR5660044_2.fastq'),
    'SRR5660045': ('SRR5660045_1.fastq', 'SRR5660045_2.fastq'), 
    }
```

Using Python code, the calculation of the number of contigs with a length > 1000 will be needed. The quality will need to be put into the log file. Furthermore, write Python code to also calculate the length of the assembly and place the value within the log file. 
```Python
from Bio import SeqIO
#Count the number of contigs with length > 1000
num_long_contigs = 0
for record in SeqIO.parse(assembly_file, "fasta"):
     #is length of the record is greater then 1000 then the value of the num long contig will add 1 
    if len(record) > 1000:
        num_long_contigs += 1

# Write the result to the log file. This will then contain the contig 
with open(log_file, "a") as log:
   #is lengths of the log file is greater then 1000 then it will put onto new line 
    log.write(f"\nNumber of contigs > 1000 bp: {num_long_contigs}\n")

    assembly_file = "Contig.fasta"
log_file = "log.txt"

# Will calculate the length of the assembly
assembly_length = 0
for record in SeqIO.parse(assembly_file, "fasta"):
   #if the length of record is greater than 1000 then add to record.
    if len(record) > 1000:
        assembly_length += len(record)

# Write the result to the log file
with open(log_file, "a") as log:
#is lengths of the log file is greater then 1000 then it will put onto new line 
    log.write(f"\nLength of assembly > 1000 bp: {assembly_length}\n")

    assembly_file = "Contigs.fasta"
```
 Finally, write Python code to retrieve the longest contig from the assembly. Specifically using the longest contig as blast+, query the nr nucleotide database to only the Betaherpesvirinae subfamily. This run should show the best alignment. Only the top 10 hits will be needed, which the information that will be pulled is listed below: 
Subject accession, Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title.
 The headers for each will be shorted to, "sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"
