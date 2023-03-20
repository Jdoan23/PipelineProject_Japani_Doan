# Project Name: PipelineProject_Japani_Doan
# Project Purpose:
   The purpose of this project is to to develop a Python wrapper to automate the execution of specfic software tools we have discussed in COMP 383 to analyze HCMV with donor 1 and 3. This project will be done completely the command line.
   HCMV is linked here: https://www.ncbi.nlm.nih.gov/pubmed/29158406
   
# NOTE for Packages Needed
Essential items to run this project and code:

Bowtie2 
Site Recommended: https://www.metagenomics.wiki/tools/bowtie2/install
This code will be used within the command line.
```Python 
conda install -c bioconda bowtie2
```
SPAdes 
Sites Recommended: https://cab.spbu.ru/files/release3.12.0/manual.html
This code will be used within the command line.
```Python 
wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz

tar -xzf SPAdes-3.12.0-Linux.tar.gz

cd SPAdes-3.12.0-Linux/bin/
```
BLAST+ 
Site Recommended: https://www.ncbi.nlm.nih.gov/books/NBK279690/
```Python 
Install:

rpm -ivh ncbi-blast-2.2.18-1.x86_64.rpm
```
SRA-Toolkit : https://hpc.nih.gov/apps/sratoolkit.html 

Entrez Direct 
Sites Recommended: https://www.ncbi.nlm.nih.gov/books/NBK179288/
```Python
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

OR

sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```
# Other Packages Needed: 
import os

import Bio

from Bio import Entrez

from Bio import SeqIO

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
 These are the files needed. Using Bowtie make .sam files. Use -1 and -2 as a way to paired-end reads. From this point, SPAdes will be needed to pair the ends of the file. In this case it wil also create a new file. 
```Python
#create an index for HCMV
os.system("bowtie2-build Genome.fasta HCMV_Ind")

#building .sam files using bowtie2 
os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660030_1.fastq  -2 SRR5660030_2.fastq -S 30Genomemap.sam")
os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660033_1.fastq  -2 SRR5660033_2.fastq -S 33Genomemap.sam")
os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660044_1.fastq  -2 SRR5660044_2.fastq -S 44Genomemap.sam")
os.system("bowtie2 --quiet -x HCMV_Ind -1 SRR5660045_1.fastq  -2 SRR5660045_2.fastq -S 45Genomemap.sam")

#SPAdes assembly to pair ends reads of the file 
# --al-conc-gz will allow alignment in the file
os.system('bowtie2 -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S 30Genomemap.sam --al-conc-gz 30Genome_MAP_%.fq.gz')
os.system('bowtie2 -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S 33Genomemap.sam --al-conc-gz 33Genome_MAP_%.fq.gz')
os.system('bowtie2 -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S 44Genomemap.sam --al-conc-gz 44Genome_MAP_%.fq.gz')
os.system('bowtie2 -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S 45Genomemap.sam --al-conc-gz 45Genome_MAP_%.fq.gz')
```
# Step 3
This will assemble all four transcriptomes together using SPAdes. 
```Python
os.system('spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 30Genome_MAP_1.fq.gz --pe-2 1 30Genome_MAP_2.fq.gz --pe-1 2 33Genome_MAP_1.fq.gz --pe-2 2 33Genome_MAP_2.fq.gz --pe-1 3 44Genome_MAP_1.fq.gz --pe-2 3 44Genome_MAP_2.fq.gz --pe-1 4 45Genome_MAP_1.fq.gz --pe-2 4 44Genome_MAP_2.fq.gz -o HCMV2-SRR_assembly') 
```
# Step 4
Python code will be needed to count the reads for transcriptomes for each file and then place it into the Longest.txt file. 
```Python
#Longest contig 
#How it works: 
  #Each seq will be appended and recorded as the longest contig
from Bio import SeqIO

Contig = []
with open('HCMV2-SRR_assembly/contigs.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
        Contig.append(record.seq)
        Longest = (max(Contig, key = len)) 
          #Function uses the length of each contig to determine the maximum value

#It will also write the longest into an output file. 
with open('Longest.txt', 'w') as file: 
    file.write(str(Longest)) 
    file.close()
 ``` 
 Now, we need to find the length of the assembly. This is a needed step to check if the length of will be smaller than 1000. 
 ```Python
 #Length of the assembly
#How it works
  #Using the same as before, open into the contigs.fasta file 
  #While recording in the list named "LenOfAsb", it will appened each seq 
LenOfAsb = []
with open('HCMV2-SRR_assembly/contigs.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
        LenOfAsb.append(len(record.seq))
        Start = 0
        #For statement that if the length is smaller then 1000, then add to count by increments of 1. After that it will continue the list. 
        for StartingLength in LenOfAsb:
            if StartingLength > 1000:
                Start += 1
                Largest = [StartingLength]
#Contig Count 
#How it works 
  #Open file and see what is the number of contings. Then write into file 
with open('Counts.txt', 'w') as outfile:
    Number = (sum(Largest))
    Start  = (str(Start))
    outfile.write(Start + ' > 1000 bp in the assembly')
    print('Basepair Counts: ', sum(Largest))
    print('Contig Counts: ', Start)
    #This will print as needed for the log file. 
```
# Step 5 
For the last step of the project, a local database of just sequences from the Betaherpesvirinae subfamily will be needed. For this to happen, we need to call in the family using specfic tools. 
```Python
Query = '"Betaherpesvirinae"[Organism] OR Betaherpesvirinae[All Fields]'
Command = f'esearch -db Nucleotide -query {Query} | efetch -format fasta > BetaherpesvirinaeData.txt'
#The command will use esearch to read the query. The -db will specfiy the datatbase. Query will connect the search with esearch.
#efetch will also work with esearch. The format will then format the file into an output file, which will be BetaherpesvirinaeData.txt

os.system(Command)

InFile = 'BetData.txt'
OutFile = 'Betaherpesvirinae'
database_type = 'nucl'
title = 'Betaherpesvirinae'

os.system(f'makeblastdb -in {input_file} -out {output_name} -title {title} -dbtype {database_type}')
# In will be for the input file. Out will be for the output file. Title will be the name of the database while dbtyep of fpr the type that was used. 

#Blasting
#Uses longest contig
inputFile = 'Longest.txt'
outputFile = 'BetResults.csv'
Blast = 'blastn -query ' + inputFile + ' -db Betaherpesvirinae -out ' + outputFile + ' -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitles" > blastout.tsv '
print(Blast)

os.system(Blast)

#making table of blast command outputs
os.system('echo "sacc pident length qstart qend sstart send bitscore evalue stitle" | cat - BetResults.csv > BetResults.tsv')
os.system('head -n 11 BetResults.tsv  > BetResults2.tsv')
```
