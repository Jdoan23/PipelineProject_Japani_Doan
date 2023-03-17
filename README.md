# Project Name: PipelineProject_Japani_Doan
# Project Purpose:
   The purpose of this project is to to develop a Python wrapper to automate the execution of specfic software tools we have discussed in COMP 383 to analyze HCMV with donor 1 and 3. This project will be done completely the command line.
   HCMV is linked here: https://www.ncbi.nlm.nih.gov/pubmed/29158406
# How this was done:
   First, "!pip install Bio" needs to be downloaded first. This will download BioPython, which will be needed in coding other parts of this pipeline 
   Next, we need to import os, import Bio, and from Bio import Entrez. Import os will run python as script. Import Bio will be used for addtional assistance in  biological computation. Finally, Bio import Entrez will assist in online search into the NCBI site.
   To start the code, the user will need to obtain the AWS links by using the links.
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363 
Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374 
Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
   Once this is done, using "wget" on each link will retreive and download the file URL into the command line.
   The "fastq-dump -I --split-files" is then used with each file obtained about using "wget" to split the file in half. This will create file names ending in 1 and 2.
   To reads map to the HCMV genome the Bio import Entrez will be used retrive the HCMV file. This will record and then close the file. It will be stored in the command line. 
   Then, the index of the DNA sequence will be opened wirh python code to prepare it for BowTie 
   BowTie is a tool meant to work through The Burrows-Wheeler matrix. Bowtie works by using "seed" substrings to see any specfic matches within the genome. It is then aligned and ensure that there are no gaps. In this case we will be using BowTie2. "os.system("bowtie2-build fasta HCMV")" will be used to download BowTie2 to build the HCMV fasta. 
  Sequence Alignment/Map will be create in which it will so multiple sequence alignments within a genome. For more information, use this link:  http://samtools.github.io/hts-specs/SAMv1.pdf
 These are the files needed. Using Bowtie make .sam files. Use -1 and -2 as a way to paired-end reads.
 Python code will be needed to count the reads for transcriptomes for each file. 
 More python code will then be needed to compare Donor 1 and 3 to HCMV to see the possible similaries. 
 Bowtie2 is then needed to assemble all four transcriptomes together. For this to be possible, 1 assembly via SPAdes will be done. For more information about SPAdes it will be linked here: http://cab.spbu.ru/software/spades/
 Write a SPAdes command to make a log file. 
 Using Python code, the calculation of the number of contigs with a length > 1000 will be needed. The quality will need to be put into the log file. Furthermore, write Python code to also calculate the length of the assembly and place the value within the log file. 
 Finally, write Python code to retrieve the longest contig from the assembly. Specifically using the longest contig as blast+, query the nr nucleotide database to only the Betaherpesvirinae subfamily. This run should show the best alignment. Only the top 10 hits will be needed, which the information that will be pulled is listed below: 
Subject accession, Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title.
 The headers for each will be shorted to, "sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"
