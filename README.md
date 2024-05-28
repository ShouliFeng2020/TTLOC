# Citation


# Installing TTLOC 

## Software Dependencies
TTLOC has been tested on the Ubuntu 11.4.0 Linux:
The following programs need to be installed and the executable commands should be in $PATH of system.
* fastp (Version = "0.23.4")
* flash (Version = "1.2.11")
* BWA (Version ="0.7.17")
* Samtools (Version ="1.19")
* perl (Version = "v5.32.1")
* blast+ (Version = "2.15.0")

# Using TTLOC

## Step 1 - Merge genome reference and tDNA reference
The name of tDNA reference sequence shoud be: "tDNA"

`cat genome.fa tDAN.fa > genome.tDNA.fa`

### Example:

An example data set is provided with this repository.

Running the following example code will create a merged reference:

`cat `

## Step 2 - Identify the T-DNA integration sites
  
### Usage: 

`perl TTLOC.pl -g genome.tDNA.fa -1 forward.fq -2 reverse.fq -t tDNA.fa -s sample`

### Parameters:

* REQUIRED -1 the paired read file 1
* REQUIRED -2 the paired read file 2
* REQUIRED -t the T-DNA sequence file in fasta format
* REQUIRED -g the genome sequence file in fasta format
* REQUIRED -p the name of your project (output files will be placed in a directory with the name you provide)
* OPTIONAL:
	* -@ cpu number for BWA and SAMTOOLS [default 8]
	* -a the window size of clustering soft clipped reads [default:3]
	* -b the length of library fragment in NGS data [default:500]

### Example:

An example data set is provided with this repository.

Running the following example code will create a project directory named 'tdna' relative to where you run the command, and will produce example output:

`python tdnascan.py -1 mt4_chr1_20x_mut_tdna_1.fq -2 mt4_chr1_20x_mut_tdna_2.fq -t t-dna_elison.fa -g mt4_chr1_2Mb.fa -p tdna`

## Step 2 - Annotate complete and truncated T-DNA insertions

### Usage: 

`python tdnaAnnot.py -i tdna_insertion.bed -f ref.gff3 -o tdna_insertion_annot.bed`

### Parameters:

* REQUIRED -i T-DNA BED file
* REQUIRED -f gff3 annotation file
* REQUIRED -o annotated insertion file

### Example:

`python tdnaAnnot.py -i tdna/5.tdna_insertion.bed -f Athaliana_447_Araport11.gene.gff3 -o ./tdna/5.tdna_insertion_annot.bed`


### Output:

### BED file
TDNAscan produces a single BED file which contains all unique deletions that were identified.

The output is placed in ./tdna (i.e. in the directory named after your project)

Running the above example code (Step 1) would produce the following BED file:

* ./tdna/**5.tdna_insertion.bed**

Annotated BED file (Step 2):

* ./tdna/**5.tdna_insertion_annot.bed**

### Output file structure

* Chr: Chromosome number;	
* Position: Start position of insertions (~ represents insertion position nearby);
* SuppRead: CLR represents the clipped reads number; DIR represents discordant reads number;
* SuppRead: tdna_st and tdna_end represent the start and end position of T-DNA sequence truncated when inserted to reference genome.
* Orientation: forward or reverse T-DNA inserted to reference genome;
* Freq: Insertion frequency;
* Genes (optional): This column will only show genes if deletions cover.


# Contact

* Dr. Liang Sun    (sunliang@udel.edu)
* Yinbing  Ge  (yinge@noble.org)
