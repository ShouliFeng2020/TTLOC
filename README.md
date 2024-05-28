
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

`cat genome.fa tDNA.fa > genome.tDNA.fa`


## Step 2 - Identify the T-DNA integration sites
  
### Usage: 

`perl TTLOC.pl --genome /path/to/genome.tDNA.fa -1 input_1.fq -2 input_2.fq --prefix sampleID --tDNA /path/to/tDNA.fa `

### Parameters:

* REQUIRED --genome the merged genome and tDNA sequence file in fasta format
* REQUIRED -1 the paired read file 1
* REQUIRED -2 the paired read file 2
* REQUIRED --prefix the prefix and sampleID
* REQUIRED --tDNA the tDNA sequence file in fasta format
* OPTIONAL:
	* -t cpu threads for BWA [default 1]
	* --LB the length of tDNA Left Repeat [default:500]
	* --RB the length of tDNA Right Repeat [default:500]


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
