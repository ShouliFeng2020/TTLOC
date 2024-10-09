
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


## Step 2 - Built index for Merged genome

`bwa index  genome.fa tDNA.fa`

## Step 3 - Identify the T-DNA integration sites
  
## Usage: 

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


## Output:
### summary file
TTLOC produces a single summary file which contains all candidate tDNA integration sites that were identified.

The output is "prefix".tDNA.summary

#### summary file structure

* Ref:Breakpoint: 	The tDNA integration Chromosome number and integration site;	
* RefSide: 		The side of tDNA integration genome that were identified;
* tDNA: 		The tDNA repeat region (Left or Right repeat) integrated into the genome;
* tDNA_break: 		The site of tDNA break were identified.
* Direction: 		The direction of tDNA integration into genome;
* Split_supportN: 	The split support reads number;
* Discort_supportN: 	The discort support reads number.
* Sample: 		The sampleID, determined by parameter "prefix".

### support reads and flank sequence in fasta format
The TTLOC provides support reads based on BWA alignment,  which are saved in files named "prefix"_split_support_reads.fa and "prefix"_discort_support_reads.fa, respectively.

The TTLOC also provides support reads for BLASTN alignment, which are saved in files named "prefix".displit.extra.fa.

The 2000bp sequences flank the insertion site is saved in the "prefix"_flank_seq.fa.

### BLASTN alignment
The split and discort support reads were aligned to plant genome using blastn and the alignment results are saved in "prefix".displit.extraVSGenome.txt file in blast "outfmt 6" format

 	1.  qseqid      query or source (gene) sequence id

  	2.  sseqid      subject or target (reference genome) sequence id

   	3.  pident      percentage of identical positions

   	4.  length      alignment length (sequence overlap)

   	5.  mismatch    number of mismatches

   	6.  gapopen     number of gap openings

   	7.  qstart      start of alignment in query

   	8.  qend        end of alignment in query

   	9.  sstart      start of alignment in subject

 	10.  send        end of alignment in subject

 	11.  evalue      expect value

	12.  bitscore    bit score



# Contact

* Dr. Shouli Feng    (fengshouli@xhlab.ac.cn)

