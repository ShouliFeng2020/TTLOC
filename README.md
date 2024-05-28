
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

### summary file
TTLOC produces a single summary file which contains all candidate tDNA integration sites that were identified.

The output is "prefix".tDNA.summary

### summary file structure

* Ref:Breakpoint: The tDNA integration Chromosome number and integration site;	
* RefSide: The side of tDNA integration genome that were identified;
* tDNA: The tDNA repeat region (Left or Right repeat) integrated into the genome;
* tDNA_break: The site of tDNA break were identified.
* Direction: The direction of tDNA integration into genome;
* Split_supportN: The split support reads number;
* Discort_support: The discort support reads number.
* Sample: The sampleID, determined by parameter "prefix".



# Contact

* Dr. Shouli Feng    (fengshouli@xhlab.ac.cn)

