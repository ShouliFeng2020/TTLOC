#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $version_num = "1.0";

my $help_message = Full_help($version_num);

#die if no any input info:
unless($ARGV[0]){
	die "$help_message" ;
}



#######################################DECLARE VARIABLES##########################################

my ($genome,$input1,$input2,$tDNA,$prefix);

my ($help,$version);


my $threads = 1;
my $LB= 500;
my $RB = 500;


############################### OPTIONS #######################################################
#get options:
GetOptions ('genome=s' => \$genome,
			'1=s' => \$input1,
			'2=s' => \$input2,
			'prefix=s' => \$prefix,
			'tDNA=s' => \$tDNA,

			't=i' => \$threads,
			'LB=i' => \$LB,
			'RB=i' => \$RB,
			'version' => \$version,
			'help' => \$help);

# If help, print and quit:
if ($help){
	die "$help_message\n";
}

# If version, print version and quit:
if ($version){
	die "TTLOC version $version_num.\n";
}

unless ($genome or $tDNA){
	die "No genome reference or tDNA reference.\n";
}

unless ($input1 or $input2){
	die "No NGS input reads.\n";
}

unless ($prefix){
	die "No output prefix.\n";
}



# unless(-e "$genome.amb"){
#	`bwa index $genome`;

# }


`fastp -p -i $input1 -I $input2 -o $input1.fastp -O $input2.fastp -h $prefix.fastp.html`;

`flash $input1.fastp $input2.fastp -o $prefix.flash`;

`bwa mem -t $threads -R "\@RG\\tID:id\\tSM:sample\\tLB:lib" $genome $input1.fastp $input2.fastp -o $prefix.bwa.sam`;

`samtools view -h $prefix.bwa.sam | samblaster --addMateTags --maxSplitCount 2 --minNonOverlap 20 -o $prefix.samblaster.sam -d $prefix.discord.sam -s $prefix.split.sam`;


open my $split,"$prefix.split.sam" or die;
open my $splitIN,"$prefix.split.sam" or die;
open my $discordIN,"$prefix.discord.sam" or die;

my $tDNA_L;
my %split;
my %support;


while(<$splitIN>){
		chomp;
		if(/SN:tDNA\tLN:(\d+)/){
				$tDNA_L = $1;
			}
		next if /^@/;
		my @mem = split;
		next if $mem[2] =~ /tDNA/;
		next unless /SA:Z:tDNA/;

### remove PCR duplicates

		my $flag = sprintf("%012b",$mem[1]);

		my @flag = split//,$flag;

		next if $flag[1];

### remove PCR duplicates
		my $ref_st = $flag[7] ? "-":"+";

		my $ref = $mem[2];
		my ($ref_side,$break,$slip);

##  Get insertion site based on split reads
		if($mem[5] =~ /\b(\d+)[SH]\d+M\b/ ){
	 		$ref_side = "Right";
	 		$break = $mem[3];
	 		$slip = $1;
		}elsif($mem[5] =~ /\b(\d+)M(\d+)[HS]\b/){
	 		$ref_side = "Left";
	 		$break = $mem[3] + $1 - 1;
	 		$slip = $2;
		}elsif($mem[5] =~ /\b(\d+)[HS](\d+)M(\d+)[SH]\b/ ){
			if($1 > $3){
			 	$ref_side = "Right";
		 		$break = $mem[3];
		 		$slip = $1;
			}elsif($1 < $3){
		 		$ref_side = "Left";
		 		$break = $mem[3] + $2 -1;
		 		$slip = $3;
			}
		}

### Get insertion site based on split reads
next unless $break;

/SA:Z:([\w,\-\+]+);/;

my @tDNA = split/,/,$1;

### Get tDNA break site
my ($tDNAr,$tbreak,$tside,$type,$tmap);
	if($tDNA[3] =~ /\d+[HS](\d+)M/){
		$tside = "Right";
		$tbreak = $tDNA[1];
		$tmap = $1;
	}elsif($tDNA[3] =~ /(\d+)M\d+[HS]/){
		$tside = "Left";
		$tbreak = $tDNA[1] + $1 - 1;
		$tmap = $1;
	}elsif($tDNA[3] =~ /(\d+)[HS](\d+)M(\d+)[HS]/){
		if($1 > $3){
			$tbreak = $tDNA[1];
			$tside = "Right";
		}elsif($1 < $3){
			$tbreak = $tDNA[1] + $2 -1;
			$tside = "Left";
		}
			$tmap = $2;
	}
### Get tDNA break site

my $direct = $ref_side eq $tside ? "Reverse":"Forward";


my $insert = $slip - $tmap;


## judgment the tDNA repeat
my $RBl = $tDNA_L - $RB;
	if($tbreak < $LB){
		$tDNAr = "LB";
	}elsif($tbreak < $tDNA_L && $tbreak > $RBl){
		$tDNAr = "RB";
	}else{
		$tDNAr = "Outer";
	}


	### retain tDNA break site
	#my $key = "$ref:$break\t$ref_side\t$tDNA\t$tbreak\t$direct";

	## ignore tDNA break site
	my $key = "$ref:$break\t$ref_side\t$tDNAr\t$tbreak\t$direct";

## ignore multipealign reads
	if($flag[3]){
		next;
	}else{
		$split{$key}++;
	}

	if(exists $support{$key}){
		$support{$key} .= ";$mem[0]";
	}else{
		$support{$key} = $mem[0];

	}

}


my %discord;
#my %discord_m;
my %dis_support;
while(<$discordIN>){
	chomp;
	next if /^@/;
	my @mem = split;
	next unless $mem[6] eq "tDNA";
	my $flag = sprintf("%012b",$mem[1]);
	my @flag = split//,$flag;
	next if $flag[1];
	#print "$_\n";
	foreach(sort keys %split){
		my @key = split/\t/,$_;
		my @chr = split/\:/,$key[0];
		if($mem[2] eq $chr[0]){
			if($key[1] eq "Left" && $mem[3] < $chr[1]){
				if($key[4] eq "Forward" && $key[3] < $mem[7]){
					  	$discord{$_}++;
				}elsif($key[4] eq "Reverse" && $key[3] > $mem[7]){
								$discord{$_}++;
				}else{
					next;
				}

			}elsif($key[1] eq "Right" && $mem[3] > $chr[1]){
				if($key[4] eq "Forward" && $key[3] > $mem[7]){
								$discord{$_}++;
				}elsif($key[4] eq "Reverse" && $key[3] < $mem[7]){

							  $discord{$_}++;
					}else{
						next;
					}
				}else{
					next;
				}
				if(exists $dis_support{$_}){
					$dis_support{$_} .= ";$mem[0]";
				}else{
					$dis_support{$_} = "$mem[0]";
				}
			}
		}
}

open my $fq1,"$input1.fastp" or die;
open my $fq2,"$input2.fastp" or die;
open my $fla,"$prefix.flash.extendedFrags.fastq" or die;

my %fq1;
my %fq2;
my %flash;

while(my $fi = <$fq1>){
	my $se = <$fq1>;
	chomp $se;
	<$fq1>;
	<$fq1>;
	chomp $se;
	my @ar = split/\s+/,$fi;
	$ar[0] =~ s/^@//;
	$fq1{$ar[0]} = $se;
}

while(my $fi = <$fq2>){
	my $se = <$fq2>;
	chomp $se;
	<$fq2>;
	<$fq2>;
	chomp $se;
	my @ar = split/\s+/,$fi;
	$ar[0] =~ s/^@//;
	$fq2{$ar[0]} = $se;
}
while(my $fi = <$fla>){
	my $se = <$fla>;
	chomp $se;
	<$fla>;
	<$fla>;
	my @ar = split/\s+/,$fi;
	$ar[0] =~ s/^@//;
	$flash{$ar[0]} = $se;
}

open my $out,">$prefix.tDNA.summary";


####loading reference genome

use Bio::SeqIO;
my %genome;
my $catchseq_obj = Bio::SeqIO -> new(-file=>"$genome",-format=>"fasta");

while(my $seq_obj = $catchseq_obj -> next_seq){
	my $name = $seq_obj -> display_name;
	my $seq = $seq_obj -> seq;
	$genome{$name} = $seq;
}

# print  header
#print $out  "Ref:Breakpoint\tRefSide\ttDNA\ttDNA_break\tDirection\tSplit_UsupportN\tDiscort_UsupportN\tSplit_MsupportN\tDiscord_Msupport\tSample\n"; # tSplit_supportID\tSplit_supportSeq\tFlank_seq\n";
print $out  "Ref:Breakpoint\tRefSide\ttDNA\ttDNA_break\tDirection\t\tSplit_supportN\tDiscort_supportN\tSample\n"; # tSplit_supportID\tSplit_supportSeq\tFlank_seq\n";


my $reads;
my @fi_su;

open my $sup_out,">${prefix}_split_support_reads.fa" or die;

open my $dis_supO,">${prefix}_discord_support_reads.fa" or die;

open my $flank,">${prefix}_flank_seq.fa" or die;
foreach(sort keys %split){
	my $disU = exists $discord{$_}? $discord{$_} : 0;
	#my $disM = exists $discord_m{$_}?$discord{$_} : 0;
	#my $splM = exists $split_m{$_}? $split_m{$_} : 0;

	my @mem = split/\t/,$_;
	my $supN = $split{$_} + $disU;
	my $key = $_;
	if($supN >= 2){
		my @support = split/;/,$support{$_};
		foreach my $su (@support){

		my @reads = split/_/,$su;
		@fi_su = @reads;
			if(exists $flash{$reads[0]}){
				 $reads = $flash{$reads[0]};
			}
			my $spsr = $fi_su[1] == 1 ? $fq1{$fi_su[0]} : $fq2{$fi_su[0]};
			print $sup_out ">$mem[0].$su\n$spsr\n";

		}
		unless($reads){
			$reads = $fi_su[1]== 1? $fq1{$fi_su[0]} : $fq2{$fi_su[0]};
		}

		my @dis_support = split/;/,$dis_support{$_};
		foreach my $diss (@dis_support){
			print $dis_supO ">$_\t$diss\ R1\n$fq1{$diss}\n>$_\t$diss\ R2\n$fq2{$diss}\n";

		}




	my @loci = split/\:/,$mem[0];
	my $star = $loci[1] - 1000;
	my $end = $loci[1] + 1000;
	my $subseq = substr($genome{$loci[0]},$star,2000);
	print $flank ">${mem[0]}_${mem[1]}\n$subseq\n";
	print $out "$mem[0]\t$mem[1]\t$mem[2]\t$mem[3]\t$mem[4]\t$split{$_}\t$disU\t$prefix\n" if $split{$key} > 0; #
	#print $out "$mem[0]\t$mem[1]\t$mem[2]\t$mem[3]\t$mem[4]\t$split{$key}\t$disU\t$splM\t$disM\t$ARGV[3]\n" if $split{$key} > 0; # t$support{$key}\t$reads\t$subseq\n" if $split{$key} > 0;



	}

}

`rm $prefix.flash.*`;


`samtools fasta $prefix.discord.sam -N -1 $prefix.dis_1.fa -2 $prefix.dis_2.fa`;

`samtools fasta $prefix.split.sam -N -s $prefix.split_s.fa`;

`cat $prefix.dis_1.fa $prefix.dis_2.fa $prefix.split_s.fa > $prefix.displit.fa`;

`makeblastdb -in $tDNA -dbtype nucl -out $tDNA`;

`blastn  -db $tDNA -word_size 50 -query $prefix.displit.fa -out $prefix.displitVStDNA.txt -max_target_seqs 10 -evalue 0.00001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore'`;

open my $dis1,"$prefix.dis_1.fa";
open my $dis2,"$prefix.dis_2.fa";
open my $spl,"$prefix.split_s.fa";


my (%dis_1,%dis_2,%split_s);


	&Getfasta($dis1,\%dis_1);
	&Getfasta($dis2,\%dis_2);
	&Getfasta($spl,\%split_s);

open my $blastn,"$prefix.displitVStDNA.txt";
my %uniq;

open my $dis1out,">$prefix.dis_1.extra.fa";
open my $dis2out,">$prefix.dis_2.extra.fa";
open my $splout,">$prefix.split_s.extra.fa";


while(<$blastn>){
chomp;
my @mem = split;

#$mem[0] =~ s/\/\d$//;
next if exists $uniq{$mem[0]};

print $dis1out ">$mem[0]\n$dis_1{$mem[0]}\n" if exists $dis_1{$mem[0]};
print $dis2out ">$mem[0]\n$dis_2{$mem[0]}\n" if exists $dis_2{$mem[0]};
print $splout ">$mem[0]\n$split_s{$mem[0]}\n" if exists $split_s{$mem[0]};
#print $spl2out ">$mem[0]\n$split_2{$mem[0]}\n" if exists $split_2{$mem[0]};

$uniq{$mem[0]}++;

}

`cat $prefix.dis_1.extra.fa $prefix.dis_2.extra.fa $prefix.split_s.extra.fa  > $prefix.displit.extra.fa`;

`makeblastdb -in $genome -dbtype nucl -out $genome`;

`blastn -num_threads $threads -db $genome -query $prefix.displit.extra.fa -out $prefix.displit.extraVSGenome.txt -max_target_seqs 10 -evalue 0.00001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore'`;


`rm $prefix*.fastp`;
`rm $prefix.dis_*.fa`;
`rm $prefix*.sam`;
`rm $prefix.displit.fa`;
`rm $prefix.displitVStDNA.txt`;
`rm $prefix.fastp.html`;
`rm $prefix.split_s*.fa`;
`rm fastp.json`;

sub Getfasta {
	my ($file,$hash) = @_;
	while(my $id = <$file>){
		my $seq = <$file>;
		chomp $id;
		chomp $seq;
		$id =~ s/^>//;
		$$hash{$id} = $seq;
	}
}


################################# HELP INFO #######################################################

sub Full_help {

    my($ver) = @_;
    my $help_info = "\nTTLOC.pl version $ver

########################## TTLOC_v1.0.pl ##########################

##Running TTLOC from the following command:

Usage: perl TTLOC.pl --genome <genome_file.fa> -1 <input_1.fq> -2 <input_2.fq> --prefix <prefix> --tDNA <tDNA.fa> [options]

The followings are the detailed descriptions of the arguments and options in the use of TTLOC:

Arguments:

	--genome <string>. Supply merged genome sequence in FASTA format as reference sequences.

	-1 <string>. Input READ1 file.

	-2 <string>. Input READ2 file.

	--prefix <string>. Prefix of output.

	--tDNA <string>. Supply tDNA sequence in FASTA format as reference.

	-t <int>. Number of  threads. The default is 1.

	--LB <int>. Length of tDNA Left Repeat. The default is 500 nt.

	--RB <int>. Length of tDNA Right Repeat. The default is 500 nt.

	--help. Print the help message and quit.

	--version. Print TTLOC version number and quit.

	Type \'TTLOC.pl --help\' for full list of options

";
    return $help_info;
}
