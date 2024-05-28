#!/usr/bin/perl
use warnings;
use strict;
unless(@ARGV == 5){
	print "\n\tUsage:perl $0 <genome.ref> <q1> <q2> <prefix> <tDNA.ref>\n\n";
	exit;

}

unless(-e "$ARGV[0].amb"){
	`bwa index $ARGV[0]`;

}

`fastp -p -i $ARGV[1] -I $ARGV[2] -o $ARGV[1].fastp -O $ARGV[2].fastp -h $ARGV[3].fastp.html`;

`flash $ARGV[1].fastp $ARGV[2].fastp -o $ARGV[3].flash`;

`bwa mem -R "\@RG\\tID:id\\tSM:sample\\tLB:lib" $ARGV[0] $ARGV[1].fastp $ARGV[2].fastp -o $ARGV[3].bwa.sam`;

`samtools view -h $ARGV[3].bwa.sam | samblaster --addMateTags --maxSplitCount 2 --minNonOverlap 20 -o $ARGV[3].samblaster.sam -d $ARGV[3].discord.sam -s $ARGV[3].split.sam`;


open my $split,"$ARGV[3].split.sam" or die;


open my $splitIN,"$ARGV[3].split.sam" or die;
open my $discordIN,"$ARGV[3].discord.sam" or die;

my $tDNA_L;
my %split;
my %support;
#my %split_m;

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
my ($tDNA,$tbreak,$tside,$type,$tmap);
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
#my $flip = $ref_st eq $tDNA[2] ? "unflip":"flip";

my $insert = $slip - $tmap;


## Hard judgment the tDNA repeat
my $RB = $tDNA_L - 500;
	if($tbreak < 500){
		$tDNA = "LB";
	}elsif($tbreak < $tDNA_L && $tbreak > $RB){
		$tDNA = "RB";
	}else{
		$tDNA = "Outer";
	}


	### retain tDNA break site
	#my $key = "$ref:$break\t$ref_side\t$tDNA\t$tbreak\t$direct";

	## ignore tDNA break site
	my $key = "$ref:$break\t$ref_side\t$tDNA\t$tbreak\t$direct";


#	if($flag[3]){
#		$split_m{$key}++;
#	}else{
#		$split{$key}++;
#	}
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

open my $fq1,"$ARGV[1].fastp" or die;
open my $fq2,"$ARGV[2].fastp" or die;
open my $fla,"$ARGV[3].flash.extendedFrags.fastq" or die;

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

open my $out,">$ARGV[3].tDNA.summary";


####loading reference genome

use Bio::SeqIO;
my %genome;
my $catchseq_obj = Bio::SeqIO -> new(-file=>"$ARGV[0]",-format=>"fasta");

while(my $seq_obj = $catchseq_obj -> next_seq){
	my $name = $seq_obj -> display_name;
	my $seq = $seq_obj -> seq;
	$genome{$name} = $seq;
}

# print  header
#print $out  "Ref:Breakpoint\tRefSide\ttDNA\ttDNA_break\tDirection\tSplit_UsupportN\tDiscort_UsupportN\tSplit_MsupportN\tDiscord_Msupport\tSample\n"; # tSplit_supportID\tSplit_supportSeq\tFlank_seq\n";
print $out  "Ref:Breakpoint\tRefSide\ttDNA\tDirection\tSample\n"; # tSplit_supportID\tSplit_supportSeq\tFlank_seq\n";


my $reads;
my @fi_su;

open my $sup_out,">${ARGV[3]}_split_support_reads.fa" or die;

open my $dis_supO,">${ARGV[3]}_discord_support_reads.fa" or die;

open my $flank,">${ARGV[3]}_flank_seq.fa" or die;
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
	print $out "$mem[0]\t$mem[1]\t$mem[2]\t$mem[4]\t$ARGV[3]\n" if $split{$key} > 0; #
	#print $out "$mem[0]\t$mem[1]\t$mem[2]\t$mem[3]\t$mem[4]\t$split{$key}\t$disU\t$splM\t$disM\t$ARGV[3]\n" if $split{$key} > 0; # t$support{$key}\t$reads\t$subseq\n" if $split{$key} > 0;



	}

}

`rm $ARGV[3].flash.*`;


`samtools fasta $ARGV[3].discord.sam -N -1 $ARGV[3].dis_1.fa -2 $ARGV[3].dis_2.fa`;

`samtools fasta $ARGV[3].split.sam -N -s $ARGV[3].split_s.fa`;

`cat $ARGV[3].dis_1.fa $ARGV[3].dis_2.fa $ARGV[3].split_s.fa > $ARGV[3].displit.fa`;

`makeblastdb -in $ARGV[4] -dbtype nucl -out $ARGV[4]`;

`blastn -db $ARGV[4] -word_size 50 -query $ARGV[3].displit.fa -out $ARGV[3].displitVStDNA.txt -max_target_seqs 10 -evalue 0.00001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore'`;

open my $dis1,"$ARGV[3].dis_1.fa";
open my $dis2,"$ARGV[3].dis_2.fa";
open my $spl,"$ARGV[3].split_s.fa";


my (%dis_1,%dis_2,%split_s);


	&Getfasta($dis1,\%dis_1);
	&Getfasta($dis2,\%dis_2);
	&Getfasta($spl,\%split_s);

open my $blastn,"$ARGV[3].displitVStDNA.txt";
my %uniq;

open my $dis1out,">$ARGV[3].dis_1.extra.fa";
open my $dis2out,">$ARGV[3].dis_2.extra.fa";
open my $splout,">$ARGV[3].split_s.extra.fa";


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

`cat $ARGV[3].dis_1.extra.fa $ARGV[3].dis_2.extra.fa $ARGV[3].split_s.extra.fa  > $ARGV[3].displit.extra.fa`;

`makeblastdb -in $ARGV[0] -dbtype nucl -out $ARGV[0]`;

`blastn -db $ARGV[0] -query $ARGV[3].displit.extra.fa -out $ARGV[3].displit.extraVSGenome.txt -max_target_seqs 10 -evalue 0.00001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore'`;


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
