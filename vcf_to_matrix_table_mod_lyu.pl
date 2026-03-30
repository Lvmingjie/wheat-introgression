#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui_taas@126.com; chenrui.taas@gmail.com>
# Program Date:		2021.02.07	For Indel marker 
# Modified:			2021.12.07	Uniq Indel marker for each sample 
# Modified:			2022.02.13	For filtering wheat VCFs
# Modified:			2022.02.14	vct to matrix
############################################################

my $genome_fa = shift;
my $vcf = shift;
my $window = shift;
my $out	= shift;

die "
perl   xxx.pl   genome_fa   vcf.gz   10[M]   out_prefix

Input:
/aglab200T/aglab/cr/3_genomes/Wheat_IWGSC/V2.1/Assembly/iwgsc_refseqv2.1_assembly.fa
snp_XJ166_filt_merge.vcf.gz
10
matrix_table_XJ166

Output:
matrix_table_XJ166_10M_xxxxxx.xls


" if !defined $out;

my $date = `date +%Y%b%d%H%M`;
$date =~ s/[\r\n]+//g;

my($count, $d_num, $r_num, $total_number);
my(@t, $i, $count_zero, $z, $w); 
my($chr, $site, $ref, $target, $var, @m, $geno_d, $geno_r);
my(@ID_group, $tag, %chr_list, %seq, %chr_length, $j, $step);
my($head_line, %out_line, $end, %matrix, %count_region, %color);
my($region, $k, %cc, $gg, @nn, $dp);
############################################################


## Record genome_fa
open (IN, "< $genome_fa") or die $!;
while (<IN>) {
    $_ =~ s/[\r\n]+//g;
	
	if(/^>(\S+)/){  #ChrUnknown
		print "Recording the Chr: $1\r";
		$tag = $1;
		$chr_list{$tag} ++;
	}else{
		$seq{$tag} .= $_;
	}
}
close IN;
print "\n";


## Record vcf and output
if($vcf =~ /\.gz/){
	open(IN, "zcat $vcf |") or die $!;	##T630
}else{
	open (IN, "< $vcf") or die $!;
}

while (<IN>) {
    $_ =~ s/[\r\n]+//g;
	
	$count ++;
	print "Recording vcf lines:\t$count\r" if $count % 10000 == 0;
	
	undef @t;
	@t = split("\t", $_);
	
	if(/^##/){
		#print OUT $_."\n";
	}elsif(/^#CHROM/){

		if($t[9] =~ /^Donor/ ){
			$d_num = 9;
		}else{
			print "Check the Donor ID!\n";
		}
		if($t[10] =~ /^SHI4185/ ){
			$r_num = 10;
		}else{
			print "Check the Receptor ID!\n";
		}
		#print OUT $_."\n";
		
		#$total_number = $#t +1 -11;
		@ID_group = @t;
		
		#$head_line = "Region\tColor";
		#for $j(9..$#t){
		#	$head_line .= "\t".$t[$j];
		#}
		
	}elsif(/^Chr/){
		
		$chr 	= $t[0];
		$site 	= $t[1];
		#$ref 	= $t[3];
		#$var	= $t[4];
		
		$step = int($site/($window*1000000)) +1;
		$region = $t[0]."_".$step;
		$count_region{$region} ++;
		
		if($chr =~ /Chr\dA/){
			$color{$chr} = 7;
		}elsif($chr =~ /Chr\dB/){
			$color{$chr} = 8;
		}elsif($chr =~ /Chr\dD/){
			$color{$chr} = 9;
		}
		
		undef @m;
		@m = split(":", $t[$d_num]);
		$geno_d = $m[0];
		
		if($geno_d eq "0/0" || $geno_d eq "0|0") {
			for $j(9..$#t){
				undef @m;
				@m = split(":", $t[$j]);
				$dp = $m[2];
				
				if($m[0] eq "0/0" || $m[0] eq "0|0") {
					
					if($dp >= 2){
						push @{$matrix{$region."_".$ID_group[$j]}}, "1";
					}else{
						push @{$matrix{$region."_".$ID_group[$j]}}, "N";
					}
				
				}elsif($m[0] eq "1/1" || $m[0] eq "1|1"){
					if($dp >= 2){
						push @{$matrix{$region."_".$ID_group[$j]}}, "0";
					}else{
						push @{$matrix{$region."_".$ID_group[$j]}}, "N";
					}
					
				}elsif($m[0] eq "./." || $m[0] eq ".|."){
					push @{$matrix{$region."_".$ID_group[$j]}}, "N";
				}else{
					#push @{$matrix{$region."_".$ID_group[$j]}}, "0";
				}
			}
			
		}elsif($geno_d eq "1/1" || $geno_d eq "1|1"){
			for $j(9..$#t){
				undef @m;
				@m = split(":", $t[$j]);
				$dp = $m[2];
				
				if($m[0] eq "1/1" || $m[0] eq "1|1") {
					
					if($dp >= 2){
						push @{$matrix{$region."_".$ID_group[$j]}}, "1";
					}else{
						push @{$matrix{$region."_".$ID_group[$j]}}, "N";
					}
					
				}elsif($m[0] eq "0/0" || $m[0] eq "0|0"){
					
					if($dp >= 2){
						push @{$matrix{$region."_".$ID_group[$j]}}, "0";
					}else{
						push @{$matrix{$region."_".$ID_group[$j]}}, "N";
					}
					
				}elsif($m[0] eq "./." || $m[0] eq ".|."){
					push @{$matrix{$region."_".$ID_group[$j]}}, "N";
				}else{
					#push @{$matrix{$region."_".$ID_group[$j]}}, "0";
				}
			}
		}
		#print OUT "\n";
		
	}
}
close IN;
print "\n";
undef $count;


#open (OUT, "| gzip > $out") or die $!;
open (OUT, "> $out"."_".$window."M_".$date.".xls") or die $!;
open (OUT2, "> $out"."_".$window."M_".$date."_detailed.xls") or die $!;

print OUT "Region\tColor";
print OUT2 "Region\tColor";
for $i(9..$#ID_group){
	print OUT "\t".$ID_group[$i];
	print OUT2 "\t".$ID_group[$i];
}
print OUT "\n";
print OUT2 "\n";

foreach $i(sort keys %chr_list){
	
	if($i =~ /^Chr\d/){
	$chr_length{$i} = length($seq{$i});
	print $i."\t".$chr_length{$i}."\n";
	
	$end = int($chr_length{$i}/($window*1000000)) +1;
	for $j(1..$end){
		$region = $i."_".$j;
		print OUT $region."\t".$color{$i};
		print OUT2 $region."\t".$color{$i};
		
		if($count_region{$region}){
			
			for $k(9..$#ID_group){
				undef %cc;
				undef $gg;
				foreach $z(@{$matrix{$region."_".$ID_group[$k]}}){
					$cc{$z} ++;
				}
				
				if(!defined %cc){ #使用  !defined keys %cc  不报错但输出少了一些列
					#print OUT "\t0";
					print OUT "\tN";
					print OUT2 "\tN";
				}else{
				
					foreach $w (sort {$cc{$b} <=> $cc{$a}} keys %cc){
						$gg .= $w."_".$cc{$w}.";";
					}
					
					##策略1(原始): 输出占比最多的;
					##策略2(选用): 最多的如果是N, 则替换成第二多的非N;
					##策略3(待定): 优先级 1 > 0 > N, 多元情况下输出优先级高的;
					##策略4（）：1最多，输出1；0最多，1>=2 时，输出1；N最多，如果N+1 > count0, 输出1，如果0 >= 2,输出0，否则，输出N
					
					undef @nn;
					@nn = sort {$cc{$b} <=> $cc{$a}} keys %cc;
					
					if (scalar @nn > 1){
						if($nn[0] eq "N"){
							if ($cc{"1"} >=1 && $cc{"N"}+$cc{"1"} > $cc{"0"}){
								print OUT "\t1";
								print OUT2 "\t1".":".$gg;
							}
							elsif ($cc{"0"} >= 2){
								print OUT "\t0";
								print OUT2 "\t0".":".$gg;
							}
							else{
								print OUT "\tN";
								print OUT2 "\tN".":".$gg;
							}
						}elsif ($nn[0] eq "1"){
							print OUT "\t".$nn[0];
							print OUT2 "\t".$nn[0].":".$gg;
						}
						else{
							if ($cc{"1" >= 2}){
								print OUT "\t".$nn[0];
								print OUT2 "\t".$nn[0].":".$gg;

							}else{
							print OUT "\t".$nn[0];
							print OUT2 "\t".$nn[0].":".$gg;
								}
							}
					}else{
						print OUT "\t".$nn[0];
						print OUT2 "\t".$nn[0].":".$gg;
					}
					
					#foreach $w (sort {$cc{$b} <=> $cc{$a}} keys %cc){
					#	print OUT "\t".$w;					##0与1一样多时, 输出的是1
					#	print OUT2 "\t".$w.":".$gg; 		##0与1一样多时, 输出的是1
					#	last;
					#}
					
					
					##只要出现1次即认定为渗入;
					#if((scalar keys %cc) > 1){
					#	print OUT "\t1";
					#}
					
					##屏幕输出1个windw内的多个结果, 辅助判断; 
					#if((scalar @{$matrix{$region."_".$ID_group[$k]}}) > 1){
					#	print $region;
					#	foreach $w (@{$matrix{$region."_".$ID_group[$k]}}){
					#		print "\t".$w;
					#	}
					#	print "\n";
					#}
				}
			}
			print OUT "\n";
			print OUT2 "\n";
			
		}else{
			print OUT "\t1\t0";
			print OUT2 "\t1\t0";
			#print OUT "\tN\tN";
			
			for $j(11..$#ID_group){
				#print OUT "\t0";
				print OUT "\tN";
				print OUT2 "\tN";
			}
			print OUT "\n";
			print OUT2 "\n";
		}
	}
	
	}
}
close OUT;
print "\n";

