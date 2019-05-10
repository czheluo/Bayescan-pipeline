#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$pop,$single);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"pop:s"=>\$pop,
	"single:s"=>\$single,
			) or &USAGE;
&USAGE unless ($vcf and $out and $pop);
$vcf=ABSOLUTE_DIR($vcf);
$pop=ABSOLUTE_DIR($pop);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$dsh||="$out/work_sh";
mkdir $dsh if (!-d $dsh);
$single||=1;
if ($single eq 1) {
	open SH,">$dsh/bayes.fst.sh";
	open In,$pop;
	my %group;
	my %filehand;
	my %pophand;
	while (<In>) {
		chomp;
		next if ($_ eq ""|| /^$/);
		my ($id,$gid)=split(/\s+/,$_);
		if (!exists $filehand{$gid}) {
			open $filehand{$gid},">$out/$gid.list";
		}
		if (!exists $pophand{$gid}) {
			open $pophand{$gid},">$out/$gid.pophand";
		}
		print {$filehand{$gid}} "$id\n";
		print {$pophand{$gid}} "$id\t$gid\n";
	}
	close In;
	my @groid=sort keys %filehand;
	for (my $i=0;$i<@groid;$i++) {
		open SPID,">$out/$groid[$i].pid";
		print SPID "PARSER_FORMAT=VCF\n";
		print SPID "VCF_PARSER_QUAL_QUESTION=\n";
		print SPID "VCF_PARSER_POP_FILE_QUESTION=$out/$groid[$i].pophand\n";
		print SPID "VCF_PARSER_PLOIDY_QUESTION=DIPLOID\n";
		print SPID "VCF_PARSER_POP_QUESTION=true\n";
		print SPID "VCF_PARSER_GTQUAL_QUESTION=\n";
		print SPID "VCF_PARSER_MONOMORPHIC_QUESTION=false\n";
		print SPID "VCF_PARSER_IND_QUESTION=\n";
		print SPID "VCF_PARSER_REGION_QUESTION=\n";
		print SPID "VCF_PARSER_READ_QUESTION=\n";
		print SPID "VCF_PARSER_PL_QUESTION=false\n";
		print SPID "VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false\n";
		print SPID "WRITER_FORMAT=GESTE_BAYE_SCAN\n";
		print SPID "GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP\n";
		close SPID; 
		print SH "java -jar $Bin/bin/PGDSpider2-cli.jar -inputformat vcf -outputformat GESTE_BAYE_SCAN -spid $out/$groid[$i].pid -inputfile $vcf -outputfile $out/$groid[$i].bayscan && ";
		print SH "/mnt/ilustre/users/dna/.env/bin/bayescan_2.1 $out/$groid[$i].bayscan -snp -od $out -o $groid[$i]  -threads 8 -pr_odds 10 -out_pilot -out_freq -pilot 500 -burn 5000 && ";
		print SH "Rscript $Bin/bin/bayes.R --infile $out/$groid[$i]_fst.txt --outfile $out/$groid[$i].bayes\n";
	}
	close SH;
}else{
	open SH,">$dsh/fst.sh";
	open In,$pop;
	my %group;
	my %filehand;
	my %pophand;
	while (<In>) {
		chomp;
		next if ($_ eq ""|| /^$/);
		my ($id,$gid)=split(/\s+/,$_);
		if (!exists $filehand{$gid}) {
			open $filehand{$gid},">$out/$gid.list";
		}
		if (!exists $pophand{$gid}) {
			open $pophand{$gid},">$out/$gid.pophand";
		}
		print {$filehand{$gid}} "$id\n";
		print {$pophand{$gid}} "$id\t$gid\n";
	}
	close In;
	my @groid=sort keys %filehand;
	for (my $i=0;$i<@groid;$i++) {
		for (my $j=$i+1;$j<@groid;$j++) {
			open SPID,">$out/$groid[$i]-$groid[$j].pid";
			print SPID "PARSER_FORMAT=VCF\n";
			print SPID "VCF_PARSER_QUAL_QUESTION=\n";
			print SPID "VCF_PARSER_POP_FILE_QUESTION=$out/$groid[$i]-$groid[$j].pophand\n";
			print SPID "VCF_PARSER_PLOIDY_QUESTION=DIPLOID\n";
			print SPID "VCF_PARSER_POP_QUESTION=true\n";
			print SPID "VCF_PARSER_GTQUAL_QUESTION=\n";
			print SPID "VCF_PARSER_MONOMORPHIC_QUESTION=false\n";
			print SPID "VCF_PARSER_IND_QUESTION=\n";
			print SPID "VCF_PARSER_REGION_QUESTION=\n";
			print SPID "VCF_PARSER_READ_QUESTION=\n";
			print SPID "VCF_PARSER_PL_QUESTION=false\n";
			print SPID "VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false\n";
			print SPID "WRITER_FORMAT=GESTE_BAYE_SCAN\n";
			print SPID "GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP\n";
			close SPID; 
			print SH "cat $out/$groid[$i].pophand $out/$groid[$j].pophand > $out/$groid[$i]-$groid[$j].pophand && ";
			print SH "java -jar $Bin/bin/PGDSpider2-cli.jar -inputformat vcf -outputformat GESTE_BAYE_SCAN -spid $out/$groid[$i]-$groid[$j].pid -inputfile $vcf -outputfile $out/$groid[$i]-$groid[$j].bayscan && ";
			print SH "/mnt/ilustre/users/dna/.env/bin/bayescan_2.1 $out/$groid[$i]-$groid[$j].bayscan -snp -od $out -o $groid[$i]-$groid[$j]  -threads 8 -pr_odds 10 -out_pilot -out_freq && ";
			print SH "Rscript $Bin/bin/bayes.R --infile $out/$groid[$i]-$groid[$j]_fst.txt --outfile $out/$groid[$i]-$groid[$j].bayes\n";
		}
	}
	close SH;
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        men.luo\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl bayescan.fst.pl -vcf pop.recode.vcf -out ./ -pop group.list -single 1

Usage:
  Options:

  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -pop	<str>	group list
  -single   for the single group (default 1)
  -h         Help

USAGE
        print $usage;
        exit;
}
