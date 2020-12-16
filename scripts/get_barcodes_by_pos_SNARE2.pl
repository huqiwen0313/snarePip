# attach barcode to atac fastq files
# version 2, auto read barcode prefix
# Dinh Diep & Qiwen Hu

#!/usr/bin/perl -w
use strict;

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';

sub revComp{
  my $seq = shift;
  my $rcSeq='';
  for(my $i=0; $i<length($seq); $i++){
    $rcSeq = $rcTable{uc(substr($seq,$i,1))} . $rcSeq;
  }
  return $rcSeq;
}

my $read1_file = $ARGV[0];
my $read2_file = $ARGV[1];
my $read3_file = $ARGV[2];
my $prefix = $ARGV[3];
my $output = $ARGV[4];
my $link_path = $ARGV[5];
my %links;

# read link table
open(LINK, $link_path) or die;
while(<LINK>){
	chomp;
	next if(/^Experiment_ID/);
	next if(/SNARE2-R/ or /SPL-R/);
	my @line = split(/\t/, $_);
	if($line[2] eq "BICCN-H"){
	  $line[0] =~ s/-SNARE2-AC//;
	  $line[0] =~ s/-SPL-AC//;
	  $line[0] =~ s/_SPL-AC//;
	  $links{$line[0]} = $line[1];
	  #print "$line[0]\t$line[1]\n";
	}
	else{
	  $line[0] =~ s/SNARE2-AC//;
	  $line[0] =~ s/_P[\d]*_//;
	  $line[0] =~ s/__/_/g;
	  $links{$line[0]} = $line[1];
	  #print "$line[0]\t$line[1]\n";
	}
}
close LINK;

my $samplePrefix = $prefix;
$samplePrefix =~ s/-SPL.*Ad1/_/;
$samplePrefix =~ s/.Ad1//;
$samplePrefix =~ s/_S.*//;
$samplePrefix =~ s/.P*_Ad1//;
my $barcodePrefix = "";

if(exists($links{$samplePrefix})){
  $barcodePrefix = $links{$samplePrefix};
}

#print "$samplePrefix\t$barcodePrefix";
#exit;

open(R1, "zcat $read1_file |") or die("error reading $read1_file\n");
open(R2, "zcat $read2_file |") or die("error reading $read2_file\n");
open(R3, "zcat $read3_file |") or die("error reading $read3_file\n");

open(R1_OUT, ">>$output/$prefix.R1.fastq") or die("error writing to $prefix.R1.fastq");
open(R3_OUT, ">>$output/$prefix.R3.fastq") or die("error writing to $prefix.R3.fastq");

while(my $r1_1 = <R1>){

  my $r1_2 = <R1>;
  my $r1_3 = <R1>;
  my $r1_4 = <R1>;
  my $r2_1 = <R2>;
  my $r2_2 = <R2>;
  my $r2_3 = <R2>;
  my $r2_4 = <R2>;
  my $r3_1 = <R3>;
  my $r3_2 = <R3>;
  my $r3_3 = <R3>;
  my $r3_4 = <R3>;
  
  my @fields = split / /, $r1_1;
  
  my $barcode1 = substr($r2_2, 0, 8);
  my $barcode2 = substr($r2_2, 38, 8);
  my $barcode3 = substr($r2_2, 76, 8);
  my $umi = substr($r2_2, 84, 10);
  #print("$barcode1\t$barcode2\t$barcode3\n");
  #exit; 

  $r1_1=~ s/@//g;
  $r3_1=~ s/@//g;
  
  #my $cell_barcode = $barcode1 . $barcode2 . $barcode3 . $barcode4;
  my $cell_barcode = $barcode1 . $barcode2 . $barcode3;
  my $rc_cell_barcode = revComp($cell_barcode);

  print R1_OUT "@", $barcodePrefix, "_", $rc_cell_barcode, ":$umi:", $r1_1, $r1_2, $r1_3, $r1_4;
  print R3_OUT "@", $barcodePrefix, "_", $rc_cell_barcode, ":$umi:", $r3_1, $r3_2, $r3_3, $r3_4;

}
