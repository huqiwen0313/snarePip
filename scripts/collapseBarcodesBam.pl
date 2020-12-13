#!/usr/bin/perl -w
# Dinh Diep
use strict;
use File::Basename;
my $dirname = dirname(__FILE__);
require "$dirname/enumerate_mm_edit.pl";

my $whitelist = $ARGV[0];
my $bam_file = $ARGV[1];
my $out_file = $ARGV[2];
my $max_distance = 2;

my %barcodes;
my %unplaced_barcodes;

open(WHITELIST, "$whitelist") || die("Error reading $whitelist\n");
while(my $line = <WHITELIST>){
  my @fields = split " ", $line;
  foreach my $s1 (@fields){
    $barcodes{$s1} = 1;
    my @enum = enumerate($s1, $max_distance);
    foreach my $en (@enum){
      my ($string, $poss) = split "\t", $en;
      push(@{$unplaced_barcodes{$string}}, $poss);
    }
  }
}
close(WHITELIST);

my $m = 0;
open(BAMFILE, "samtools view -H $bam_file |") || die("Error reading $bam_file\n");
open(OUTFILE, ">$out_file.sam") || die("Error writing to $out_file\n");
while(my $header_line = <BAMFILE>){
  print OUTFILE $header_line;
}
close(BAMFILE);

open(BAMFILE, "samtools view $bam_file |") || die("Error reading $bam_file\n");
while(my $sam_line = <BAMFILE>){
  chomp($sam_line);
  $m++;
  if($m%1000000==0){
    print "$m reads in bam\n";
  }
  my @fields = split "\t", $sam_line;
  my ($cell, $read_id) = split ":", $fields[0];
  $cell =~ s/_2/\.2/;
  $cell =~ s/_N/\.N/;
  $cell =~ s/_S/\.S/;
  my ($pool_id, $cell_barcode) = split "_", $cell;
  my $barcode1 = substr($cell_barcode, 0, 8);
  my $barcode2 = substr($cell_barcode, 8, 8);
  my $barcode3 = substr($cell_barcode, 16, 8);
  #print $barcode1, ",", $barcode2, ",", $barcode3, "\n";
  if(!exists($barcodes{$barcode1}) and exists($unplaced_barcodes{$barcode1})){
    my @Poss = @{$unplaced_barcodes{$barcode1}};
    next if(scalar(@Poss) > 1);
    my ($query, $edit) = split ":", $Poss[0];
    $barcode1 = $query;
  }
  if(!exists($barcodes{$barcode2}) and exists($unplaced_barcodes{$barcode2})){
    my @Poss = @{$unplaced_barcodes{$barcode2}};
    next if(scalar(@Poss) > 1);
    my ($query, $edit) = split ":", $Poss[0];
    $barcode2 = $query;
  }
  if(!exists($barcodes{$barcode3}) and exists($unplaced_barcodes{$barcode3})){
    my @Poss = @{$unplaced_barcodes{$barcode3}};
    next if(scalar(@Poss) > 1);
    my ($query, $edit) = split ":", $Poss[0];
    $barcode3 = $query;
  }
  if(exists($barcodes{$barcode1}) and exists($barcodes{$barcode2}) and exists($barcodes{$barcode3})){
    $cell = $pool_id . "_" . $barcode1 . $barcode2 . $barcode3;
    $fields[0] = $cell . ":" . $read_id;
    print OUTFILE join("\t", @fields), "\n";
  }
}
close(BAMFILE);
close(OUTFILE);
print "total reads in bam: ", $m, "\n";
system("samtools view -b $out_file.sam | samtools sort -n - > $out_file.bam");
system("rm $out_file.sam");
#print "total unplaced: ", scalar(keys %unplaced_barcodes), "\n";
#print "total barcodes: ", scalar(keys %barcodes), "\n";

