#!/usr/local/bin/perl -w
# Dinh Diep
use strict;

#($#ARGV == 1) or die "Usage: $0 <sequence> <edit distance>\n";

#my ($s1, $n) = (@ARGV);

#print "The number of sequence with substitution edit distance $n is: " . enumerate($s1, $n) . "\n";



# Return the number of sequence with distance n 
# 
#
sub enumerate
{
    # $s1 is the sequence
    # $len1 is the length of the sequence
    #
    my ($s1, $n) = @_;
    my ($len1) = length $s1;
    $s1 = uc($s1);
    return 0 if ($len1 == 0 or $n == 0);

    my %mat;
    $mat{$s1} = "NA";
    my @dict = ("A", "T", "C", "G");
    for(my $i = 0; $i < $n; $i++){
      my %enumerated;
      foreach my $string (keys %mat){
        my @edits = split ":", $mat{$string};
	my %mods;
	foreach my $mod (@edits){
          $mods{$mod} = 1;
	}
	if(scalar(@edits) == $i or $i == 0){
          for(my $j = 0; $j < $len1; $j++){
            my $cur_base = substr($s1, $j, 1);
	    next if($mods{$j});
            foreach my $alt_base (@dict){
              next if($alt_base eq $cur_base);
              my $cur_seq = substr($string, 0, $j) . $alt_base . substr($string, $j+1, $len1-$j);
	      #print $cur_seq, ",", $mat{$string}, ":$j\n";
	      if($mat{$string} eq "NA"){
                $mat{$cur_seq} = $j;
	      }else{
	        $mat{$cur_seq} = $mat{$string}.":".$j;
	      } 
	    }
	  }
	}
      }
    }
    my @enum;
    foreach my $string (keys %mat){
      $mat{$string} =~ s/NA://g;
      next if($string eq $s1);
      my @edits = split ":", $mat{$string};
      push(@enum, $string . "\t" . $s1 . ":". scalar(@edits));
    }
    return @enum;
}

1;
