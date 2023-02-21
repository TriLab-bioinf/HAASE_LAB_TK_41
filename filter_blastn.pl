#!/usr/bin/perl
use strict;

#pTet_pBac_LE35_LE35    f799ca50-ac73-4904-b442-56dcce778329    95.109  1472    30  31  1   1451    3904    2454    1451    4093

my %NR; # keep track of redundant read ids to avoid duplicated reads

while(<>){
	chomp; 
	next if m/^#/; 
    my @x=split /\t/; 
	
    my ($end5, $end3);
	my $PIGGYBAC_5 = 2;
	my $PIGGYBAC_3 = 1450;
	my $WINDOW = 50;
	my $strand = "+";
	if ($x[8] > $x[9]){
        #($x[8],$x[9]) = ($x[9],$x[8]);
		$strand = "-";
	}
	
	if ($x[3]>700){
        	if ($x[6] <= $PIGGYBAC_5 and (($x[8] >= 100 and $strand eq '+') or ( ($x[11] - $x[8]) >= 100 and $strand eq '-') ) ){
                $NR{$x[1]}++; next if $NR{$x[1]} > 1;
                # Extract region 5' to the insertion of PiggyBac
		        if ($strand eq '+'){
                    $end5 = $x[8] - 50;
                    $end3 = $x[8];
                    print "$x[1]\t$end5\t$end3\t$x[8]\n";
                } else  {
                    # subject on negative strand
                    $end5 = $x[8];
                    $end3 = $x[8] + 50;
                    print "$x[1]\t$end5\t$end3\t$x[8]\n";
                }

        	} elsif ($x[7] >= $PIGGYBAC_3 and ((($x[11] - $x[9])>= 100 and $strand eq '+') or ($x[9] >= 100 and $strand eq '-') ) ){
                $NR{$x[1]}++; next if $NR{$x[1]} > 1;    
                # Extract region 3' to the insertion of PiggyBac
                if ($strand eq '+'){
                        $end5 = $x[9];
                        $end3 = $x[9] + 50;
                        print "$x[1]\t$end5\t$end3\t$x[9]\n";
                } else {
                        # subject on negative strand
                        $end5 = $x[9] - 50;
                        $end3 = $x[9];
                        print "$x[1]\t$end5\t$end3\t$x[9]\n";
                }        
	        }
        }            
}

