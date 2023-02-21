#!/usr/bin/perl
use strict;

# Cluster_8 4_1258:1508-1607    -
# Cluster_8 4_2692:2372-2471    -
# Cluster_8 4_3494:3632-3731    +
# Cluster_8 4_4012:1883-1982    +
# Cluster_8 4_4106:1508-1607    -

my $usage = "$0 -c <cluster_file> -f <fasta file> -p output_prefix [def=cdhit]\n\n";
my %arg = @ARGV;
die $usage unless $arg{-f} and $arg{-c};

my $PREFIX = $arg{-p} || "cdhit";

# Load cluster info
open(CLUSTER, "<$arg{-c}") || die "ERROR, I cannot open $arg{-c}: $!\n\n";
my ($id, %cluster_ids, %cluster_strands);
while(<CLUSTER>){
    chomp;
    if(m/^>(Cluster.\d+)$/){
            $id = $1;
            $id =~ s/\s/_/;
    }elsif(m/>(\d_\d+:.+)\.{3}.+(([\+-])|(\*))/){
            push @{$cluster_ids{$id}}, $1;
            my $strand = $2 ;
            $strand =~ s/\*/\+/;
            push @{$cluster_strands{$id}}, $strand;
    }
}
close CLUSTER;

# load sequences
my (%seqs, $id); 
open(SEQ, "<$arg{-f}") || die "ERROR, I cannot open $arg{-f}: $!\n\n";
while(<SEQ>){
    chomp;
    if(m/^>(.+)/){
        $id = $1;
    }elsif(m/^([ATGCN]+)$/){
        $seqs{$id} = $1;
    }
}
close SEQ;


# print cluster sequences
foreach my $cluster (keys %cluster_ids){
    open(CLUSTEROUT,">$PREFIX"."_$cluster.fasta");    
    my $idx = -1;    
    foreach my $seqid (@{$cluster_ids{$cluster}}){
        $idx++;    
        # get sequence
        my $seq = $seqs{$seqid};
        my $strand = $cluster_strands{$cluster}[$idx];
        if ($strand eq '-'){
            $seq = &revcomp($seq);
        }
        print CLUSTEROUT ">$seqid\n$seq\n";
    }
    close CLUSTEROUT;
}  
exit(0);

sub revcomp {
    my $seq = shift;
    my $rev=reverse($seq); 
    $rev=~tr/ACGTN/TGCAN/;
    return $rev;
}




