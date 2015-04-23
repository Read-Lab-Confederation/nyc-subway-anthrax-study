#!/usr/bin/perl

# # # # # #
# extractPartOfFasta.pl
# written by Linnéa Smeds 3 Feb 2011
# ====================================================
# Extract a specific region from a sequence. Takes a 
# fasta file (can contain several sequences), a sequence
# header and start and end positions.
# NB! The sequence header must not contain spaces. If
# it does, give only the first tab, eg "contig001" if 
# header is ">contig001 chr=1 length=1000..."
# ====================================================
# Usage: extractPartOfFasta.pl <fastafile> <seq name>
#									<start> <end>
#
# Example: extractPartOfFasta.pl mySeq.fa contig4 200 \
#													299 >contig4_pos200-299.fa

use strict;
use warnings;


# Input parameters
my $fasta = $ARGV[0];
my $name = $ARGV[1];
my $start = $ARGV[2];
my $end = $ARGV[3];

my $noOfBases = $end-$start+1;

open(FAS, $fasta);
my ($seq, $head, $seqFlag) = ("", "", "off");
while(<FAS>) {
	if(/>/) {
		if($seqFlag eq "on" && $seq ne "") {
			my $len = length($seq);
			print "$head orig_len=$len, extract $start-$end ($noOfBases bp)\n";
			my $substr = substr($seq, $start-1, $noOfBases);
			my @seqParts = split(/(.{100})/, $substr);
			for my $seqs (@seqParts) {
				unless($seqs eq "") {
					print $seqs."\n";
				}
			}
			$seq="";
			$seqFlag="off";
		}
		my @tab = split(/\s+/, $_);
		if($tab[0] eq ">$name") {
			$seqFlag = "on";
			$head = $tab[0];
		}
	}
	else {
		if($seqFlag eq "on") {
			chomp($_);
			$seq.=$_;
		}
	}
}
if($seqFlag eq "on" && $seq ne "") {
	my $len = length($seq);
	print "$head orig_len=$len, extract $start-$end ($noOfBases bp)\n";
	my $substr = substr($seq, $start-1, $noOfBases);
	my @seqParts = split(/(.{100})/, $substr);
	for my $seqs (@seqParts) {
		unless($seqs eq "") {
			print $seqs."\n";
		}
	}
	$seq="";
	$seqFlag="off";
}
