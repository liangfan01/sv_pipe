#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
my $nanopore=shift;
my $list=shift;

my %hash;
open LIST,$list;
while(<LIST>){
	chomp;
	$hash{$_}=1;
}
close LIST;

if($nanopore=~/gz$/){
	open IN,"zcat $nanopore |";
}else{
	open IN,$nanopore;
}
while(<IN>){
	chomp;
	my $head=$_;
	my $seq=<IN>;
	my $strand=<IN>;
	my $q=<IN>;
	my @list=split /\s+/,$head;
	$list[0]=~s/^@//;
	if(exists $hash{$list[0]}){
		print "$head\n$seq$strand$q";
	}
}
close IN;
