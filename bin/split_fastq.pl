#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
die "perl $0 <seq file> <seq num>\n" unless (@ARGV==2);
my $seq=shift;
my $seq_num=shift;
if($seq=~/gz$/){
	open IN,"zcat $seq|";
}else{
	open IN,$seq;
}
my $num=0;

while(<IN>){
	chomp;
		
	$num++;
}
close IN;
