#/bin/sh
bam=$1
s=$2
outfile=$3
export PATH=/home/yangqi/workdir/softwares/Sniffles-1.0.6/bin/sniffles-core-1.0.6:${PATH}
sniffles -t 20 -l 50 -s $s -m $bam -v $outfile
