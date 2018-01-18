#prefix=$1
if [ $# = 0 ]  
	then   
	echo usage:
	echo sh $0 list outdir
	exit 1  
fi 
list=$1
out_dir=$2

#set pipe home
sv_pipe=/data/liangfan/pipeline/sv/
mkdir $out_dir;

outdir=`readlink -f  $out_dir`
#creat outdir
mkdir $outdir/01_split_fastq/
mkdir $outdir/02_alignment/
mkdir $outdir/03_sv/

cat $list |while read line
do 
	indir=`dirname $line`
	prefix=`basename $line |sed 's/.fq.gz//'`
	echo $sv_pipe/bin/fastqDear -cs 10000 -f1 $indir/$prefix.fq.gz -o $outdir/01_split_fastq/$prefix >> $outdir/01_split_fastq/run_split.sh
	echo $outdir/01_split_fastq/$prefix >> $outdir/01_split_fastq/all.split.list
done

echo SPLIT START `date`
$sv_pipe/bin/qsub_mgr.pl -s $outdir/01_split_fastq/run_split.sh -l 1 -jp split -m 1g -q all.q -r 
echo SPLIT DONE `date`

cat $outdir/01_split_fastq/all.split.list |while read line 
do 
	ls $line/*fq |while read sub_line
	do 
		sub_prefix=`basename $sub_line |sed 's/.fq//'`
		sub_indir=`dirname $sub_line`
		sub_outdir=$outdir/02_alignment/$sub_prefix
		echo $outdir/02_alignment/$sub_prefix/$sub_prefix.split.sort.bam >> $outdir/02_alignment/total.bam.list
		mkdir $sub_outdir
		echo sh $sv_pipe/scripts/run_alignment.sh $sub_prefix $sub_indir $sub_outdir >> $outdir/02_alignment/run_alignment.sh
	done
done

echo alignment START `date`
$sv_pipe/bin/qsub_mgr.pl -s $outdir/02_alignment/run_alignment.sh -l 1 -jp aln -m 1g -q all.q -P 4 -r -mj 100
echo alignment DONE `date`

echo $sv_pipe/bin/samtools merge -h $sv_pipe/database/sam.head -b $outdir/02_alignment/total.bam.list $outdir/02_alignment/total.split.bam  > $outdir/03_sv/run_sv.sh
echo perl $sv_pipe/bin/nanosv/nanosv.pl  -t 40 -c 1 -sambamba $sv_pipe/bin/sambamba $outdir/02_alignment/total.split.bam \> $outdir/03_sv/nanosv.vcf >> $outdir/03_sv/run_sv.sh

echo call SV START `date`
$sv_pipe/bin/qsub_mgr.pl -s $outdir/03_sv/run_sv.sh -l 2 -jp sv -m 100g -q all.q -P 40 -r 
echo call SV DONE `date`

