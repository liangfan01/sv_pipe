prefix=$1
indir=$2
outdir=$3
sv_pipe=/data/liangfan/pipeline/sv/
echo ngmlr START `date`
$sv_pipe/bin/ngmlr -t 4 -x ont -r $sv_pipe/database/ngmlr_index/ucsc.hg19.fasta -q $indir/$prefix.fq | samtools sort -m 4G  -o $outdir/$prefix.ngmlr.sort.bam -
echo ngmlr END `date`
$sv_pipe/bin/samtools view $outdir/$prefix.ngmlr.sort.bam |awk '$0~/^@/ || $0~/SA/{print $1}'  |sort |uniq > $outdir/$prefix.split.list
$sv_pipe/bin/samtools view -f 4  $outdir/$prefix.ngmlr.sort.bam|awk '{print $1}' |sort |uniq > $outdir/$prefix.unmapped.list
cat $outdir/$prefix.unmapped.list $outdir/$prefix.split.list > $outdir/$prefix.split_unmapped.list
perl $sv_pipe/bin/get_nanopore_reads.pl $indir/$prefix.fq $outdir/$prefix.split_unmapped.list > $outdir/$prefix.split.fq

echo LAST START `date`
head -40000 $outdir/$prefix.split.fq  > $outdir/$prefix\_10000.fq
$sv_pipe/bin/last-train -Q 1 $sv_pipe/database/last_db/humdb $outdir/$prefix\_10000.fq > $outdir/$prefix\_10000.params
$sv_pipe/bin/lastal -P 4 -Q1 $sv_pipe/database/last_db/humdb  $outdir/$prefix.split.fq $outdir/$prefix\_10000.params > $outdir/$prefix.split.maf
$sv_pipe/bin/maf-convert -d sam $outdir/$prefix.split.maf |samtools view -bS - |samtools sort -m 4G -o $outdir/$prefix.split.sort.bam -
echo LAST END `date`
