#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
nanopore SV pipe
"""

# ---------
# Change Logs:
#
# ---------
__email__ = 'liangfan@grandomics.com'
__version__ = '0.0.1'
__status__ = 'Beta'
import os
import sys
import argparse
import logging
import time
import re
from sv.lib.system import *
from datetime import datetime


class Sv_pipe(object):
	def __init__(self, sv_pipe_home, prefix, output_dir):
		self.sv_pipe_home = sv_pipe_home
		self.prefix =  prefix
		self.output_dir =  output_dir
	
	def create_ngmlr_align_scripts(self):
		cmd = """echo ngmlr START `date`
{sv_home}/ngmlr \
 -t 4 \
 -x ont \
 -r {sv_home}/../database/ngmlr_index/ucsc.hg19.fasta \
 -q {1}/02_alignment/NGMLR/data/{0}.fq \
 | {sv_home}/samtools sort \
 -m 4G  \
 -o {1}/02_alignment/NGMLR/work/{0}.ngmlr.sort.bam -
echo ngmlr END `date`
""".format(self.prefix, self.output_dir, sv_home=self.sv_pipe_home)
		return cmd

	def create_last_align_scripts(self):
		cmd = """echo last START `date`
{sv_home}/lastal \
 -P 8 \
 -Q1 {sv_home}/../database/last_db/humdb  \
 {1}/02_alignment/LAST/data/{0}.fq \
 -p {sv_home}/../config/nanopore.params \
 > {1}/02_alignment/LAST/work/{0}.last.maf
{sv_home}/maf-convert -d sam {1}/02_alignment/LAST/work/{0}.last.maf \
 | {sv_home}/samtools view -bS - \
 | {sv_home}/samtools sort -m 4G -o {1}/02_alignment/LAST/work/{0}.last.sort.bam - \
 && rm  {1}/02_alignment/LAST/work/{0}.last.maf
echo last END `date`
""".format(self.prefix, self.output_dir, sv_home=self.sv_pipe_home)
		return cmd

	def create_breakpoint_reads_scripts(self):
		cmd = """echo find reads with break point START `date`
{sv_home}/samtools view {1}/02_alignment/NGMLR/work/{0}.ngmlr.sort.bam \
 |awk '$0~/^@/ || $0~/SA/{{print $1}}'  \
 |sort \
 |uniq \
 > {1}/02_alignment/LAST/data/{0}.split.list
{sv_home}/samtools view -f 4  {1}/02_alignment/NGMLR/work/{0}.ngmlr.sort.bam \
 | awk '{{print $1}}' |sort |uniq > {1}/02_alignment/LAST/data/{0}.unmapped.list
cat {1}/02_alignment/LAST/data/{0}.split.list {1}/02_alignment/LAST/data/{0}.unmapped.list \
 > {1}/02_alignment/LAST/data/{0}.split_unmapped.list
perl {sv_home}/get_nanopore_reads.pl {1}/02_alignment/NGMLR/data/{0}.fq {1}/02_alignment/LAST/data/{0}.split_unmapped.list \
 > {1}/02_alignment/LAST/data/{0}.fq
echo find reads with break point DONE `date`
""".format(self.prefix, self.output_dir, sv_home=self.sv_pipe_home)
		return cmd

def create_fq_link(prefix, indir, outdir):
	cmd = """ln -s {1}/{0}.fq {2}/{0}.fq
""".format(prefix, indir, outdir)
	return cmd

def create_qsub_scripts(sv_pipe_home, l, jp, mv, smp, shell):
	cmd = """{0}/qsub_mgr.pl \
    -s {5} \
    -l {1} \
    -jp {2} \
    -m {3} \
    -q all.q \
    -P {4} \
    -mj 200
""".format(sv_pipe_home, l, jp, mv, smp, shell)
	return cmd

def read_file_list(file_list):
	file={}
	with open(file_list, 'r') as f:
		for line in f:
			line = line.strip()
			[file_dir, file_basename] = os.path.split(os.path.realpath(line))
			file_dir = os.path.realpath(file_dir)
			file[file_basename] = file_dir
	f.close()
	return file

def create_split_fq_scripts(file, sv_pipe_home, output_dir):
	#prefix = os.path.split(os.path.realpath(file))[1]
	prefix = re.sub('.fq.gz',"",os.path.split(os.path.realpath(file))[1])
	cmd = """{1}/fastqDear \
    -cs 400000 \
    -f1 {0} \
    -o {2}/01_split_fastq/{3}
""".format(file, sv_pipe_home, output_dir, prefix)
	return cmd


#.
#├── 01_split_fastq
#│   ├── all.split.list
#│   ├── BNP17L0008-A1-1206.pass_reads/
#│   │   └── BNP17L0008-A1-1206.pass_reads.01.fq
#│   ├── BNP17L0008-A2-1206.pass_reads/
#│   │   └── BNP17L0008-A2-1206.pass_reads.120.fq
#│   ├── run_split.sh
#├── 02_alignment
#│   ├── NGMLR
#│   │   ├── run_ngmlr_alignment.sh
#│   │   ├── work
#│   │   ├── total.ngmlr.bam
#│   │   ├── total.ngmlr.bam.list
#│   │   └── data
#│   └── LAST
#│       ├── run_last_alignment.sh
#│       ├── work
#│       ├── total.last.bam
#│       ├── total.ngmlr.bam.list
#│       └── data
#└── 03_sv
# ├── nanosv.vcf
# ├── sniffles.vcf
# ├── run_sniffles.sh
# └── run_nanosv.sh


#def run_job(file):

def write_scripts(file, cmd):
	f = open(file, "w")
	f.write(cmd)
	f.close()

def file_name(file_dir):   
	L=[]   
	for root, dirs, files in os.walk(file_dir):  
		for file in files:  
			if os.path.splitext(file)[1] == '.fq':  
				L.append(os.path.join(root, file))  
	return L  

def mkdir_tree(method, output_dir):
	mkdir (output_dir)
	mkdir (output_dir + '/01_split_fastq')
	mkdir (output_dir + '/02_alignment')
	mkdir (output_dir + '/03_sv')
	mkdir (output_dir + '/workflow')
	if re.match('ngmlr', method):
		mkdir (output_dir + '/02_alignment/NGMLR')
		mkdir (output_dir + '/02_alignment/NGMLR/data')
		mkdir (output_dir + '/02_alignment/NGMLR/work')
	if re.match(r'.*last', method):
		mkdir (output_dir + '/02_alignment/LAST')
		mkdir (output_dir + '/02_alignment/LAST/data')
		mkdir (output_dir + '/02_alignment/LAST/work')

def create_merge_bam_scripts(sv_pipe_home, file, outbam):
	cmd = """{0}/samtools merge \
    -h {0}/../database/sam.head \
    -b {1} \
    {2}
""".format(sv_pipe_home, file, outbam)
	return cmd

def call_sv_scripts(sv_pipe_home, inbam, method, outvcf):
	if method == "nanosv":
		cmd = """{0}/nanosv/nanosv.pl \
    -t 20 \
    -c 1 \
    -sambamba {0}/sambamba {1} \
    > {2}
""".format(sv_pipe_home, inbam, outvcf)
	if method == "sniffles":
		cmd = """{0}/sniffles \
    --report_BND \
	--ignore_sd \
	-q 0 \
	--genotype \
	-n  10 \
	-t 20 \
    -l 50 \
    -s 1 \
    -m {1} \
    -v {2}
""".format(sv_pipe_home, inbam, outvcf)
	return cmd
	#nanosv,sniffie

def get_args():
	BASE_DIR = os.path.split(os.path.realpath(sys.argv[0]))[0]
	parser = argparse.ArgumentParser(prog='nanopore SV pipeline')
	parser.add_argument('--list','-l', required=True, help='fastq list, each sample a list file')
	parser.add_argument('--outdir','-o', required=True, help='output dirctory')
	parser.add_argument('--sv_pipe_home', '-s', default=BASE_DIR, help='sv pipeline home dirctory')
	parser.add_argument('--method','-m',default='ngmlr-last',choices=['ngmlr','last','ngmlr-last','ngmlr+last'], help='ngmlr, only use ngmlr; last, only use last; ngmlr+last, \
		use both ngmlr and last, finally two resluts; ngmlr-last, use ngmlr to find reads wiht break points and use last to align them to the reference.')
	#  parser.add_argument("--verbose", help="increase output verbosity",
	#  action="store_false")
	return parser.parse_args()

def main():
	args = get_args()
	SV_PIPE_HOME = args.sv_pipe_home
	outdir = os.path.realpath(args.outdir)
	method = args.method
	fq_list = args.list
	#get file list
	fq_files = read_file_list(fq_list)
	##mkdir 
	mkdir_tree(method, outdir)

	##create cmd
	split_cmd = ''
	ngmlr_cmd = ''
	split_cmd = ''
	last_cmd = ''
	ln_ngmlr_cmd = ''
	ln_last_cmd = ''
	split_fq_files = {}
	last_bam_list = ''
	ngmlr_bam_list = ''
	ngmlr_align_cmd = ''
	last_align_cmd = ''
	break_point_cmd = ''
	step1_cmd = ''
	workflow_sniffies_cmd = ''
	workflow_nanosv_cmd = ''

	for prefix in fq_files:
		file = fq_files[prefix] + '/' + prefix
		prefix = re.sub('.fq.gz',"",prefix)
		split_cmd += create_split_fq_scripts(file, SV_PIPE_HOME, outdir)

####split fastq files
	write_scripts(outdir + '/01_split_fastq/run_split.sh' , split_cmd)
	step1_cmd += create_qsub_scripts(SV_PIPE_HOME, '1', 's_fq', '3G', '1', outdir + '/01_split_fastq/run_split.sh')
	write_scripts(outdir + '/workflow/run_step1.sh' , step1_cmd)
	os.system('sh ' + outdir + '/workflow/run_step1.sh > ' + outdir + '/workflow/run_step1.sh.log 2>&1 ')
	print ('sh ' + outdir + '/workflow/run_step1.sh > ' + outdir + '/workflow/run_step1.sh.log 2>&1 ')

####get split fastq list and create scripts
	for prefix in fq_files:
		prefix = re.sub('.fq.gz',"",prefix)
		split_fq_list = file_name(outdir + "/01_split_fastq/"+ prefix)
		for split_fq_file in split_fq_list:
			sub_basename = os.path.split(os.path.realpath(split_fq_file))[1]
			sub_dirname = os.path.split(os.path.realpath(split_fq_file))[0]
			sub_prefix =  re.sub('.fq',"",sub_basename)
			sv_pipe_molecular = Sv_pipe(SV_PIPE_HOME, sub_prefix, outdir)
			#ngmlr_cmd += sv_pipe_molecular.create_ngmlr_align_scripts()
			#last_cmd += sv_pipe_molecular.create_last_align_scripts()
			break_point_cmd += sv_pipe_molecular.create_breakpoint_reads_scripts() 
			ln_ngmlr_cmd += create_fq_link(sub_prefix, sub_dirname, outdir + "/02_alignment/NGMLR/data/")
			ln_last_cmd += create_fq_link(sub_prefix, sub_dirname,  outdir + "/02_alignment/LAST/data/")
			ngmlr_align_cmd += sv_pipe_molecular.create_ngmlr_align_scripts()
			last_align_cmd += sv_pipe_molecular.create_last_align_scripts()
			ngmlr_bam_list += outdir + "/02_alignment/NGMLR/work/" + sub_prefix + ".ngmlr.sort.bam\n"
			last_bam_list +=  outdir + "/02_alignment/LAST/work/" + sub_prefix + ".last.sort.bam\n"
	merge_last_align_cmd = create_merge_bam_scripts(SV_PIPE_HOME, outdir + "/02_alignment/total_LAST_bam.list", outdir + "/02_alignment/total_last.bam")
	merge_ngmlr_align_cmd = create_merge_bam_scripts(SV_PIPE_HOME, outdir + "/02_alignment/total_ngmlr_bam.list", outdir + "/02_alignment/total_ngmlr.bam") 
#call sv_method sniffles,nanosv
	nanosv_cmd = call_sv_scripts(SV_PIPE_HOME, outdir + "/02_alignment/total_last.bam", "nanosv",  outdir + "/03_sv/nanosv.vcf")
	sniffles_cmd = call_sv_scripts(SV_PIPE_HOME, outdir + "/02_alignment/total_ngmlr.bam", "sniffles",  outdir + "/03_sv/sniffles.vcf")

	if re.match(r'.*last', method, flags=0):
		write_scripts(outdir + '/02_alignment/LAST/run_LAST_align.sh', last_align_cmd)
		write_scripts(outdir + "/02_alignment/total_LAST_bam.list", last_bam_list)
		write_scripts(outdir + "/02_alignment/run_merge_LAST.sh", merge_last_align_cmd)
		write_scripts(outdir + '/03_sv/run_nanosv.sh', nanosv_cmd)



	if re.match('ngmlr', method, flags=0):
		write_scripts(outdir + '/02_alignment/NGMLR/run_NGMLR_align.sh', ngmlr_align_cmd)
		write_scripts(outdir + "/02_alignment/total_ngmlr_bam.list", ngmlr_bam_list)
		write_scripts(outdir + "/02_alignment/run_merge_ngmlr.sh", merge_ngmlr_align_cmd)
		write_scripts(outdir + '/03_sv/run_sniffles.sh', sniffles_cmd)

	if re.match('ngmlr', method, flags=0):
		write_scripts(outdir + '/02_alignment/NGMLR/ln.sh' , ln_ngmlr_cmd)
		workflow_sniffies_cmd += 'sh ' + outdir + "/02_alignment/NGMLR/ln.sh\n"
		workflow_sniffies_cmd += create_qsub_scripts(SV_PIPE_HOME, '3', 'ngm', '10G', '4', outdir + '/02_alignment/NGMLR/run_NGMLR_align.sh')
		workflow_sniffies_cmd += create_qsub_scripts(SV_PIPE_HOME, '1', 'mngm', '10G', '1', outdir + "/02_alignment/run_merge_ngmlr.sh")
		workflow_sniffies_cmd += "qsub -cwd -l vf=100G -q all.q -pe smp 20 " + outdir + "/03_sv/run_sniffles.sh\n"
		if re.match(r'.*-last', method, flags=0):
			write_scripts(outdir + '/02_alignment/LAST/ln.sh' , break_point_cmd)
			workflow_sniffies_cmd += create_qsub_scripts(SV_PIPE_HOME, '6', 'g_fq', '3G', '1', outdir + '/02_alignment/LAST/ln.sh')
			workflow_sniffies_cmd += create_qsub_scripts(SV_PIPE_HOME, '4', 'last', '10G', '8', outdir + '/02_alignment/LAST/run_LAST_align.sh')
			workflow_sniffies_cmd += create_qsub_scripts(SV_PIPE_HOME, '1', 'mngm', '10G', '1', outdir + "/02_alignment/run_merge_LAST.sh")
			workflow_sniffies_cmd += "qsub -cwd -l vf=100G -q all.q -pe smp 20 " + outdir + "/03_sv/run_nanosv.sh\n"	
		elif re.match(r'.*\+last', method, flags=0):
			write_scripts(outdir + '/02_alignment/LAST/ln.sh' , ln_last_cmd)
			workflow_nanosv_cmd += 'sh ' + outdir + "/02_alignment/LAST/ln.sh\n"
			workflow_nanosv_cmd += create_qsub_scripts(SV_PIPE_HOME, '4', 'last', '10G', '8', outdir + '/02_alignment/LAST/run_LAST_align.sh')
			workflow_nanosv_cmd += create_qsub_scripts(SV_PIPE_HOME, '1', 'mngm', '10G', '1', outdir + "/02_alignment/run_merge_LAST.sh")
			workflow_nanosv_cmd += "qsub -cwd -l vf=100G -q all.q -pe smp 20 " + outdir + "/03_sv/run_nanosv.sh\n"				

	else:
		write_scripts(outdir + '/02_alignment/LAST/ln.sh' , ln_last_cmd)
		workflow_nanosv_cmd += 'sh ' + outdir + "/02_alignment/LAST/ln.sh\n"
		workflow_nanosv_cmd += create_qsub_scripts(SV_PIPE_HOME, '4', 'last', '10G', '8', outdir + '/02_alignment/LAST/run_LAST_align.sh')
		workflow_nanosv_cmd += create_qsub_scripts(SV_PIPE_HOME, '1', 'mngm', '10G', '1', outdir + "/02_alignment/run_merge_LAST.sh")
		workflow_nanosv_cmd += "qsub -cwd -l vf=100G -q all.q -pe smp 20 " + outdir + "/03_sv/run_nanosv.sh\n"			
	write_scripts(outdir + '/workflow/run_sv_sniffies.sh' , workflow_sniffies_cmd)
	os.system('sh ' + outdir + '/workflow/run_sv_sniffies.sh > ' + outdir + '/workflow/run_sv_sniffies.sh.log 2>&1 &')
	print ('sh ' + outdir + '/workflow/run_sv_sniffies.sh > ' + outdir + '/workflow/run_sv_sniffies.sh.log 2>&1 &')
	if workflow_nanosv_cmd:
		write_scripts(outdir + '/workflow/run_sv_nanosv.sh' , workflow_nanosv_cmd)
		os.system('sh ' + outdir + '/workflow/run_sv_nanosv.sh > ' + outdir + '/workflow/run_sv_nanosv.sh.log 2>&1 &')
		print ('sh ' + outdir + '/workflow/run_sv_nanosv.sh > ' + outdir + '/workflow/run_sv_nanosv.sh.log 2>&1 &')


if __name__ == '__main__':
	main()
