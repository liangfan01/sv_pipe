A pipeline for sv calling
##step 1 Add lib to PYTHONPATH
```shell
export PYTHONPATH=/path/pipeline/:$PYTHONPATH
```shell
##step 2 run run_sv_pipe.py
```shell
python3 run_sv_pipe.py --help
usage: nanopore SV pipeline [-h] --list LIST --outdir OUTDIR
                            [--sv_pipe_home SV_PIPE_HOME]
                            [--method {ngmlr,last,ngmlr-last,ngmlr+last}]

optional arguments:
  -h, --help            show this help message and exit
  --list LIST, -l LIST  fastq list, each sample a list file
  --outdir OUTDIR, -o OUTDIR
                        output dirctory
  --sv_pipe_home SV_PIPE_HOME, -s SV_PIPE_HOME
                        sv pipeline home dirctory
  --method {ngmlr,last,ngmlr-last,ngmlr+last}, -m
{ngmlr,last,ngmlr-last,ngmlr+last}
                        ngmlr, only use ngmlr; last, only use last;
                        ngmlr+last, use both ngmlr and last, finally two
                        resluts; ngmlr-last, use ngmlr to find reads wiht
                        break points and use last to align them to the
                        reference.

```
