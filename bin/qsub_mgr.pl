#!/usr/bin/perl -w
=head1 Name
  qsub.pl

=head1 Description
  This program is based on the qsub-sge.pl written by Fan Wei and Hu Yujie, Many Thanks to both of them.
  1) cpu efficiency is no less than 0.3, [cpu time/uptime]
  2) the maximum overflowed_mem is 2G and the maximum uptime of this kind of jobs is 30 minutes

=head1 Usage
  perl qsub.pl [options]
  -s/shell      <str> * the job shell
  -l/lines      <int>   number of lines to form a job [1], nonblank lines
  -jp/jobprefix <str>   the prefix for qsubed jobs [work]
  -mj/maxjobs   <int>   the maximum number of jobs to throw out [30]
  -md/mindisk   <str>   the minimum limit of disk space [200G]
  -m/memory     <str> * the required resource used in qsub
  -q/queue      <str> * the queue to use
  -P/project    <str> * -pe smp 
  -r/reqsub             don't reqsub the unfinished jobs, default reqsub
  -i/internal   <int>   the interval time of checking by qstat
  -silence              don't output the process check information, default output
  -ce/cpu_eff   <float> cpu efficiency, cpu time/uptime [0]
  -h/help               help information

=head1 Example
  perl qsub.pl -s run.sh -q bc.q -P phartest -m 0.1m

=head1 Version 
  version: 2.0
  modified: 20121119, 20120814, 20120801
=cut
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
# get options
my ($lines,$job_prefix,$shell,$memory,$queue,$project,$maxjobs,$reqsub,$interval,$mindisk,$silence,$cpu_efficiency,$help);
GetOptions(
  "s|shell:s"=>\$shell,
  "l|lines:i"=>\$lines,
  "jp|jobprefix:s"=>\$job_prefix,
  "mj|maxjobs:s"=>\$maxjobs,
  "md|mindisk:s"=>\$mindisk,
  "m|memory:s"=>\$memory,
  "q|queue:s"=>\$queue,
  "P|project:s"=>\$project,
  "r|reqsub"=>\$reqsub,
  "i|interval:i"=>\$interval,
  "silence:s"=>\$silence,
  "ce|cpu_eff:f"=>\$cpu_efficiency,
  "h|help:s"=>\$help
);
#die`pod2text $0`if(!$shell || !$queue || !$project || !$memory || defined $help);
die`pod2text $0`if(!$shell || !$queue ||  !$memory || defined $help);
# default parameters
#
$project ||=1;
$lines||=1;
$interval||=295;
$maxjobs||=30;
$mindisk||="200G";
$mindisk=&TransUnit($mindisk);
$job_prefix||="work";
my @memory=split /,/,$memory;
$shell=abs_path($shell);
#print $shell,"\n";
die"file does not exist!\n" if(!-e $shell);
my $outdir=$shell."\.".$$.".qsub";
#print $outdir,"\n";
my $outlog=$shell."\.".$$.".log";
open LOG,"> $outlog"||die"failed to open $outlog\n";
system("mkdir -p $outdir"); #                                                +--- overflow
my $overflow_mem=2; # 2G                                                     |--- path
my $time='`date +%F'."'  '".'%H:%M`'; #                                      |--- status (-2-dr -1-qdel 0-error 1-r 1.5-hqw/qw/t 2-hold 3-finish)
my $max_circle=30; #                                          +--- jobid  ---|--- diskcheck
$cpu_efficiency||=0;      # cpu time / total time             |              |--- cpu
my $cpu_check=30;  # minutes                                  |              +--- submit
my $disk_check=720; # minutes                          Jobs---|
my $disk_mem=&CheckDisk; #                                    |              +--- circle
my %Jobs; #                                                   +--- script ---|--- lines
&Qsub; #                                                                     +--- memory 
                                                   
sub Debug{
  foreach my $script(keys %{$Jobs{'script'}}){
    print "$script\t".$Jobs{'script'}{$script}{'memory'}."\t".$Jobs{'script'}{$script}{'circle'}."\t".$Jobs{'script'}{$script}{'lines'}."\n";
  }
  foreach my $jobid(keys %{$Jobs{'jobid'}}){
    print $jobid,"\t",$Jobs{'jobid'}{$jobid}{'status'},"\t",$Jobs{'jobid'}{$jobid}{'cpu'},"\t",$Jobs{'jobid'}{$jobid}{'diskcheck'},"\t",$Jobs{'jobid'}{$jobid}{'uptime'},"\t",$Jobs{'jobid'}{$jobid}{'path'},"\n";
  }
}

sub Qsub{
  &ReadShell;
  my $run_num=0;
  chdir($outdir);
  foreach my $script(keys %{$Jobs{'script'}}){
    $run_num++;
    while(1){
      $run_num=&CheckJobs if($run_num>$maxjobs);
      if($disk_mem>$mindisk && $run_num<=$maxjobs){
        &CastTask($script);
        $Jobs{'script'}{$script}{'circle'}=1;
        last;
      }else{
        if($disk_mem<=$mindisk){
          &ShowLog("no enough disk space, wait for next throwing.") if(!defined $silence);
          sleep $interval;
          $disk_mem=&CheckDisk;
        }else{
          &ShowLog("task to limit, wait for next throwing.") if(!defined $silence);
          sleep $interval;
        }
      }
    }
  }
  while($run_num!=0){
    sleep $interval;
    $run_num=&CheckJobs;
  }
  my $flag=0;
  print LOG "\n****   task state   ****\n"; 
  foreach my $jobid(keys %{$Jobs{'jobid'}}){
    if($Jobs{'jobid'}{$jobid}{'status'}!=3){
      print LOG "$Jobs{'jobid'}{$jobid}{'path'} is uncompleted!\n";
      $flag++;
    }
  }
  if($flag==0){
    print LOG "All jobs completed!\n";
  }
  close LOG;
# &Debug;
}

sub CheckJobs{
  my $run_num=0;
  my %run_jobs;
  my $user=$ENV{"USER"} || $ENV{"USERNAME"};
  my @qstat=`qstat -u $user | awk '{if(/$user/)print \$1"\t"\$5}'`;
  for(@qstat){
    chomp;
    my @tmp=split /\t/,$_;
    next if(!exists $Jobs{'jobid'}{$tmp[0]});
    $run_jobs{$tmp[0]}="";
    if($tmp[1] eq "r"){
      $Jobs{'jobid'}{$tmp[0]}{'status'}=1;
    }elsif($tmp[1] eq "hqw"  || $tmp[1] eq "qw" || $tmp[1] eq "t"){
      $Jobs{'jobid'}{$tmp[0]}{'status'}=1.5;
    }elsif($tmp[1] eq "dr"){
      $Jobs{'jobid'}{$tmp[0]}{'status'}=-2;
    }else{
      $Jobs{'jobid'}{$tmp[0]}{'status'}=0;
    }
  }
  sleep 5;
  my $now_time=time;
  foreach my $jobid(keys %{$Jobs{'jobid'}}){
    if($now_time-$Jobs{'jobid'}{$jobid}{'diskcheck'}>=$disk_check*60){
      $disk_mem=&CheckDisk;
      $Jobs{'jobid'}{$jobid}{'diskcheck'}=$now_time;
    }
    if(($Jobs{'jobid'}{$jobid}{'status'}==1 || $Jobs{'jobid'}{$jobid}{'status'}==1.5 || $Jobs{'jobid'}{$jobid}{'status'}==2) && !exists $run_jobs{$jobid}){
      # check .o
      if(-f "$Jobs{'jobid'}{$jobid}{'path'}.o$jobid"){
        chomp(my $check_log=`grep -c This-Work-is-Completed $Jobs{'jobid'}{$jobid}{'path'}.o$jobid`);
        if($check_log==$Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'lines'}){
          $Jobs{'jobid'}{$jobid}{'status'}=3;
        }else{
          if(!defined $reqsub && $Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'circle'}<$max_circle){
            &Reqsub($jobid);
            ShowLog("<reqsub> *.o is uncompleted") if(!defined $silence);
            $run_num++;
          }
        }
      }else{
        if($Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'circle'}<$max_circle){
          &Reqsub($jobid);
          $run_num++;
        }
      }
    }elsif(($Jobs{'jobid'}{$jobid}{'status'}==1 || $Jobs{'jobid'}{$jobid}{'status'}==1.5) && exists $run_jobs{$jobid}){
      if($disk_mem<=$mindisk){
        `qhold $jobid`;
        &ShowLog("<hold> disk space is not enough") if(!defined $silence);
        $Jobs{'jobid'}{$jobid}{'status'}=2;
        $run_num++;
      }else{
        if($now_time-$Jobs{'jobid'}{$jobid}{'cpu'}>=$cpu_check*60){
          if(&CheckCPU($jobid)){
            $run_num++;
          }else{
            if(!defined $reqsub && $Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'circle'}<$max_circle){
              &Reqsub($jobid);
              $run_num++;
            }else{
              `qdel $jobid`;
            }
          }
          $Jobs{'jobid'}{$jobid}{'cpu'}=$now_time;
        }else{
          $run_num++;
        }
      }
    }elsif($Jobs{'jobid'}{$jobid}{'status'}==2 && exists $run_jobs{$jobid}){
      $disk_mem=&CheckDisk if($disk_mem<=$mindisk);
      if($disk_mem>$mindisk){
        `qrls $jobid`;
        $Jobs{'jobid'}{$jobid}{'status'}=1;
      }
      $run_num++;
    }elsif($Jobs{'jobid'}{$jobid}{'status'}==-1){
      my $qdel=`qdel $jobid`;
      delete $Jobs{'jobid'}{$jobid} if($qdel=~/denied/ || $qdel=~/does not exist/ || $qdel=~/deletion/);
    }elsif($Jobs{'jobid'}{$jobid}{'status'}==-2){
      Reqsub($jobid);
    }elsif($Jobs{'jobid'}{$jobid}{'status'}==0){
      if($Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'circle'}<$max_circle){
        &Reqsub($jobid);
        $run_num++;
      }else{
        `qdel $jobid`;
      }
    }
  }
  return $run_num;
}

sub Reqsub{
  my $jobid=shift;
  `qdel $jobid`;
  &CastTask($Jobs{'jobid'}{$jobid}{'path'});
  $Jobs{'jobid'}{$jobid}{'status'}=-1;
  system("rm -rf $Jobs{'jobid'}{$jobid}{'path'}.o$jobid $Jobs{'jobid'}{$jobid}{'path'}.e$jobid");
  $Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'circle'}++;
}

sub CastTask{
  my $script=shift;
  while(1){
#    my $job_info=`qsub -cwd -S /bin/sh -q $queue -P $project -l vf=$Jobs{'script'}{$script}{'memory'} $script`;
      my $job_info=`qsub -cwd -q $queue -pe smp $project -l vf=$Jobs{'script'}{$script}{'memory'} $script`;
	  print "$job_info\n";
      if($job_info=~/Your job (\d+)/){
      my $jobid=$1;
      $Jobs{'jobid'}{$jobid}{'status'}=1;
      $Jobs{'jobid'}{$jobid}{'diskcheck'}=$Jobs{'jobid'}{$jobid}{'cpu'}=$Jobs{'jobid'}{$jobid}{'submit'}=time;
      $Jobs{'jobid'}{$jobid}{'path'}=$script;
      last;
    }else{
      sleep $interval;
      redo;
    }
  }
}

sub CheckCPU{
  my $jobid=shift;
  my @cpu_info=split /\s+/,`qstat -j $jobid | grep usage`;
  if(@cpu_info){
    my $vmem=$1 if($cpu_info[-2]=~/vmem=(.*),$/);
    if($vmem=~/\d+.*/){
      $vmem=TransUnit($vmem);
      if($vmem-TransUnit($Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'memory'})>=$overflow_mem && exists $Jobs{'jobid'}{$jobid}{'overflow'}){
        ShowLog("<killed> job: $jobid   submit_mem: $Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'memory'}   vmem: $vmem"."G") if(!defined $silence);
        $Jobs{'script'}{$Jobs{'jobid'}{$jobid}{'path'}}{'memory'}=$vmem."G";
        return 0;
      }else{
        $Jobs{'jobid'}{$jobid}{'overflow'}="";
      }
    }
    my $cpu_time=$1 if($cpu_info[2]=~/cpu=(.*),$/);
	my @cpu_time=split /\:/,$cpu_time;
    my $now_time=time;
	my $minutes=0;
    if(@cpu_time>4){
      return 1;
    }else{
      if(@cpu_time==2){
        $minutes=$cpu_time[0];
      }elsif(@cpu_time==3){
        $minutes=$cpu_time[0]*60+$cpu_time[1];
      }elsif(@cpu_time==4){
        $minutes=$cpu_time[0]*12*60+$cpu_time[1]*60+$cpu_time[2];
      }
      if($minutes<($now_time-$Jobs{'jobid'}{$jobid}{'submit'})*$cpu_efficiency/60){
        &ShowLog("<killed> job: $jobid   cpu time: $minutes   uptime: ".int(($now_time-$Jobs{'jobid'}{$jobid}{'submit'})/60)) if(!defined $silence);
        return 0;
      }else{
        return 1;
      }
    }
  }else{
    return 1;
  }
}

sub ReadShell{
  my $job_mark="00000";
  my $line_mark=0;
  my $shell_mark=0;
  open IN,"< $shell"||die"failed open $shell";
  while(<IN>){
    chomp;
    next if($_=~/^\s+$/ || $_=~/\s+\n$/ || $_ eq ""); # the real lines in this file
    if($line_mark % $lines==0){
      $job_mark++;
      open OUT,"> $outdir/$job_prefix\_$job_mark\.sh"||die"failed creat $outdir/$job_prefix\_$job_mark\.sh";
      if(@memory>1){
        if($memory[$shell_mark]){
          $Jobs{'script'}{"$outdir/$job_prefix\_$job_mark\.sh"}{'memory'}=&TransUnit($memory[$shell_mark]).'G';
          $shell_mark++;
        }else{
          $Jobs{'script'}{"$outdir/$job_prefix\_$job_mark\.sh"}{'memory'}=&TransUnit($memory[-1]).'G';
        }
      }else{
        $Jobs{'script'}{"$outdir/$job_prefix\_$job_mark\.sh"}{'memory'}=&TransUnit($memory).'G';
      }
      $Jobs{'script'}{"$outdir/$job_prefix\_$job_mark\.sh"}{'lines'}=0;
      print OUT "echo start at time: $time \n";
    }
    s/;\s*$//;
    s/;\s*;/;/g;
    if($_=~/[\s+]?#+.*/ || $_=~/^\s+$/ || $_=~/\s+\n$/ || $_ eq ""){
      print OUT "echo \"$_\"\n";
    }else{
      $Jobs{'script'}{"$outdir/$job_prefix\_$job_mark\.sh"}{'lines'}++;
      print OUT $_.' &&  echo This-Work-is-Completed!'." Line ".$Jobs{'script'}{"$outdir/$job_prefix\_$job_mark\.sh"}{'lines'}."\n";
    }
    if($line_mark % $lines==$lines-1){
      print OUT "echo finish at time: $time\n";
      close OUT;
    }
    $line_mark++;
  }
  close IN;
  if($line_mark % $lines!=0){
    print OUT "echo finish at time: $time\n";
#    &ShowLog("<warning> There must be lose some rows of $shell");
  }
}
#sub CheckDisk{
# my $disk_info=`df -h $outdir`;
#  die"Cannot get the disk space info!\n" if($disk_info !~ /Filesystem.+Size.+Used.+Avail/i);
#  $disk_info=(split /\n/,$disk_info)[2];
#  $disk_info=(split /\n/,$disk_info)[1];
#  $disk_info=~s/^\s+//g;
#  $disk_info=~s/\s+$//g;
#  # $disk_info=(split /\s+/,$disk_info)[2];
#  $disk_info=(split /\s+/,$disk_info)[3];
#  $disk_info=&TransUnit($disk_info);
#  return $disk_info;
#}
sub CheckDisk{
	my $disk_info=10000;
	return $disk_info;
}
sub TransUnit{
  my $mem=shift;
  my ($num,$unit)=($mem=~/^([\d\.]+)([mgt])$/i);
  die"Illegal memory type <0.2G> $mem\n"unless($num && $unit);
  if($unit eq "g" || $unit eq "G"){
    return &DeleteZero($num);
  }elsif($unit eq "m" || $unit eq "M"){
    my $m=sprintf("%.3f",$num/1024);
    return &DeleteZero($m);
  }elsif($unit eq "t" || $unit eq "T"){
    my $t=sprintf("%.3f",$num*1024);
    return &DeleteZero($t);
  }
}

sub DeleteZero{
  my $num=shift;
  $num=~s/\.0+$// if($num=~/\.0+$/);
  $num=~s/0+$// if($num=~/\./);# delete the last zero
  $num=0.01 if($num<0.01);
  return $num;
}

sub ShowLog{
  my $info=shift;
  my @times = localtime; # sec, min, hour, day, month, year
  printf LOG ("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}
