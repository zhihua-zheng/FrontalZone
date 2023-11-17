#!/bin/bash -l
#PBS -A UMCP0020
#PBS -N run_h11_Q000_W022_D000_St1
#PBS -o Logs/run_h11_Q000_W022_D000_St1.log
#PBS -j oe
#PBS -l walltime=12:00:00
#PBS -q casper
#PBS -l select=1:ncpus=1:ngpus=1:mem=6GB
#PBS -l gpu_type=v100
#PBS -M zhihua@umd.edu
#PBS -m abe

### Clear and load all the modules needed
module --force purge
module load ncarenv

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run simulation
proj_dir=/glade/u/home/zhihuaz/Projects/TRACE-SEAS/FrontalZone
time julia --check-bounds=no --project=$proj_dir Simulations/frontal_zone.jl "h11_Q000_W022_D000_St1"
###2>&1 | tee Logs/h11_Q000_W022_D000_St1.out

### Overwrite previous log file
LOG=$proj_dir/Logs/h11_Q000_W022_D000_St1.log
if [ -f "$LOG" ]; then
    rm -f $LOG
fi
mv $proj_dir/Logs/run_h11_Q000_W022_D000_St1.log $LOG

qstat -f $PBS_JOBID >> Logs/h11_Q000_W022_D000_St1.log
