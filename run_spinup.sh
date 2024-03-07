#!/bin/bash -l
### Job Name
#PBS -N run_spinup
### Project Code Allocation
#PBS -A UMCP0020
### Resources
#PBS -l select=1:ncpus=1:ngpus=1:mem=8GB
### Run Time
#PBS -l walltime=6:00:00
### Type of GPU
#PBS -l gpu_type=v100
### To the Casper queue
#PBS -q casper
### Log file
#PBS -o Logs/run_spinup.log
### Join output and error streams into single file
#PBS -j oe
### Email
#PBS -M zhihua@umd.edu
### Send email on abort, begin and end
#PBS -m abe

### Clear and load all the modules needed
module --force purge
module load ncarenv

export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Run spinup simulation 
proj_dir=/glade/u/home/zhihuaz/Projects/TRACE-SEAS/FrontalZone
#--project=<...> activates julia environment
julia --project=$proj_dir ./Simulations/spinup.jl

### Overwrite previous log file
LOG=$proj_dir/Logs/spinup.log
if [ -f "$LOG" ]; then
    rm -f $LOG
fi
mv $proj_dir/Logs/run_spinup.log $LOG

qstat -f $PBS_JOBID >> Logs/spinup.log
