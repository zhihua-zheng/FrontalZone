#!/bin/bash -l
### Job Name
#PBS -N run_eady_mwe
### Project Code Allocation
#PBS -A UMCP0020
### Resources
#PBS -l select=1:ncpus=1:ngpus=1:mem=4GB
### Run Time
#PBS -l walltime=2:00:00
### Type of GPU
#PBS -l gpu_type=v100
### To the Casper queue
#PBS -q casper
### Log file
#PBS -o Production/run.log
### Join output and error streams into single file
#PBS -j oe
### Email
#PBS -M zhihua@umd.edu
### Send email on abort, begin and end
#PBS -m abe

### Clear and load all the modules needed
module --force purge
module load ncarenv 

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run job
proj_dir=$HOME/Projects/TRACE-SEAS/FrontalZone
#--project=<...> activates julia environment
julia --project=$proj_dir ./Production/eady_mwe.jl

### Overwrite previous log file
LOG=$proj_dir/Production/eady_mwe.log
if [ -f "$LOG" ]; then
    rm -f $LOG
fi
mv $proj_dir/Production/run.log $LOG
