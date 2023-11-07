#!/bin/bash -l
### Job Name
#PBS -N run_frontal_zone_wind
### Project Code Allocation
#PBS -A UMCP0020
### Resources
#PBS -l select=1:ncpus=1:ngpus=1:mem=6GB
### Run Time
#PBS -l walltime=4:00:00
### Type of GPU
#PBS -l gpu_type=v100
### To the Casper queue
#PBS -q casper
### Log file
#PBS -o Production/run_wind.log
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
julia --project=$proj_dir ./Production/frontal_zone_wind.jl

### Overwrite previous log file
LOG=$proj_dir/Production/frontal_zone_wind.log
if [ -f "$LOG" ]; then
    rm -f $LOG
fi
mv $proj_dir/Production/run_wind.log $LOG
