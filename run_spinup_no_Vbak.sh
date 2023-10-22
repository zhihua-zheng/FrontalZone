#!/bin/bash -l
### Job Name
#PBS -N run_spinup_no_Vbak
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
#PBS -o Production/run.log
### Join output and error streams into single file
#PBS -j oe
### Email
#PBS -M zhihua@umd.edu
### Send email on abort, begin and end
#PBS -m abe

### Clear and load all the modules needed
module purge
module load julia
module load cuda/11.4.0
# module load netcdf/4.8.1
# module load ncarenv/1.3 gnu/10.1.0 ncarcompilers/0.5.0
# module load openmpi/4.1.1

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run spinup simulation 
proj_dir=$HOME/Projects/TRACE-SEAS/FrontalZone
#--project=<...> activates julia environment
julia --project=$proj_dir ./Production/spinup_no_Vbak.jl

### Overwrite previous log file
preLOG=$proj_dir/Production/spinup_no_Vbak.log
if [ -f "$preLOG" ]; then
    rm -f $preLOG
fi
mv $proj_dir/Production/run.log $preLOG
