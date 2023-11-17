#!/usr/bin/env python3
from os import system

casenames = [#"f11_Q000_W000_D000_St0",

             # convection + front
            # "f11_Q005_W000_D000_St0",
            # "f11_Q010_W000_D000_St0",

             # wind + front
            # "f11_Q000_W022_D000_St0",
            # "f11_Q000_W022_D045_St0",
            # "f11_Q000_W022_D090_St0",
            # "f11_Q000_W022_D135_St0",
            # "f11_Q000_W022_D180_St0",
            # "f11_Q000_W022_D225_St0",
            # "f11_Q000_W022_D270_St0",
            # "f11_Q000_W022_D315_St0",

             # wind + wave 
             "h11_Q000_W022_D000_St0",
             "h11_Q000_W022_D000_St1",
            ]


pbs_script = \
"""#!/bin/bash -l
#PBS -A UMCP0020
#PBS -N run_{0}
#PBS -o Logs/run_{0}.log
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
time julia --check-bounds=no --project=$proj_dir Simulations/frontal_zone.jl "{0}"
###2>&1 | tee Logs/{0}.out

### Overwrite previous log file
LOG=$proj_dir/Logs/{0}.log
if [ -f "$LOG" ]; then
    rm -f $LOG
fi
mv $proj_dir/Logs/run_{0}.log $LOG

qstat -f $PBS_JOBID >> Logs/{0}.log
"""

verbose = 1 
pbs_filename = "run_frontal_zone.sh"
for casename in casenames:
    pbs_case = pbs_script.format(casename)

    with open(pbs_filename, "w") as f:
        f.write(pbs_case)

    cmd_run = f"qsub {pbs_filename}"
    if verbose: print(cmd_run)
    system(cmd_run)

    print()
