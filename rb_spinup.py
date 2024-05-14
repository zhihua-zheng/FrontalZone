#!/usr/bin/env python3
from os import system

casenames = [# "",
             # "_NoVbak",
             "_ShortFront"
            ]

with open('Templates/run_spinup.template', 'r') as f:
    pbs_script = f.read()

verbose = 1
pbs_filename = 'run_spinup.sh'
for casename in casenames:
    pbs_case = pbs_script.format(casename)

    with open(pbs_filename, 'w+') as f:
        f.write(pbs_case)

    cmd_run = f'qsub {pbs_filename}'
    if verbose: print(cmd_run)
    system(cmd_run)

    print()
