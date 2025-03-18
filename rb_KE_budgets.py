#!/usr/bin/env python3
from os import system

casenames = [# front
            # "f11_Q000_W000_D000_St0",

             # front + oscillatory wind
            # "f11_Q000_W022_O090_St0",

             # short front + wind
             "s11_M036_Q000_W037_D270_St0",
             #"s11_M036_Q-135_W037_D270_St0",
             #"s11_M036_Q000_W444_D090_St0",

             #"s11_M018_Q000_W148_D090_St0",
             #"s11_M018_Q000_W444_D090_St0",

             #"s11_M009_Q000_W444_D090_St0",

             "s11_M009_Q000_W148_D000_St0",
             "s11_M009_Q000_W148_D090_St0",
             "s11_M009_Q000_W148_D180_St0",
             "s11_M009_Q000_W148_D270_St0",
             "s11_M009_Q135_W148_D090_St0",

             "s11_M003_Q000_W444_D270_St0",

             #"s22_M036_Q000_W037_D270_St0",
             #"s22_M003_Q000_W444_D270_St0",
             #"s22_M009_Q000_W148_D000_St0",
             #"s22_M009_Q000_W148_D180_St0",
             #"s22_M009_Q000_W148_D090_St0",
             #"s22_M009_Q135_W148_D090_St0",

             # wind
             #"s11_M000_Q000_W148_D180_St0",
             #"s11_M000_Q000_W444_D180_St0",
            ]

with open('Templates/run_KE_budgets.template', 'r') as f:
    pbs_script = f.read()

verbose = 1
pbs_filename = 'run_KE_budgets.sh'
for casename in casenames:
    pbs_case = pbs_script.format(casename)

    with open(pbs_filename, 'w+') as f:
        f.write(pbs_case)

    cmd_run = f'qsub {pbs_filename}'
    if verbose: print(cmd_run)
    system(cmd_run)
    print()
