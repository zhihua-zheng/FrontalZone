#!/usr/bin/env python3
from os import system

casenames = [# front
            # "f11_Q000_W000_D000_St0",

             # front + wave
            # "f11_Q000_W000_D000_St1",

             # front + convection
            # "f11_Q005_W000_D000_St0",
            # "f11_Q010_W000_D000_St0",

             # front + wind
            # "f11_Q000_W370_D270_St0",
            # "f11_Q000_W037_D000_St0",
            # "f11_Q000_W037_D090_St0",
            # "f11_Q000_W037_D180_St0",
            # "f11_Q000_W037_D270_St0",
            # "f11_Q000_W010_D000_St0",
            # "f11_Q000_W010_D090_St0",
            # "f11_Q000_W010_D180_St0",
            # "f11_Q000_W010_D270_St0",

             # front + oscillatory wind
            # "f11_Q000_W022_O090_St0",

              # front + wind + wave
            # "f11_Q000_W037_D000_St1",
            # "f11_Q000_W037_D090_St1",
            # "f11_Q000_W037_D180_St1",
            # "f11_Q000_W037_D270_St1",
            # "f11_Q000_W009_D000_St1",
            # "f11_Q000_W009_D090_St1",
            # "f11_Q000_W009_D180_St1",
            # "f11_Q000_W009_D270_St1",

             # short front
            # "s11_Q000_W000_D000_St0",

             # short front + wind
            # "s11_Q000_W370_D000_St0",
            # "s11_Q000_W370_D090_St0",
            # "s11_Q000_W370_D180_St0",
            # "s11_Q000_W370_D270_St0",
            # "s11_M300_Q000_W147_D000_St0",
            # "s11_M300_Q000_W147_D090_St0",
            # "s11_M300_Q000_W147_D180_St0",
            # "s11_M300_Q000_W147_D270_St0",
            # "s11_M030_Q000_W147_D090_St0_Vg0",
             "s11_M030_Q000_W147_D270_St0_Vg0",
             "s11_M030_Q000_W147_D090_St0_Vg1",
             "s11_M030_Q000_W147_D270_St0_Vg1",
            # "s11_Q000_W037_D000_St0",
            # "s11_Q000_W037_D090_St0",
            # "s11_Q000_W037_D180_St0",
            # "s11_Q000_W037_D270_St0",
            # "s11_Q000_W010_D000_St0",
            # "s11_Q000_W010_D090_St0",
            # "s11_Q000_W010_D180_St0",
            # "s11_Q000_W010_D270_St0",

             # wind
            # "n11_Q000_W037_D000_St0",
            # "n11_Q000_W009_D000_St0",

             # wind + wave
            # "n11_Q000_W037_D000_St1",
            # "n11_Q000_W009_D000_St1",

             # wave
            # "n11_Q000_W000_D000_St1",
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
