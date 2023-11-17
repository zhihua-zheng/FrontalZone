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
            # "n11_Q000_W022_D000_St0",
             "n11_Q000_W022_D000_St1",
            ]

with open("run_frontal_zone.template", "r") as f:
    pbs_script = f.read()

verbose = 1 
pbs_filename = "run_frontal_zone.sh"
for casename in casenames:
    pbs_case = pbs_script.format(casename, "")

    with open(pbs_filename, "w+") as f:
        f.write(pbs_case)

    cmd_run = f"qsub {pbs_filename}"
    if verbose: print(cmd_run)
    system(cmd_run)

    print()
