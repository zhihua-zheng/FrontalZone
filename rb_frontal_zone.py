#!/usr/bin/env python3
#from os import system
import subprocess

# Vg0 doesn't work as there is no pressure gradient from the background buoyancy
casenames = [# front
            # "f11_Q000_W000_D000_St0",

             # front + wave
            # "f11_Q000_W000_D000_St1",

             # front + convection
            # "f11_Q005_W000_D000_St0",
            # "f11_Q010_W000_D000_St0",

             # front + wind
            # "f11_Q000_W147_D270_St0",
            # "f11_Q000_W037_D000_St0",
            # "f11_Q000_W037_D090_St0",
            # "f11_Q000_W037_D180_St0",
            # "f11_Q000_W037_D270_St0",
            # "f11_Q000_W009_D000_St0",
            # "f11_Q000_W009_D090_St0",
            # "f11_Q000_W009_D180_St0",
            # "f11_Q000_W009_D270_St0",

             # front + oscillatory wind
            # "f11_Q000_W022_O090_St0",

              # front + wind + wave
            # "f11_Q000_W037_D000_St1",
            # "f11_Q000_W037_D090_St1",
            # "f11_Q000_W037_D180_St1",
            # "f11_Q000_W037_D270_St1",

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

             "s11_M030_Q000_W147_D090_St0",
             "s11_M030_Q000_W147_D270_St0",

            # "s11_M030_Q000_W037_D000_St0",
            # "s11_M030_Q000_W037_D090_St0",
            # "s11_M030_Q000_W037_D180_St0",
            # "s11_M030_Q000_W037_D270_St0",

             # wind
            # "n11_Q000_W037_D000_St0",
            # "n11_Q000_W009_D000_St0",

             # wind + wave
            # "n11_Q000_W037_D000_St1",
            # "n11_Q000_W009_D000_St1",

             # wave
            # "n11_Q000_W000_D000_St1",
            ]

with open('Templates/run_frontal_zone.template', 'r') as f:
    main_script = f.read()
with open('Templates/run_KE_budgets.template', 'r') as f:
    post_script = f.read()

verbose = 1
main_filename = 'run_frontal_zone.sh'
post_filename = 'run_KE_budgets.sh'
for casename in casenames:
    pbs_main = main_script.format(casename, '')#--nTf 6.1
    pbs_post = post_script.format(casename, '')

    with open(main_filename, 'w+') as f:
        f.write(pbs_main)
    with open(post_filename, 'w+') as f:
        f.write(pbs_post)

    cmd_run = f'JOBID1=$(qsub -h {main_filename}); qsub -W depend=afterok:$JOBID1 {post_filename}; qrls $JOBID1'
    #cmd_run = f'qsub {pbs_filename}'
    if verbose: print(cmd_run)
    #system(cmd_run)
    result = subprocess.run(cmd_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout.decode())

    print()
