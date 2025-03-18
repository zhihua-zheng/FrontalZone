#!/usr/bin/env python3
import subprocess

# Vg0 doesn't work as there is no pressure gradient from the background buoyancy
casenames = [# front
            # "f11_Q000_W000_D000_St0",
            # "f11_M009_Q000_W148_D180_St0",

             # front + oscillatory wind
            # "f11_Q000_W022_O090_St0",

             # short front + wind
            # "s11_M036_Q000_W037_D270_St0",
            # "s11_M036_Q-135_W037_D270_St0",
            # "s11_M036_Q000_W444_D090_St0",

            # "s11_M018_Q000_W148_D090_St0",
            # "s11_M018_Q000_W444_D090_St0",

            # "s11_M009_Q000_W444_D090_St0",

            # "s11_M009_Q000_W148_D000_St0",
            # "s11_M009_Q000_W148_D090_St0",
            # "s11_M009_Q000_W148_D180_St0",
            # "s11_M009_Q000_W148_D270_St0",
            # "s11_M009_Q135_W148_D090_St0",

            # "s11_M003_Q000_W444_D270_St0",

            # "s22_M036_Q000_W037_D270_St0",
            # "s22_M003_Q000_W444_D270_St0",
            # "s22_M009_Q000_W148_D000_St0",
            # "s22_M009_Q000_W148_D180_St0",
            # "s22_M009_Q000_W148_D090_St0",
            # "s22_M009_Q135_W148_D090_St0",

             # wind
            # "s11_M000_Q000_W148_D180_St0",
             "s11_M000_Q000_W444_D180_St0",
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
    if verbose: print(cmd_run)
    result = subprocess.run(cmd_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout.decode())

    print()
