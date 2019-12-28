#!/usr/bin/env python3

import fire
import glob
import subprocess
import run
import os
import shutil

def resub():
    outfiles = glob.glob("sim10000-*")
    for file in outfiles:
        with open(file, "r") as fh:
            lines = fh.readlines()
            if lines[-1].startswith("Simulation Complete"):
                pass
            else:
                name = file.split(".")[0]
                seed = name.split("-")[1]
                seed_dir = "out-sp2-gen4-loc10000-len1000/seed{}-reps1/".format(seed) 
                try:
                    shutil.rmtree(seed_dir)
                except:
                    pass
                cmd = "myqsub -N sim10000-{seed} -t 48:00:00 -m 6gb \"./simulate.py {seed} configs/loc10000-len1000.yml configs/eco-config.yml\"".format(seed=seed) 
                subprocess.call(cmd, shell=True)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE, shell=True)
                stdout, stderr = p.communicate()
                if p.returncode == 0:
                    os.remove(file)

if __name__ == "__main__":
    fire.Fire(resub)
