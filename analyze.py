#!/usr/bin/env python

import fire
import glob
import random
import os
import yaml
import time
import subprocess
import random
import math
from shutil import copyfile

def qsub(dir, jobname, walltime, script, ssh=False):
    qsub = [
        "/cm/shared/apps/torque/6.1.1.1.h3/bin/qsub",
        "-m", "n",
        "-N", jobname,
        "-l", "nodes=1:ppn=1",
        "-l", "mem=1G",
        "-l", "walltime={}".format(walltime),
        "-j", "oe",
        "-d", dir,
        "-o", dir,
        ]
    if time.strptime(walltime, "%H:%M:%S") > time.strptime("1:0:0", "%H:%M:%S"):
        rand = random.random()
        # if rand < .5:
        if rand > 1:
            qsub.extend([
            "-q", "general",
            "-W", "group_list=jro0014_lab",
            "-W", "x=FLAGS:ADVRES:jro0014_lab",
            ])
        else:
            qsub.extend([
                "-q", "gen28",
                "-W", "group_list=jro0014_lab",
                "-W", "x=FLAGS:ADVRES:jro0014_s28",
            ])
    ssh_script = ""
    if ssh:
        ssh_script = "ssh kac0070@hopper.auburn.edu "

    cmd = "{ssh}{qsub} <<EOF\n{script}\nEOF".format(
        ssh=ssh_script, qsub=" ".join(qsub), script=script)
    proc = subprocess.call(cmd, shell=True)

def run_starbeast(dir, rep, nchains, rng, ssh=False):
    for chain in range(1, nchains+1):
        xml = os.path.join(dir, "starbeast.xml")
        chain_dir = os.path.join(dir, "starbeast-chain-{}".format(chain))
        os.mkdir(chain_dir)
        seed = str(rng.randint(1,1000000))
        with open(os.path.join(chain_dir, "seed.txt"), "w") as fh:
            fh.write(seed)
        starbeast_script = [
            "module load beast \n",
            "module load beagle \n",
            "cd {}\n".format(chain_dir),
            "beast",
            "-beagle",
            "-seed", seed,
            xml]
        qsub(
            dir=chain_dir,
            jobname="star-{}-{}".format(rep, chain),
            walltime="8:00:00",
            script=" ".join(starbeast_script),
            ssh=ssh)
 
def run_ecoevolity(dir, rep, nchains, rng, ssh=False):
    for chain in range(1, nchains+1):
        chain_dir = os.path.join(dir, "ecoevo-chain-{}".format(chain))
        os.mkdir(chain_dir)
        config = os.path.join(dir, "ecoevolity-config.yml")
        new_config = os.path.join(chain_dir, "ecoevolity-config.yml")
        copyfile(config, new_config)
        seed = str(rng.randint(1,1000000))
        with open(os.path.join(chain_dir, "seed.txt"), "w") as fh:
            fh.write(seed)
        eco_script = [
            "cd {}\n".format(chain_dir),
            "ecoevolity-ar",
            "--seed", seed,
            "--relax-constant-sites",
            new_config]
        qsub(
            dir=chain_dir,
            jobname="ecoevo-{}-{}".format(rep, chain),
            walltime="1:00:00",
            script=" ".join(eco_script),
            ssh=ssh)

def run_analyses(dir, ssh=False):
    dir = os.path.abspath(dir)
    config_path = os.path.join(dir, "config.yml")
    config = yaml.safe_load(open(config_path))
    seed = config["seed"]
    nreps = config["nreps"]
    eco_chains = config["ecoevolity_chains"]
    star_chains = config["starbeast_chains"]
    rng = random.Random(seed)

    digits = int(math.log10(nreps))+1
    for rep in range(0, nreps):
        rep = str(rep).zfill(digits)
        rep_dir = os.path.join(dir, "replicate-{}".format(rep))
        run_ecoevolity(rep_dir, rep, eco_chains, rng, ssh=ssh)
        run_starbeast(rep_dir, rep, star_chains, rng, ssh=ssh)
        

if __name__ == "__main__":
    fire.Fire(run_analyses)
