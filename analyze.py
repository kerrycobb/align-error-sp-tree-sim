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

def qsub(dir, jobname, walltime, script):
    qsub = [
        "/cm/shared/apps/torque/6.1.1.1.h3/bin/qsub",
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
        if rand < .5:
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
    cmd = "{qsub} <<EOF\n{script}\nEOF".format(qsub=" ".join(qsub), script=script)
    proc = subprocess.call(cmd, shell=True)

def run_analyses(dir):
    dir = os.path.abspath(dir)
    dir_name = os.path.basename(dir)
    config_path = os.path.join(dir, "config.yml")
    config = yaml.safe_load(open(config_path))
    seed = config["seed"]
    nreps = config["nreps"]
    ecoevolity_chains = config["ecoevolity_chains"]
    starbeast_chains = config["starbeast_chains"]
    rng = random.Random(seed)

    digits = int(math.log10(nreps))+1
    for rep in range(0, nreps):
        rep = str(rep).zfill(digits)
        rep_dir = os.path.join(dir, "replicate-{}".format(rep))
        eco_config_path = os.path.join(rep_dir, "ecoevolity-config.yml")
        align_path = os.path.join(rep_dir, "alignment.nex")
        xml_path = os.path.join(rep_dir, "starbeast.xml")

        # Ecoevolity chains
        for chain in range(1, ecoevolity_chains+1):
            seed = str(rng.randint(1,1000000))
            seed_path = os.path.join(rep_dir, "ecoevo-{}-{}-seed.txt".format(rep, chain))
            with open(seed_path, "w") as fh:
                fh.write(seed)
            eco_script = [
                "cd {}\n".format(rep_dir),
                "ecoevolity-ar",
                "--seed", seed,
                "--relax-constant-sites",
                eco_config_path]
            qsub(
                dir=rep_dir,
                jobname="ecoevo-{}-{}".format(rep, chain),
                walltime="1:00:00",
                script=" ".join(eco_script))

        # Starbeast chains
        for chain in range(1, starbeast_chains+1):
            chain_dir = os.path.join(
                rep_dir,
                "starbeast-chain-{}".format(chain))
            os.mkdir(chain_dir)
            seed = str(rng.randint(1,1000000))
            seed_path = os.path.join(rep_dir, "star-{}-{}-seed.txt".format(rep, chain))
            with open(seed_path, "w") as fh:
                fh.write(seed)
            starbeast_script = [
                "module load beast \n",
                "module load beagle \n",
                "cd {}\n".format(chain_dir),
                "beast",
                "-beagle",
                "-seed", seed,
                xml_path
            ]
            qsub(
                dir=rep_dir,
                jobname="star-{}-{}".format(rep, chain),
                walltime="8:00:00",
                script=" ".join(starbeast_script))

if __name__ == "__main__":
    fire.Fire(run_analyses)
