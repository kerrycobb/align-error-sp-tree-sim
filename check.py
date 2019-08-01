#!/usr/bin/env python

import fire
import os
import yaml
import glob
import math
from analyze import qsub

def check_output(pattern, ending):
    complete = False
    error = "{} output file ".format(os.path.basename(pattern))
    outfiles = sorted(glob.glob(pattern), reverse=True)
    if len(outfiles) != 0:
        fh = open(outfiles[0])
        lines = fh.readlines()
        fh.close()
        if len(lines) != 0:
            pos = -1
            if lines[pos].startswith("Warning: Permanently added"):
                pos -= 1
            if lines[pos].startswith(ending):
                complete = True
            else:
                print(error + "is incomplete")
        else:
            print(error + "is empty")
    else:
        print(error + "does not exist")
    return complete

def check(dir, rerun=False):
    dir = os.path.abspath(dir)
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    nreps = config["nreps"]
    ecoevolity_chains = config["ecoevolity_chains"]
    starbeast_chains = config["starbeast_chains"]
    total_eco = nreps * ecoevolity_chains
    total_star = nreps * starbeast_chains
    completed_eco = 0
    completed_star = 0

    digits = int(math.log10(nreps))+1
    for rep in range(0, nreps):
        rep = str(rep).zfill(digits)
        rep_dir = os.path.join(dir, "replicate-{}".format(rep))

        # Check starbeast outputs
        for chain in range(1, starbeast_chains+1):
            star_pattern = os.path.join(rep_dir, "star-{}-{}.o*".format(rep, chain))
            if check_output(star_pattern, "End likelihood:"):
                completed_star += 1
            else:
                if rerun:
                    seed_path = os.path.join(
                        rep_dir, "star-{}-{}-seed.txt".format(rep, chain))
                    with open(seed_path) as fh:
                        seed = fh.readline()
                    chain_dir = os.path.join(
                        rep_dir,
                        "starbeast-chain-{}".format(chain))
                    starbeast_script = [
                       "module load beast \n",
                       "module load beagle \n",
                       "cd {}\n".format(chain_dir),
                       "beast",
                       "-beagle",
                       "-overwrite",
                       "-seed", seed,
                       os.path.join(rep_dir, "starbeast.xml")]
                    qsub(
                        dir=rep_dir,
                        jobname="star-{}-{}".format(rep, chain),
                        walltime="8:00:00",
                        script=" ".join(starbeast_script))

        for chain in range(1, ecoevolity_chains+1):
            # Check Ecoevolity outputs
            eco_pattern = os.path.join(rep_dir, "ecoevo-{}-{}.o*".format(rep, chain))
            if check_output(eco_pattern, "Runtime:"):
                completed_eco += 1
            else:
                if rerun:
                    seed_path = os.path.join(
                        rep_dir, "ecoevo-{}-{}-seed.txt".format(rep, chain))
                    with open(seed_path) as fh:
                        seed = fh.readline()
                    eco_script = [
                        "cd {}\n".format(rep_dir),
                        "ecoevolity-ar",
                        "--seed", seed,
                        "--relax-constant-sites",
                        os.path.join(rep_dir, "ecoevolity-config.yml")]
                    qsub(
                        dir=rep_dir,
                        jobname="ecoevo-{}-{}".format(rep, chain),
                        walltime="1:00:00",
                        script=" ".join(eco_script))

    print("{} out of {} ecoevolity chains complete".format(completed_eco, total_eco))
    print("{} out of {} starbeast chains complete".format(completed_star, total_star))
    if rerun:
        print("Incompleted chains restarted")

if __name__ == "__main__":
    fire.Fire(check)
