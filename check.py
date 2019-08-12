#!/usr/bin/env python

import fire
import os
import yaml
import glob
import math
from analyze import qsub
import pandas as pd

def check_output(pattern, ending, show_paths):
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
                if show_paths:
                  print(error + "is incomplete")
        else:
            if show_paths:
                print(error + "is empty")
    else:
        if show_paths:
            print(error + "does not exist")
    return complete

def check_eco_output(log, nsamples, show_paths):
    try:
        df = pd.read_csv(log, sep="\t")
        if len(df.index):
            return True
        else:
            if show_paths:
                print("{} does not have expected length".format(log))
            return False
    except:
        if show_paths:
            print("unable to open {}".format(log))
        return False

def check(dir, rerun=False, show_paths=False):
    dir = os.path.abspath(dir)
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    eco_config = yaml.safe_load(open(os.path.join(dir, "eco-config.yml")))
    nreps = config["nreps"]
    ecoevolity_chains = config["ecoevolity_chains"]
    starbeast_chains = config["starbeast_chains"]
    total_eco = nreps * ecoevolity_chains
    total_star = nreps * starbeast_chains
    completed_eco = 0
    completed_star = 0
    nsamples = int(eco_config["mcmc_settings"]["chain_length"] / eco_config["mcmc_settings"]["sample_frequency"])
    digits = int(math.log10(nreps))+1
    for rep in range(0, nreps):
        rep = str(rep).zfill(digits)
        rep_dir = os.path.join(dir, "replicate-{}".format(rep))

        # Check starbeast outputs
        for chain in range(1, starbeast_chains+1):
            star_pattern = os.path.join(rep_dir, "star-{}-{}.o*".format(rep, chain))
            if check_output(star_pattern, "End likelihood:", show_paths):
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
            eco_out = os.path.join(rep_dir, "ecoevolity-config-state-run-{}.log".format(chain))
            if check_output(eco_pattern, "Runtime:", show_paths) and check_eco_output(eco_out, nsamples, show_paths):
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
