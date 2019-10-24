#!/usr/bin/env python

import fire
import os
import yaml
import glob
import math
from analyze import qsub
import pandas as pd

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

def check_starbeast(dir, rep, chains, rerun, ssh=False):
    completed = 0
    for chain in range(1, chains+1):
        chain_dir = os.path.join(dir, "starbeast-chain-{}".format(chain))
        star_pattern = os.path.join(chain_dir, "star-{}-{}.o*".format(rep, chain))
        if check_output(star_pattern, "End likelihood:"):
            completed += 1
        else:
            if rerun:
                with open(os.path.join(chain_dir, "seed.txt")) as fh:
                    seed = fh.readline()
                starbeast_script = [
                   "module load beast \n",
                   "module load beagle \n",
                   "cd {}\n".format(chain_dir),
                   "beast",
                   "-beagle",
                   "-overwrite",  
                   "-seed", seed,
                   os.path.join(dir, "starbeast.xml")]
                qsub(
                    dir=chain_dir,
                    jobname="star-{}-{}".format(rep, chain),
                    walltime="1:00:00",
                    script=" ".join(starbeast_script),
                    ssh=ssh)
    return completed

def check_eco_output(log, nsamples):
    try:
        df = pd.read_csv(log, sep="\t")
        if len(df.index) == nsamples + 1:
            return True
        else:
            print("{} does not have expected length".format(log))
            return False
    except:
        print("unable to open {}".format(log))
        return False

def check_ecoevolity(dir, rep, chains, rerun, nsamples, ssh=False):
    completed = 0
    for chain in range(1, chains+1):
        chain_dir = os.path.join(dir, "ecoevo-chain-{}".format(chain))
        eco_pattern = os.path.join(chain_dir, "ecoevo-{}-{}.o*".format(rep, chain))
        eco_out = os.path.join(chain_dir, "ecoevolity-config-state-run-1.log")
        if check_output(eco_pattern, "Runtime:") and check_eco_output(eco_out, nsamples):
            completed += 1
            pass
        else:
            if rerun:
                try:
                    os.remove(os.path.join(chain_dir, "ecoevolity-config-operator-run-1.log"))
                    os.remove(os.path.join(chain_dir, "ecoevolity-config-state-run-1.log"))
                except:
                    pass
                with open(os.path.join(chain_dir, "seed.txt")) as fh:
                    seed = fh.readline()
                eco_script = [
                    "cd {}\n".format(chain_dir),
                    "ecoevolity-ar",
                    "--seed", seed,
                    "--relax-constant-sites",
                    os.path.join(chain_dir, "ecoevolity-config.yml")]
                qsub(
                    dir=chain_dir,
                    jobname="ecoevo-{}-{}".format(rep, chain),
                    walltime="1:00:00",
                    script=" ".join(eco_script),
                    ssh=ssh)
    return completed

def check(dir, rerun=False, ssh=False, method="all"):
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
        rep_dir = os.path.join(dir, "rep-{}".format(rep))
        if method in ["all", "starbeast"]:
            completed_star += check_starbeast(
              dir=rep_dir, 
              rep=rep, 
              chains=starbeast_chains, 
              rerun=rerun, 
              ssh=ssh)
        if method in ["all", "ecoevolity"]:
            completed_eco += check_ecoevolity(
                dir=rep_dir, 
                rep=rep, 
                chains=ecoevolity_chains, 
                rerun=rerun, 
                nsamples=nsamples,
                ssh=ssh)
    if method in ["all", "ecoevolity"]:
        print("{} out of {} ecoevolity chains complete".format(completed_eco, total_eco))
    if method in ["all", "starbeast"]:
        print("{} out of {} starbeast chains complete".format(completed_star, total_star))
    if rerun:
        print("Incompleted chains restarted")

if __name__ == "__main__":
    fire.Fire(check)
