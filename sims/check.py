#!/usr/bin/env python

import fire
import os
import yaml
import glob
import math
from run import run_starbeast, run_ecoevolity
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
            if lines[pos].startswith("su: cannot set user id:"):
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

def check_starbeast(dir, job, rep, chains, rerun, ssh=False):
    completed = 0
    for chain in range(1, chains+1):
        chain_dir = os.path.join(dir, "starbeast-chain-{}".format(chain))
        # star_pattern = os.path.join(chain_dir, "{}-star-{}-{}.o*".format(job, rep, chain))
        star_pattern = os.path.join(chain_dir, "*.o*")
        if check_output(star_pattern, "End likelihood:"):
            completed += 1
        else:
            if rerun:
                xml_path = os.path.join(dir, "starbeast.xml")
                with open(os.path.join(chain_dir, "seed.txt")) as fh:
                    seed = fh.readline()
                jobname = "{}-star-{}-{}".format(job, rep, chain)
                run_starbeast(chain_dir, seed, jobname, xml_path,
                        rerun=True, ssh=ssh)
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

def check_ecoevolity(dir, job, rep, chains, rerun, nsamples, ssh=False):
    completed = 0
    for chain in range(1, chains+1):
        chain_dir = os.path.join(dir, "ecoevo-chain-{}".format(chain))
        # eco_pattern = os.path.join(chain_dir, "{}-eco-{}-{}.o*".format(job, rep, chain))
        eco_pattern = os.path.join(chain_dir, "*.o*")
        eco_out = os.path.join(chain_dir, "ecoevolity-config-state-run-1.log")
        if check_output(eco_pattern, "Runtime:") and check_eco_output(eco_out, nsamples):
            completed += 1
            pass
        else:
            if rerun:
                with open(os.path.join(chain_dir, "seed.txt")) as fh:
                    seed = fh.readline()
                eco_config_path = os.path.join(chain_dir,
                        "ecoevolity-config.yml")
                jobname = "{}-eco-{}-{}".format(job, rep, chain)
                run_ecoevolity(chain_dir, seed, jobname, eco_config_path,
                        rerun=True, ssh=ssh )
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
    job_basename = "{}-{}".format(os.path.basename(dir), config["locus_length"])
    for rep in range(0, nreps):
        rep_dir = os.path.join(dir, "rep-{}".format(rep))
        if method in ["all", "starbeast"]:
            completed_star += check_starbeast(
              dir=rep_dir,
              job=job_basename,
              rep=rep,
              chains=starbeast_chains,
              rerun=rerun,
              ssh=ssh)
        if method in ["all", "ecoevolity"]:
            completed_eco += check_ecoevolity(
                dir=rep_dir,
                job=job_basename,
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
