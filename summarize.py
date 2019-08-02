#!/usr/bin/env python

import fire
import glob
import os
import yaml
import math
import pycoevolity as pe
import pandas as pd
from check import check_output
import dendropy as dp
from statistics import mean, stdev

def calc_stats(method, rep, true, chains):
    combined = [i for j in chains for i in j]
    lower_quant, upper_quant = pe.stats.get_hpd_interval(combined)
    if len(chains) > 1:
        red_fac = pe.stats.potential_scale_reduction_factor(chains)
    else:
        red_fac = None
    d = dict(
        method=method,
        replicate=rep,
        true=true,
        mean=mean(combined),
        red_fac=red_fac,
        ess=pe.stats.effective_sample_size(combined),
        quant_025=pe.stats.quantile(combined, .025),
        quant_975=pe.stats.quantile(combined, .975),
        hpd_025=lower_quant,
        hpd_975=upper_quant
        )
    return d

def summarize(dir, eco_burnin=501, star_burnin=201):
    dir = os.path.abspath(dir)
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    nreps = config["nreps"]
    eco_chains = config["ecoevolity_chains"]
    star_chains = config["starbeast_chains"]
    digits = int(math.log10(nreps))+1

    theta_dicts = []
    time_dicts = []
    for rep in range(0, nreps):
        rep_id = str(rep).zfill(digits)
        rep_dir = os.path.join(dir, "replicate-{}".format(rep_id))
        true_df = pd.read_csv(os.path.join(rep_dir, "simulated-true-values.txt"), sep='\t')
        true_thetas = dict(
            root=true_df["pop_size_root_sp1"][0],
            sp1=true_df["pop_size_sp1"][0],
            sp2=true_df["pop_size_sp2"][0])

        # Calc and store stats for ecoevolity
        eco_theta_chains = dict(root=[], sp1=[], sp2=[])
        eco_time_chains = []
        for chain in range(1, eco_chains+1):
            eco_pattern = os.path.join(rep_dir, "ecoevo-{}-{}.o*".format(rep_id, chain))
            if check_output(eco_pattern, "Runtime:"):
                eco_log = os.path.join(rep_dir, "ecoevolity-config-state-run-{}.log".format(chain))
                try:
                    est_df = pd.read_csv(eco_log, sep='\t')
                except:
                    quit("Unable to read {}".format(eco_log))
                est_df = est_df.loc[eco_burnin:]
                eco_theta_chains["root"].append(est_df["pop_size_root_sp1"].tolist())
                eco_theta_chains["sp1"].append(est_df["pop_size_sp1"].tolist())
                eco_theta_chains["sp2"].append(est_df["pop_size_sp2"].tolist())
                eco_time_chains.append(est_df["root_height_sp1"].tolist())
        for key in eco_theta_chains:
            theta_dict = calc_stats(
                method="ecoevolity", rep=rep,
                true=true_thetas[key],
                chains=eco_theta_chains[key])
            theta_dict["taxon"] = key
            theta_dicts.append(theta_dict)
        time_dicts.append(calc_stats(
            method="ecoevolity", rep=rep,
            true=true_df["root_height_sp1"][0],
            chains=eco_time_chains))

        # Calc and store stats for starbeast
        star_theta_chains = dict(root=[], sp1=[], sp2=[])
        star_time_chains = []
        for chain in range(1, star_chains+1):
            star_pattern = os.path.join(rep_dir, "star-{}-{}.o*".format(rep_id, chain))
            if check_output(star_pattern, "End likelihood:"):
                trees_path = os.path.join(
                    rep_dir,
                    "starbeast-chain-{}".format(chain),
                    "species.trees")
                trees = dp.TreeList.get(path=trees_path, schema='nexus')
                theta_chain = {key: [] for key in true_thetas}
                time_chain = []
                for tree in trees[star_burnin:]:
                    for node in tree:
                        if node.taxon:
                            theta_chain[node.taxon.label].append(float(node.annotations["dmv"].value[0]))
                            if node.taxon.label == "sp1":
                                time_chain.append(node.edge.length)
                        else:
                            theta_chain["root"].append(float(node.annotations["dmv"].value[0]))
                for key in star_theta_chains:
                    star_theta_chains[key].append(theta_chain[key])
                star_time_chains.append(time_chain)
        for key in star_theta_chains:
            theta_dict = calc_stats(
                method="starbeast", rep=rep,
                true=true_thetas[key],
                chains=star_theta_chains[key])
            theta_dict["taxon"] = key
            theta_dicts.append(theta_dict)
        time_dicts.append(calc_stats(
            method="starbeast", rep=rep,
            true=true_df["root_height_sp1"][0],
            chains=star_time_chains))

    # Output to file
    theta_df = pd.DataFrame(theta_dicts)
    theta_df.to_csv(os.path.join(dir, "summary-theta.csv"), index=False)
    time_df = pd.DataFrame(time_dicts)
    time_df.to_csv(os.path.join(dir, "summary-time.csv"), index=False)
    print("Summary Complete")

if __name__ == "__main__":
    fire.Fire(summarize)
