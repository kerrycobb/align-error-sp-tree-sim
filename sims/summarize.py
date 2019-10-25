#!/usr/bin/env python

import fire
import glob
import os
import yaml
import math
import pycoevolity as pe
import pandas as pd
from check import check_output, check_eco_output
import dendropy as dp
from statistics import mean, stdev
from tqdm import tqdm
from decimal import Decimal as D
from decimal import getcontext

def calc_stats(method, rep, true, chains, interval=0.95):
    combined = [i for j in chains for i in j]
    
    # Determine upper and lower quantiles
    getcontext().prec = 3
    quant = (D(1.0) - D(interval)) / 2
    lower_quant = 0 + quant
    upper_quant = 1 - quant

    # Calculate credible interval and hpd
    lower_hpd, upper_hpd = pe.stats.get_hpd_interval(combined, 
            interval_prob=interval)
    lower_ci = pe.stats.quantile(combined, lower_quant)
    upper_ci = pe.stats.quantile(combined, upper_quant)

    # Check if true value within hpd
    if lower_hpd < true < upper_hpd:
        hpd_covered = True 
    else:
        hpd_covered = False

    # Check if true value within credible interval
    if lower_ci < true < upper_ci:
        ci_covered = True
    else: 
        ci_covered = False

    # Calculate potential scale reduction factor
    if len(chains) > 1:
        red_fac = pe.stats.potential_scale_reduction_factor(chains)
    else:
        red_fac = None
    
    # Add stats to dictionary
    d = dict(
        method=method,
        replicate=rep,
        true=true,
        mean=mean(combined),
        red_fac=red_fac,
        ess=pe.stats.effective_sample_size(combined),
        interval=interval,
        ci_lower=lower_ci,
        ci_upper=upper_ci,
        ci_covered=ci_covered,
        hpd_lower=lower_hpd,
        hpd_upper=upper_hpd,
        hpd_covered=hpd_covered)

    return d

def summarize(dir, eco_burnin=501, star_burnin=201, method="all", interval=0.95):
    dir = os.path.abspath(dir)
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    nreps = config["nreps"]
    eco_chains = config["ecoevolity_chains"]
    star_chains = config["starbeast_chains"]
    theta_dicts = []
    time_dicts = []
    for rep in tqdm(range(0, nreps)):
        rep_dir = os.path.join(dir, "rep-{}".format(rep))

        # Get true simulated values
        sim_tree_path = os.path.join(rep_dir, "species_tree.nex")
        sim_tree = dp.Tree.get(path=sim_tree_path, schema="nexus")
        true_thetas = {}
        true_times = {}
        for node in sim_tree:
            if node.taxon:
                label = node.taxon.label
            else:
                label = "root" 
            true_thetas[label] = float(node.annotations["pop_size"].value)
            true_times[label] = float(node.edge.length)

        if method in ["all", "ecoevolity"]:
            # Calculate and store stats for ecoevolity
            eco_theta_chains = dict(root=[], T1=[], T2=[])
            eco_time_chains = []
            for chain in range(1, eco_chains+1):
                chain_dir = os.path.join(rep_dir, "ecoevo-chain-{}".format(chain))
                # eco_pattern = os.path.join(chain_dir, "ecoevo-{}-{}.o*".format(rep, chain))
                eco_log = os.path.join(chain_dir, "ecoevolity-config-state-run-1.log")
                est_df = pd.read_csv(eco_log, sep='\t').loc[eco_burnin:]
                eco_theta_chains["root"].append(est_df["pop_size_root_T1"].tolist())
                eco_theta_chains["T1"].append(est_df["pop_size_T1"].tolist())
                eco_theta_chains["T2"].append(est_df["pop_size_T2"].tolist())
                eco_time_chains.append(est_df["root_height_T1"].tolist())
            for key in eco_theta_chains:
                theta_dict = calc_stats(
                    method="ecoevolity",
                    rep=rep,
                    true=true_thetas[key],
                    chains=eco_theta_chains[key],
                    interval=interval)
                theta_dict["taxon"] = key
                theta_dicts.append(theta_dict)
            time_dicts.append(calc_stats(
                method="ecoevolity",
                rep=rep,
                true=true_times["T1"],
                chains=eco_time_chains, 
                interval=interval))

        if method in ["all", "starbeast"]:
            # Calculate and store stats for starbeast
            star_theta_chains = dict(root=[], T1=[], T2=[])
            star_time_chains = []
            for chain in range(1, star_chains+1):
                chain_dir = os.path.join(
                        rep_dir, "starbeast-chain-{}".format(chain))
                trees_path = os.path.join(chain_dir, "species.trees")
                trees = dp.TreeList.get(path=trees_path, schema='nexus')
                theta_chain = {key: [] for key in true_thetas}
                time_chain = []
                for tree in trees[star_burnin:]:
                    for node in tree:
                        if node.taxon:
                            theta_chain[node.taxon.label].append(
                                    float(node.annotations["dmv"].value[0]))
                            if node.taxon.label == "T1":
                                time_chain.append(node.edge.length)
                        else:
                            theta_chain["root"].append(float(
                                    node.annotations["dmv"].value[0]))
                for key in star_theta_chains:
                    star_theta_chains[key].append(theta_chain[key])
                star_time_chains.append(time_chain)
            for key in star_theta_chains:
                theta_dict = calc_stats(
                    method="starbeast",
                    rep=rep,
                    true=true_thetas[key],
                    chains=star_theta_chains[key],
                    interval=interval)
                theta_dict["taxon"] = key
                theta_dicts.append(theta_dict)
            time_dicts.append(calc_stats(
                method="starbeast", 
                rep=rep,
                true=true_times["T1"], 
                chains=star_time_chains,
                interval=interval))
    
    # Output stats to file
    theta_df = pd.DataFrame(theta_dicts)
    theta_df.to_csv(os.path.join(dir, "summary-theta.csv"), index=False)
    time_df = pd.DataFrame(time_dicts)
    time_df.to_csv(os.path.join(dir, "summary-time.csv"), index=False)

    # Here for testing, avoid recomputing all stats above
    # theta_df = pd.read_csv(os.path.join(dir, "summary-theta.csv"))
    # time_df = pd.read_csv(os.path.join(dir, "summary-time.csv"))

    # Get theta coverage
    theta_coverage = []
    theta_coverage_strs = [] 
    for method, group in theta_df.groupby(["method"]):
        for taxon, g in group.groupby(["taxon"]):
            hpd_cover = g["hpd_covered"].sum() / nreps
            ci_cover = g["ci_covered"].sum() / nreps
            theta_coverage.append(dict(
                interval="hpd",
                method=method,
                taxon=taxon,
                coverage=hpd_cover))
            theta_coverage.append(dict(
                interval="ci",
                method=method,
                taxon=taxon,
                coverage=ci_cover))
            theta_coverage_strs.append(
                "Theta-{}-{}-hpd= {}".format(method, taxon, hpd_cover)) 
            theta_coverage_strs.append(
                "Theta-{}-{}-ci= {}".format(method, taxon, ci_cover)) 

    # Get time coverage
    time_coverage = []
    time_coverage_strs = [] 
    for method, group in time_df.groupby(["method"]):
        hpd_cover = g["hpd_covered"].sum() / nreps
        ci_cover = g["ci_covered"].sum() / nreps
        time_coverage.append(dict(
            interval="hpd",
            method=method,
            coverage=hpd_cover))
        time_coverage.append(dict(
            interval="ci",
            method=method,
            coverage=ci_cover))
        time_coverage_strs.append(
            "Time-{}-hpd= {}".format(method, hpd_cover)) 
        time_coverage_strs.append(
            "Time-{}-ci= {}".format(method, ci_cover)) 

    # Store coverage in csv
    theta_coverage_df = pd.DataFrame(theta_coverage)
    theta_coverage_df.to_csv(os.path.join(dir, "coverage-theta.csv"), index=False)
    time_coverage_df = pd.DataFrame(time_coverage)
    time_coverage_df.to_csv(os.path.join(dir, "coverage-time.csv"), index=False)

    # Print coverage and write strings to text file
    for i in theta_coverage_strs:
        print(i)
    with open(os.path.join(dir, "coverage-theta.txt"), "w") as fh:
        fh.write('\n'.join(theta_coverage_strs))

    for i in time_coverage_strs:
        print(i)
    with open(os.path.join(dir, "coverage-time.txt"), "w") as fh:
        fh.write('\n'.join(time_coverage_strs))


    print("Summary Complete")

if __name__ == "__main__":
    fire.Fire(summarize)
