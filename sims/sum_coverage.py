#!/usr/bin/env python

import click
import os
import glob
import pandas as pd
import re

@click.command()
@click.argument("dir_name", type=click.Path())
@click.argument("alignment", type=click.Path())
def summarize_coverage(dir_name, alignment, interval="ci"):
    theta_dfs = []
    theta_cov_dfs = []
    time_dfs = []
    time_cov_dfs = []
    for i in glob.glob(os.path.join(dir_name, "seed*-reps*")):
        seed = re.split(r'/|-', i)[5] 
        try:
            # Get theta summary and coverage
            theta_df = pd.read_csv(os.path.join(i, alignment, 
                    "summary-theta.csv"))
            theta_df["seed"] = seed
            theta_dfs.append(theta_df)
            theta_cov_df = pd.read_csv(os.path.join(i, alignment,
                    "coverage-theta.csv"))
            theta_cov_dfs.append(theta_cov_df)
            # Get time summary and coverage
            time_df = pd.read_csv(os.path.join(i, alignment, 
                    "summary-time.csv"))
            time_df["seed"] = seed
            time_dfs.append(time_df)
            time_cov_df = pd.read_csv(os.path.join(i, alignment, 
                    "coverage-time.csv"))
            time_cov_dfs.append(time_cov_df)
        except:
            print("Warning: missing summary file in {}".format(i +"/"+alignment))

    # Sum coverage from dataframes and output to csv 
    cat_theta_cov_df = pd.concat(theta_cov_dfs)\
       .groupby(["method", "interval", "taxon"])["coverage"]\
       .mean().reset_index()
    cat_theta_cov_df.to_csv(os.path.join(dir_name,
            "{}-coverage-theta.csv".format(alignment)), index=False)
    
    cat_time_cov_df = pd.concat(time_cov_dfs)\
       .groupby(["method", "interval"])["coverage"]\
       .mean().reset_index()  
    cat_time_cov_df.to_csv(os.path.join(dir_name, 
            "{}-coverage-time.csv".format(alignment)), index=False)

    # Read in summary csv and output merged summary csv files 
    theta_df = pd.concat(theta_dfs)
    theta_df.to_csv(os.path.join(dir_name, 
            "{}-summary-theta.csv".format(alignment)), index=False)
    time_df = pd.concat(time_dfs)
    time_df.to_csv(os.path.join(dir_name, 
            "{}-summary-time.csv".format(alignment)), index=False)

if __name__ == "__main__":
    summarize_coverage()
