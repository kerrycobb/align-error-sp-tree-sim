#!/usr/bin/env python

import click
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import re

def set_color(row):
    if row["ess"] < 200 and row["red_fac"] > 1.2:
        return "#f25a64"        # Red
    elif row["ess"] < 200:
        return "#eebd30"        # Yellow
    elif row["red_fac"] > 1.2:
        return "#13b010"        # Green
    else:
        return "#1f77b4"        # Blue

def set_axis(ax, df, lim, interval):
    lower = "{}_lower".format(interval)
    upper = "{}_upper".format(interval)
    for ix, row in df.iterrows():
        ax.errorbar(
            x=row["true"],
            y=row["mean"],
            yerr=[[row["mean"] - row[lower]], [row[upper] - row["mean"]]],
            marker="o",
            markersize=2.5,
            markeredgecolor=row["color"],
            markerfacecolor="none",
            elinewidth=0.5,
            linestyle="",
            ecolor="#1f77b4",
            capsize=0.8
        )
    # ax.set_ylabel("Estimated Value")
    # ax.set_xlabel("True Value")
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
    # ax.set_xlim([0, ax.get_ylim()[1]])
    # ax.set_ylim([0, ax.get_ylim()[1]])
    ax.plot([0, ax.get_xlim()[1]], [0, ax.get_ylim()[1]], linestyle="--", linewidth=.5, color="#1f77b4")

@click.command()
@click.argument("dir", type=click.Path())
@click.argument("alignment", type=click.Path())
def make_plot(dir, alignment, time_lim=0.25, theta_lim=0.006, interval="ci"):
    theta_dfs = []
    theta_cov_dfs = []
    time_dfs = []
    time_cov_dfs = []
    for i in glob.glob(os.path.join(dir, "seed*-reps*")):
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
    cat_theta_cov_df.to_csv(os.path.join(dir,
            "{}-coverage-theta.csv".format(alignment)), index=False)
    
    cat_time_cov_df = pd.concat(time_cov_dfs)\
       .groupby(["method", "interval"])["coverage"]\
       .mean().reset_index()  
    cat_time_cov_df.to_csv(os.path.join(dir, 
            "{}-coverage-time.csv".format(alignment)), index=False)

    # Read in summary csv and output merged summary csv files 
    theta_df = pd.concat(theta_dfs)
    theta_df.to_csv(os.path.join(dir, 
            "{}-summary-theta.csv".format(alignment)), index=False)
    time_df = pd.concat(time_dfs)
    time_df.to_csv(os.path.join(dir, 
            "{}-summary-time.csv".format(alignment)), index=False)

    # Make sure data will fit within defined limits of plots
    upper_int = "{}_lower".format(interval)
    time_mean_greater = time_df[time_df["mean"] > time_lim].count()["mean"]
    time_upper_greater = time_df[time_df[upper_int] > time_lim].count()[upper_int]
    theta_mean_greater = theta_df[theta_df["mean"] > theta_lim].count()["mean"]
    theta_upper_greater = theta_df[theta_df[upper_int] > theta_lim].count()[upper_int]
    greater = [
        (time_mean_greater, "mean time"),
        (time_upper_greater, "time upper interval"),
        (theta_mean_greater, "mean theta"),
        (theta_upper_greater, "theta upper interval")]
    for count, param in greater:
        if count > 0:
            print("Warning: {} {} values are currently outside of plot limit"\
                .format(count, param))
    

    # Assign color to data point 
    theta_df["color"] = theta_df.apply(set_color, axis=1)
    time_df["color"] = time_df.apply(set_color, axis=1)

    # Group data for plotting 
    star_root_theta = theta_df.loc[(theta_df["method"] == "starbeast") & (theta_df["taxon"] == "root")]
    star_theta = theta_df.loc[(theta_df["method"] == "starbeast") & (theta_df["taxon"] != "root")]
    star_time = time_df.loc[(time_df["method"] == "starbeast")]
    eco_root_theta = theta_df.loc[(theta_df["method"] == "ecoevolity") & (theta_df["taxon"] == "root")]
    eco_theta = theta_df.loc[(theta_df["method"] == "ecoevolity") & (theta_df["taxon"] != "root")]
    eco_time = time_df.loc[(time_df["method"] == "ecoevolity")]

    plt.style.use("ggplot")

    plot_params = [
        (star_time, time_lim, "Starbeast Time", "starbeast-time.pdf"),
        (star_root_theta, theta_lim, "Starbeast Root Theta", "starbeast-root-theta.pdf"),
        (star_theta, theta_lim, "Starbeast Theta", "starbeast-theta.pdf"),
        (eco_time, time_lim, "Ecoevolity Time", "ecoevolity-time.pdf"),
        (eco_root_theta, theta_lim, "Ecoevolity Root Theta", "ecoevolity-root-theta.pdf"),
        (eco_theta, theta_lim, "Ecoevolity Theta", "ecoevolity-theta.pdf")]
        
    for i in plot_params:
        f, ax = plt.subplots()    
        set_axis(ax, i[0], i[1], interval) 
        # f.suptitle(i[2])
        plt.savefig(os.path.join(dir, alignment + "-" + i[3]), 
            bbox_inches='tight', pad_inches=0) 

    # f, ax = plt.subplots()
    # set_axis(ax, star_time, time_lim, interval)
    # # f.suptitle("Starbeast Time")
    # plt.savefig(os.path.join(dir, alignment + "-" + "starbeast-time.pdf"))

    # f, ax = plt.subplots()
    # set_axis(ax, star_theta, theta_lim, interval)
    # # f.suptitle("Starbeast Theta")
    # plt.savefig(os.path.join(dir, alignment + "-" + "starbeast-theta.pdf"))

    # f, ax = plt.subplots()
    # set_axis(ax, star_root_theta, theta_lim, interval)
    # # f.suptitle("Starbeast Root Theta")
    # plt.savefig(os.path.join(dir, alignment + "-" + "starbeast-root-theta.pdf"))

    # f, ax = plt.subplots()
    # set_axis(ax, eco_time, time_lim, interval)n
    # # f.suptitle("Ecoevolity Time")
    # plt.savefig(os.path.join(dir, alignment + n "ecoevolity-time.pdf"))

    # f, ax = plt.subplots()
    # set_axis(ax, eco_theta, theta_lim, intervan
    # # f.suptitle("Ecoevolity Theta")
    # plt.savefig(os.path.join(dir, alignment + "-" + "ecoevolity-theta.pdf"))

    # f, ax = plt.subplots()
    # set_axis(ax, eco_root_theta, theta_lim, interval)
    # # f.suptitle("Ecoevolity Root Theta")
    # plt.savefig(os.path.join(dir, alignment + "-" + "ecoevolity-root-theta.pdf"))

if __name__ == "__main__":
    make_plot()
