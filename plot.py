#!/usr/bin/env python

import click
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def set_color(row):
    if row["ess"] < 200 and row["red_fac"] > 1.2:
        return "#f25a64"        # Red
    elif row["ess"] < 200:
        return "#eebd30"        # Yellow
    elif row["red_fac"] > 1.2:
        return "#13b010"        # Green
    else:
        return "#1f77b4"        # Blue

def set_axis(ax, df):
    for ix, row in df.iterrows():
        ax.errorbar(
            x=row["true"],
            y=row["mean"],
            yerr=[[row["mean"] - row["hpd_025"]], [row["hpd_975"] - row["mean"]]],
            marker="o",
            markersize=2.5,
            markeredgecolor=row["color"],
            markerfacecolor="none",
            elinewidth=0.5,
            linestyle="",
            ecolor="#1f77b4",
            capsize=0.8
        )
    ax.set_ylabel("Estimated Value")
    ax.set_xlabel("True Value")
    ax.set_xlim([0, ax.get_ylim()[1]])
    ax.set_ylim([0, ax.get_ylim()[1]])
    ax.plot([0, ax.get_xlim()[1]], [0, ax.get_ylim()[1]], linestyle="--", linewidth=.5, color="#1f77b4")

@click.command()
@click.argument("dir", type=click.Path())
@click.argument("alignment", type=click.Path())
def make_plot(dir, alignment):
    theta_dfs = []
    time_dfs = []
    for i in glob.glob(os.path.join(dir, "seed*-reps*")):
        try:
            theta_dfs.append(pd.read_csv(os.path.join(i, alignment, "summary-theta.csv")))
            time_dfs.append(pd.read_csv(os.path.join(i, alignment, "summary-time.csv")))
        except:
            print("Warning: no summary .csv file in {}".format(i +"/"+alignment))
            pass
    theta_df = pd.concat(theta_dfs)
    time_df = pd.concat(time_dfs)

    theta_df["color"] = theta_df.apply(set_color, axis=1)
    time_df["color"] = time_df.apply(set_color, axis=1)

    star_root_theta = theta_df.loc[(theta_df["method"] == "starbeast") & (theta_df["taxon"] == "root")]
    star_theta = theta_df.loc[(theta_df["method"] == "starbeast") & (theta_df["taxon"] != "root")]
    star_time = time_df.loc[(time_df["method"] == "starbeast")]

    eco_root_theta = theta_df.loc[(theta_df["method"] == "ecoevolity") & (theta_df["taxon"] == "root")]
    eco_theta = theta_df.loc[(theta_df["method"] == "ecoevolity") & (theta_df["taxon"] != "root")]
    eco_time = time_df.loc[(time_df["method"] == "ecoevolity")]


    plt.style.use("ggplot")

    f, ax = plt.subplots()
    set_axis(ax, star_time)
    f.suptitle("Starbeast Time")
    plt.savefig(os.path.join(dir, alignment + "-" + "starbeast-time.svg"))

    f, ax = plt.subplots()
    set_axis(ax, star_theta)
    f.suptitle("Starbeast Theta")
    plt.savefig(os.path.join(dir, alignment + "-" + "starbeast-theta.svg"))

    f, ax = plt.subplots()
    set_axis(ax, star_root_theta)
    f.suptitle("Starbeast Root Theta")
    plt.savefig(os.path.join(dir, alignment + "-" + "starbeast-root-theta.svg"))

    f, ax = plt.subplots()
    set_axis(ax, eco_time)
    f.suptitle("Ecoevolity Time")
    plt.savefig(os.path.join(dir, alignment + "-" + "ecoevolity-time.svg"))

    f, ax = plt.subplots()
    set_axis(ax, eco_theta)
    f.suptitle("Ecoevolity Theta")
    plt.savefig(os.path.join(dir, alignment + "-" + "ecoevolity-theta.svg"))

    f, ax = plt.subplots()
    set_axis(ax, eco_root_theta)
    f.suptitle("Ecoevolity Root Theta")
    plt.savefig(os.path.join(dir, alignment + "-" + "ecoevolity-root-theta.svg"))

if __name__ == "__main__":
    make_plot()
