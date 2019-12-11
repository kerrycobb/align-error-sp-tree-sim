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
        zorder = 200
        if row["color"] != "#1f77b4":
            zorder = 300
        ax.errorbar(
            x=row["true"],
            y=row["mean"],
            yerr=[[row["mean"] - row[lower]], [row[upper] - row["mean"]]],
            marker="o",
            markersize=8.0,
            markeredgewidth = 2.5,
            markeredgecolor=row["color"],
            markerfacecolor="none",
            elinewidth=2.0,
            linestyle="",
            # ecolor="#1f77b4",
            ecolor=row["color"],
            capsize=2.0,
            zorder = zorder,
        )
    # ax.set_ylabel("Estimated Value")
    # ax.set_xlabel("True Value")
    lim_buffer = lim * 0.05
    ax.set_xlim([0 - lim_buffer, lim + lim_buffer])
    ax.set_ylim([0 - lim_buffer, lim + lim_buffer])
    # ax.set_xlim([0, ax.get_ylim()[1]])
    # ax.set_ylim([0, ax.get_ylim()[1]])
    identity_line, = ax.plot(
            [ax.get_xlim()[0], ax.get_xlim()[1]],
            [ax.get_ylim()[0], ax.get_ylim()[1]])
    plt.setp(identity_line,
            color = '0.7',
            linestyle = '-',
            linewidth = 1.5,
            marker = '',
            zorder = 100)
    # Increase size of tick labels
    ax.tick_params(axis = "both", which = "major", labelsize = 16)
    # Only show up to 5 ticks and labels
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))

def get_max_values(dir_name, statistic = "mean"):
    max_time = float("-inf")
    max_theta = float("-inf")
    max_root_theta = float("-inf")
    theta_sum_paths = glob.glob(os.path.join(dir_name, "*-summary-theta.csv"))
    time_sum_paths = glob.glob(os.path.join(dir_name, "*-summary-time.csv"))
    for theta_path in theta_sum_paths:
        theta_df = pd.read_csv(theta_path)

        current_max_root_theta = max(theta_df.loc[theta_df["taxon"] == "root"][statistic])
        current_max_true_root_theta = max(theta_df.loc[theta_df["taxon"] == "root"]["true"])
        current_max_root_theta = max(current_max_root_theta, current_max_true_root_theta)

        current_max_theta = max(theta_df.loc[theta_df["taxon"] != "root"][statistic])
        current_max_true_theta = max(theta_df.loc[theta_df["taxon"] != "root"]["true"])
        current_max_theta = max(current_max_theta, current_max_true_theta)

        if current_max_root_theta > max_root_theta:
            max_root_theta = current_max_root_theta
        if current_max_theta > max_theta:
            max_theta = current_max_theta

    for time_path in time_sum_paths:
        time_df = pd.read_csv(time_path)
        current_max_time = max(time_df[statistic])
        current_max_true_time = max(time_df["true"])
        current_max_time = max(current_max_time, current_max_true_time)
        if current_max_time > max_time:
            max_time = current_max_time
    return max_time, max_root_theta, max_theta


@click.command()
@click.argument("dir_name", type=click.Path())
@click.argument("alignment", type=click.Path())
def make_plot(dir_name, alignment, time_lim=0.25, theta_lim=0.006, interval="ci"):
    # Read in summary csv
    theta_df = pd.read_csv(os.path.join(dir_name, 
            "{}-summary-theta.csv".format(alignment)))
    time_df = pd.read_csv(os.path.join(dir_name, 
            "{}-summary-time.csv".format(alignment)))

    # Make sure data will fit within defined limits of plots
    # upper_int = "{}_lower".format(interval)
    # time_mean_greater = time_df[time_df["mean"] > time_lim].count()["mean"]
    # time_upper_greater = time_df[time_df[upper_int] > time_lim].count()[upper_int]
    # theta_mean_greater = theta_df[theta_df["mean"] > theta_lim].count()["mean"]
    # theta_upper_greater = theta_df[theta_df[upper_int] > theta_lim].count()[upper_int]
    # greater = [
    #     (time_mean_greater, "mean time"),
    #     (time_upper_greater, "time upper interval"),
    #     (theta_mean_greater, "mean theta"),
    #     (theta_upper_greater, "theta upper interval")]
    # for count, param in greater:
    #     if count > 0:
    #         print("Warning: {} {} values are currently outside of plot limit"\
    #             .format(count, param))
    max_time, max_root_theta, max_theta = get_max_values(dir_name,
            statistic = "mean")
    

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
        (star_time, max_time, "Starbeast Time", "starbeast-time.pdf"),
        (star_root_theta, max_root_theta, "Starbeast Root Theta", "starbeast-root-theta.pdf"),
        (star_theta, max_theta, "Starbeast Theta", "starbeast-theta.pdf"),
        (eco_time, max_time, "Ecoevolity Time", "ecoevolity-time.pdf"),
        (eco_root_theta, max_root_theta, "Ecoevolity Root Theta", "ecoevolity-root-theta.pdf"),
        (eco_theta, max_theta, "Ecoevolity Theta", "ecoevolity-theta.pdf")]
        
    for i in plot_params:
        plt.close('all')
        f, ax = plt.subplots()    
        set_axis(ax, i[0], i[1], interval) 
        plt.savefig(os.path.join(dir_name, alignment + "-" + i[3]), 
            bbox_inches='tight', pad_inches=0) 

if __name__ == "__main__":
    make_plot()
