#!/usr/bin/env python

import click
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import re
import numpy as np
import pycoevolity as pe

def set_color(row):
    if row["ess"] < 200 and row["red_fac"] > 1.2:
        return "#f25a64"        # Red
    elif row["ess"] < 200:
        return "#eebd30"        # Yellow
    elif row["red_fac"] > 1.2:
        return "#13b010"        # Green
    else:
        return "#1f77b4"        # Blue

# def check_interval(df, interval):
#     """
#     Make sure data will fit within defined limits of plots
#     """
#     upper_int = "{}_lower".format(interval)
#     time_mean_greater = df[df["mean"] > time_lim].count()["mean"]
#     time_upper_greater = df[df[upper_int] > time_lim].count()[upper_int]
#     theta_mean_greater = theta_df[theta_df["mean"] > theta_lim].count()["mean"]
#     theta_upper_greater = theta_df[theta_df[upper_int] > theta_lim].count()[upper_int]
#     greater = [
#         (time_mean_greater, "mean time"),
#         (time_upper_greater, "time upper interval"),
#         (theta_mean_greater, "mean theta"),
#         (theta_upper_greater, "theta upper interval")]
#     for count, param in greater:
#         if count > 0:
#             print("Warning: {} {} values are currently outside of plot limit"\
#                 .format(count, param))

def add_subplot(ax, df, lim, interval="ci"):
    lower = "{}_lower".format(interval)
    upper = "{}_upper".format(interval)
    # check_interval(df, interval)
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
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
    ax.set_aspect('equal')
    # ax.set_aspect("equal", adjustable="box")
    # ax.set_aspect(1.0, adjustable="box")
    ax.plot([0, ax.get_xlim()[1]], [0, ax.get_ylim()[1]], linestyle="--", linewidth=.5, color="#1f77b4")


len1000 = "out-sp2-gen4-loc100-len1000"
len500 = "out-sp2-gen4-loc200-len500"
len250 = "out-sp2-gen4-loc400-len250"

single_1_theta = "singleton-prob-1.0-summary-theta.csv"
single_1_time = "singleton-prob-1.0-summary-time.csv"

single_08_theta = "singleton-prob-0.8-summary-theta.csv"
single_08_time = "singleton-prob-0.8-summary-time.csv"

single_06_theta = "singleton-prob-0.6-summary-theta.csv"
single_06_time = "singleton-prob-0.6-summary-time.csv"

het_08_theta = "het-prob-0.8-summary-theta.csv"
het_08_time = "het-prob-0.8-summary-time.csv"

het_08_theta = "het-prob-0.6-summary-theta.csv"
het_08_time = "het-prob-0.6-summary-time.csv"

lengths = [len1000, len500, len250]
methods = ["starbeast", "ecoevolity"]

time_lim = 0.25
time_step = 0.05
theta_lim = 0.006
theta_step = 0.002
interval = "ci"  

# plt.style.use("ggplot")
# blue = "#1f77b4"

# fig = plt.figure(figsize=(4,6))
# gs = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[1,1], 
#         wspace=0, hspace=0)



### Time Singleton 1.0 
fig, axes = plt.subplots(nrows=3, ncols=2,figsize=(4,6))
# fig.suptitle("Divergence Time")
fig.subplots_adjust(wspace=0.01, hspace=0)

for i, length in enumerate(lengths): 
    df = pd.read_csv(os.path.join(length, single_1_time))
    df["color"] = df.apply(set_color, axis=1)
    for j, method in enumerate(methods):
        # ax = plt.subplot(gs[i,j])
        ax = axes[i,j]
        add_subplot(
            ax=ax, 
            df=df.loc[(df["method"] == method)],
            lim=time_lim)
        # ax.set_aspect(adjustable="box", aspect=1.0)
        ax.set_aspect(aspect=1.0)
        # ax.set_xticks(np.arange(0, time_lim + time_step, step=time_step))
        # ax.set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))
        # ax.tick_params(labelleft=False, labelbottom=False)


        # axes[i,j].set(adjustable='box-forced', aspect='equal')
#         # axes[i,j].set_aspect('equal')
#         # axes[i,j].set_aspect(aspect=1.0)
#         # , adjustable="box")
#         # axes[i,j].set(adjustable='box-forced', aspect='equal')
#         # axes[i,j].set_xticks(np.arange(0, time_lim + time_step, step=time_step))
#         # axes[i,j].set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))
#         # if i == 2 and j == 0:
#         #     pass
#         # else:
#         #     axes[i,j].tick_params(labelleft=False, labelbottom=False)

# Set labels
# gs[0,0].set_title("Starbeast")
# gs[0,1].set_title("Ecoevolity")
# axes[0,0].set_title("Starbeast")
# axes[0,1].set_title("Ecoevolity")

# Axis labels
# axes[2,0].set_ylabel("Estimated Divergence Time", fontsize=10)
# axes[2,0].set_xlabel('True Divergence Time', fontsize=10)

# # Row labels
# axes[0,1].text(1.05, 0.5,"1000 BP", verticalalignment="center",
#         rotation="horizontal", transform=axes[0,1].transAxes, fontsize=10)
# axes[1,1].text(1.05, 0.5,"500 BP", verticalalignment="center",
#         rotation="horizontal", transform=axes[1,1].transAxes, fontsize=10)
# axes[2,1].text(1.05, 0.5,"250 BP", verticalalignment="center",
#         rotation="horizontal", transform=axes[2,1].transAxes, fontsize=10)



# # Get rid of internal y ticks
# axes[0,1].set_yticks([])
# axes[1,1].set_yticks([])
# axes[2,1].set_yticks([])
# # X ticks
# axes[1,0].set_xticks(np.arange(0, time_lim + time_step,aspect=1.0 step=time_step))
# axes[1,1].set_xticks(np.arange(time_step, time_lim + time_step, step=time_step))
# # Y ticks
# axes[0,0].set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))
# axes[1,0].set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))
# axes[2,0].set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))
# axes[0,1].set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))
# axes[1,1].set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))
# axes[2,1].set_yticks(np.arange(time_step, time_lim + time_step, step=time_step))

# axes[0,0].tick_params(labelleft=False)
# axes[1,0].tick_params(labelleft=False)

# axes[2,1].xaxis.set_ticklabels([])

plt.savefig("divergence-time-singleton-1.0-plot.pdf")