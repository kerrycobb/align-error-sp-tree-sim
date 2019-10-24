#!/usr/bin/env python3

from scipy.stats import gamma, expon
import dendropy as dp
import matplotlib.pyplot as plt
import numpy as np
import math
import statistics
import os
import fire

def summarize_prior_sample(trees, alpha=5.0, beta=0.0004, lam=20):
    thetas = []
    times = []
    treelist = dp.TreeList.get_from_path(os.path.abspath(trees), 'nexus')
    for tree in treelist[1:]:
        for node in tree:
            thetas.append(float(node.annotations["dmv"].value[0]))
            if node.taxon:
                if node.taxon.label == "T1":
                    times.append(node.edge.length)
    print("Sample Size = {}".format(len(treelist)))

    # Gamma distribution of theta
    mean = alpha * beta
    print("Gamma mean with alpha={alpha} and beta={beta} is {mean}".format(alpha=alpha, beta=beta, mean=mean))
    print("Sample Mean={}".format(statistics.mean(thetas)))
    dist = gamma(alpha, scale=1/beta)
    dist = gamma(alpha, scale=beta)
    x = np.linspace(dist.ppf(0.0001), dist.ppf(0.9999), 1000)
    f, ax = plt.subplots()
    ax.plot (x, dist.pdf(x))
    ax.hist(thetas, 200, density=True)
    ax.set_title("Prior Sample Thetas")
    ax.set_xlabel('θ')
    ax.set_ylabel('Probability')
    plt.savefig("prior-sample-thetas.svg")

    # Exponential distribution of time
    mean = 1/lam 
    print("Exponential mean with lambda={} is {}".format(lam, mean))
    print("Sample Mean={}".format(statistics.mean(times)))
    dist = expon(scale=1/lam) 
    x = np.linspace(dist.ppf(0.0001), dist.ppf(0.9999), 1000)
    f, ax = plt.subplots()
    # ax.plot (x, dist.pdf(x))
    ax.hist(times, 200, density=True)
    ax.set_title("Prior Sample Times")
    ax.set_xlabel('θ')
    ax.set_ylabel('Probability')
    plt.savefig("prior-sample-times.svg")


if __name__ == "__main__":
    fire.Fire(summarize_prior_sample)
