#!/usr/bin/env python

import fire
import yaml
import os
import dendropy as dp
import pathlib
import shutil
import pandas as pd
from statistics import mean
import matplotlib.pyplot as plt

def alignment_filter(in_dir, prop, upper_lower):
    """
    upper_lower: < "upper" | "lower" >
        upper = Most variable
        lower = Least variable
    """
    config = yaml.safe_load(open(os.path.join(in_dir, "config.yml")))
    parent_dir = pathlib.Path(in_dir).parent
    new_dir = os.path.join(parent_dir, "filter-{}-{}".format(prop, upper_lower))
    os.mkdir(new_dir)
    shutil.copy(os.path.join(in_dir, "config.yml"), new_dir)
    shutil.copy(os.path.join(in_dir, "eco-config.yml"), new_dir)
    prop_n = int(prop * config["nloci"])
    for rep in range(0, config["nreps"]):
        rep_dir = os.path.join(in_dir, "rep-{}".format(rep))
        new_rep_dir = os.path.join(new_dir, "rep-{}".format(rep))
        os.mkdir(new_rep_dir)
        alignments = []
        seg_sites = []
        loci = []
        for locus in range(0, config["nloci"]):
            align_path = os.path.join(rep_dir, "alignment-{}.phy".format(locus))
            align = dp.DnaCharacterMatrix.get(
                path=align_path,
                schema="phylip")
            segregating = dp.calculate.popgenstat.num_segregating_sites(align)
            seg_sites.append(segregating)
            alignments.append(dict(
                seg_sites=segregating,
                locus=locus,
                alignment=align))
            loci.append(locus)
        alignments.sort(key=lambda dict: dict["seg_sites"])
        if upper_lower == "upper":
            filtered_alignments = alignments[-prop_n:]
        elif upper_lower == "lower":
            filtered_alignments = alignments[:prop_n]
        new_seg_sites = []
        for align in filtered_alignments:
            filt_align_path = os.path.join(new_rep_dir,
                    "alignment-{}.phy".format(align["locus"]))
            align["alignment"].write(path=filt_align_path, schema="phylip")
            new_seg_sites.append(align["seg_sites"])

        # Output stats to text file
        init_str = "Initial: {} loci, mean: {}, min: {}, max: {}\n".format(len(seg_sites),
                mean(seg_sites), min(seg_sites), max(seg_sites))
        filt_str = "Filtered: {} loci, mean: {}, min: {}, max: {}\n".format(len(new_seg_sites),
                mean(new_seg_sites), min(new_seg_sites), max(new_seg_sites))
        with open(os.path.join(new_rep_dir, "output.txt"), "w") as fh:
            fh.write(init_str)
            fh.write(filt_str)

        # Copy species tree
        shutil.copy(os.path.join(rep_dir, "species_tree.nex"),
                os.path.join(new_rep_dir, "species_tree.nex"))

        # Get gene tree for filtered loci and output to new directory
        gene_trees = dp.TreeList.get(path=os.path.join(rep_dir,
                "gene_trees.nex"), schema="nexus")
        new_gene_trees = dp.TreeList()
        for tree in gene_trees:
            if tree.label in loci:
                new_gene_trees.append(tree)
        new_gene_trees.write(path=os.path.join(new_rep_dir,
                "gene_trees.nex"), schema="nexus")
    print("Filtering Complete")

if __name__ == "__main__":
    fire.Fire(alignment_filter)
