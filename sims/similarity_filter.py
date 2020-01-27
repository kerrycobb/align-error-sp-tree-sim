#!/usr/bin/env python

import fire
import glob
import yaml
import fire
import subprocess
import os
import dendropy as dp
from dendropy.interop import seqgen
import random
import numpy as np
from tqdm import tqdm


# def check_similarity(align_path, similarity_thresh):
#     align = dp.DnaCharacterMatrix.get(
#         path=align_path,
#         schema="phylip")
#     seg_sites = dp.calculate.popgenstat.num_segregating_sites(align)
#     proportion_similar = (align.max_sequence_size - seg_sites) / align.max_sequence_size
#     if proportion_similar  > similarity_thresh:
#         is_similar = True
#     else:
#         is_similar = False
#     return is_similar 



def run(dir):
    dir = os.path.abspath(dir)
    basename = os.path.basename(dir)
    config_path = os.path.join(dir, "config.yml")
    config = yaml.safe_load(open(config_path))
    eco_config_path = os.path.join(dir, "eco-config.yml")
    eco_config = yaml.safe_load(open(eco_config_path))
    rng = random.Random(config["seed"])
    for rep in range(0, config["nreps"]):
        rep_dir = os.path.join(dir, "rep-{}".format(rep))
        for locus in range(0, config["nloci"]):
            align_path = os.path.join(rep_dir, "alignment-{}.phy".format(locus))
            align = dp.DnaCharacterMatrix.get(
                path=align_path,
                schema="phylip")
            seg_sites = dp.calculate.popgenstat.num_segregating_sites(align)
            print(seg_sites)
        # proportion_similar = (align.max_sequence_size - seg_sites) / align.max_sequence_size




#     """
#     concat: Concatenate loci into single alignment
#     """
#     config = yaml.safe_load(open(config_path))
#     eco_config = yaml.safe_load(open(eco_path))
#     rng = random.Random(seed)

#     # Make output directory
#     out_dir = os.path.join(
#         "out-sp{sp}-gen{gen}-loc{loc}-len{len}".format(
#             sp=config["nspecies"],
#             gen=config["ngenomes"],
#             loc=config["nloci"],
#             len=config["locus_length"]),
#         "seed{}-reps{}".format(seed, config["nreps"]),
#         "singleton-prob-1.0")
#     os.makedirs(out_dir) 

#    # Copy configs into output directory
#     config["seed"] = rng.randint(0, 1000000000)
#     new_config_path = os.path.join(out_dir, "config.yml")
#     yaml.dump(config, open(new_config_path, "w"))
#     new_eco_config_path = os.path.join(out_dir, "eco-config.yml")
#     yaml.dump(eco_config, open(new_eco_config_path, "w"))

#     for rep in tqdm(range(0, config["nreps"])):
#         rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
#         os.mkdir(rep_dir)
    
#         # Simulate species tree
#         sp_tree = simulate_species_tree(
#             n_sp=config["nspecies"],
#             birth_rate=config["birth_rate"],
#             death_rate=0,
#             pop_size_shape=config["pop_size_shape"],
#             pop_size_scale=config["pop_size_scale"],
#             rng=rng)
#         sp_tree.write(
#             path=os.path.join(rep_dir, "species_tree.nex"), 
#             schema="nexus")

#         # Simulate gene trees
#         gene_trees = simulate_gene_trees(
#             sp_tree=sp_tree,
#             nloci=config["nloci"],
#             ngenomes=config["ngenomes"],
#             rng=rng)
#         gene_trees.write(
#             path=os.path.join(rep_dir, "gene_trees.nex"), 
#             schema="nexus")

#         # Simulate alignments
#         simulate_alignment(
#             trees=gene_trees, 
#             length=config["locus_length"], 
#             rng=rng, 
#             dir=rep_dir, 
#             concat=concat)
    
#     print("Simulation Complete")

if __name__ == "__main__":
    fire.Fire(run)
