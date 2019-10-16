#!/usr/bin/env python

import fire
import glob
import yaml
import fire
import subprocess
import os
import jinja2 as jj
import dendropy as dp
import math
import random
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

def simulate_species_tree(n_sp, birth_rate, death_rate,
        pop_size_shape, pop_size_scale, rng):
    sp_tree = dp.simulate.treesim.birth_death_tree(
        birth_rate=birth_rate,
        death_rate=death_rate,
        num_extant_tips=n_sp,
        gsa_ntax=n_sp + 1,
        rng=rng)
    sp_tree.seed_node.edge.length = 0.0
    sp_tree.calc_node_ages()
    for node in sp_tree:
        node.annotations.drop()
    for node in sp_tree:
        node.edge.pop_size = rng.gammavariate(pop_size_shape, pop_size_scale)
        node.annotations['pop_size'] = node.edge.pop_size
    return sp_tree

def simulate_gene_trees(sp_tree, n_tips_per_sp, n_loci, rng):
    gene_to_species_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping(
        sp_tree.taxon_namespace,
        n_tips_per_sp)
    gene_trees = dp.TreeList()
    for i in range(n_loci):
        gene_tree = dp.simulate.treesim.contained_coalescent_tree(
                containing_tree = sp_tree,
                gene_to_containing_taxon_map = gene_to_species_map,
                edge_pop_size_attr = "pop_size",
                rng = rng)
        gene_tree.label = str("{}".format(i))
        gene_trees.append(gene_tree)
    return gene_trees

def sim_species_tree(rng, config):
    sp_tree = dp.simulate.treesim.birth_death_tree(
        birth_rate=config["birth_rate"],
        death_rate=0,
        num_extant_tips=config["nspecies"],
        gsa_ntax=config["nspecies"] + 1,
        rng=rng)
    sp_tree.seed_node.edge.length = 0.0
    sp_tree.calc_node_ages()
    for node in sp_tree:
        node.annotations.drop()
    for node in sp_tree:
        pop_size = rng.gammavariate(
            config["pop_size_shape"], 
            config["pop_size_scale"]) 
        node.edge.pop_size = pop_size
        node.annotations['pop_size'] = node.edge.pop_size
    return sp_tree

def sim_gene_trees(sp_tree, rng, config):
    taxon_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping(
        sp_tree.taxon_namespace,
        config["ngenomes"])
    gene_trees = dp.TreeList()
    for i in range(config["nloci"]):
        gene_tree = dp.simulate.treesim.contained_coalescent_tree(
                containing_tree = sp_tree,
                gene_to_containing_taxon_map = taxon_map,
                edge_pop_size_attr = "pop_size",
                rng = rng)
        gene_tree.label = str("{}".format(i))
        gene_trees.append(gene_tree)
    return gene_trees

def gen_alignments(dir, rng, config, gene_trees):
    for ix, tree in enumerate(gene_trees):
        out_path = os.path.join(dir, "alignment-{}.phy".format(ix))
        newick_string = tree.as_string("newick", suppress_rooting=True)
        seq_gen_script = [
            "seq-gen",
            "-z", str(rng.randint(0, 1000000000)),
            "-l", str(config["locus_length"]),
            "-m", "gtr",
            "-f", "0.0, 0.0, 0.5, 0.5",
            "-r", "0, 0, 0, 0, 0, 1.0", 
            "-op", 
            "-q",
            "<<<", "\"{}\"".format(newick_string),
            ">", out_path]
        subprocess.call(" ".join(seq_gen_script), shell=True)

def gen_data(seed, config_path, eco_path):
    config = yaml.safe_load(open(config_path))
    eco_config = yaml.safe_load(open(eco_path))

    rng = random.Random(seed)

    # Create output directories
    out_dir = os.path.join(
            os.path.abspath(
                "out-sp{sp}-gen{gen}-loc{loc}-len{len}".format(
                    sp=config["nspecies"],
                    gen=config["ngenomes"],
                    loc=config["nloci"],
                    len=config["locus_length"])),
            "seed{}-reps{}".format(seed, config["nreps"]),
            "singleton-prob-1.0")
    os.makedirs(out_dir) 

   # Copy configs into output directory
    config["seed"] = rng.randint(0, 1000000000)
    new_config_path = os.path.join(out_dir, "config.yml")
    yaml.dump(config, open(new_config_path, "w"))
    new_eco_config_path = os.path.join(out_dir, "eco-config.yml")
    yaml.dump(eco_config, open(new_eco_config_path, "w"))

    # Simulate data
    for rep in range(0, config["nreps"]):
        rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
        os.mkdir(rep_dir)

        # Species tree
        sp_tree = simulate_species_tree(
            n_sp=config["nspecies"],
            birth_rate=config["birth_rate"],
            death_rate=0,
            pop_size_shape=config["pop_size_shape"],
            pop_size_scale=config["pop_size_scale"],
            rng=rng)
        sp_tree.write(
            path=os.path.join(rep_dir, "species_tree.nex"), 
            schema="nexus")

        # Gene trees
        gene_trees = simulate_gene_trees(
            sp_tree=sp_tree,
            n_tips_per_sp=config["ngenomes"],
            n_loci=config["nloci"],
            rng=rng)
        gene_trees.write(
            path=os.path.join(rep_dir, "gene_trees.nex"), 
            schema="nexus")

        # # Simulate trees
        # sp_tree = sim_species_tree(rng=rng, config=config)
        # sp_tree.write(
        #     path=os.path.join(rep_dir, "sp-tree.tre"), 
        #     schema="newick", 
        #     suppress_rooting=True)

        # # Generate gene trees
        # gene_trees = sim_gene_trees(sp_tree=sp_tree, rng=rng, config=config)
        # gene_trees.write(
        #     path=os.path.join(rep_dir, "gene-trees.tre"), 
        #     schema="newick", 
        #     suppress_rooting=True)

        # Generate alignments    
        gen_alignments(
            dir=rep_dir, 
            rng=rng,
            config=config,
            gene_trees=gene_trees)

if __name__ == "__main__":
    fire.Fire(gen_data)
