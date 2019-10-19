#!/usr/bin/env python

import fire
import glob
import yaml
import fire
import subprocess
import os
import dendropy as dp
import random
import numpy as np

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

def simulate_gene_tree(sp_tree, n_tips_per_sp, rng):
    name_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping(
        sp_tree.taxon_namespace,
        n_tips_per_sp)
    gene_tree = dp.simulate.treesim.contained_coalescent_tree(
        containing_tree=sp_tree,
        gene_to_containing_taxon_map=name_map,
        edge_pop_size_attr="pop_size",
        rng=rng)
    return gene_tree

def simulate_alignment(path, rng, gene_tree, locus_length):
    newick_string = gene_tree.as_string("newick", suppress_rooting=True)
    seq_gen_script = [
        "seq-gen",
        "-z", str(rng.randint(0, 1000000000)),
        "-l", str(locus_length),
        "-m", "gtr",
        "-f", "0.0, 0.0, 0.5, 0.5",
        "-r", "0, 0, 0, 0, 0, 1.0", 
        "-op", 
        "-q",
        "<<<", "\"{}\"".format(newick_string),
        ">", path]
    subprocess.call(" ".join(seq_gen_script), shell=True)

def check_similarity(align_path, similarity_thresh):
    align = dp.DnaCharacterMatrix.get(
        path=align_path,
        schema="phylip")
    seg_sites = dp.calculate.popgenstat.num_segregating_sites(align)
    proportion_similar = (align.max_sequence_size - seg_sites) / align.max_sequence_size
    if proportion_similar  > similarity_thresh:
        is_similar = True
    else:
        is_similar = False
    return is_similar 

def gen_data(seed, config_path, eco_path, similarity_thresh=0.0):
    """
    similarity_thresh: minimum proportion of sites that must be conserved
    """
    max_loci = 1000000 
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
    discarded_species_trees = dp.TreeList()

    # Simulate data
    for rep in range(0, config["nreps"]):
        print("******************************************")
        print("Rep: " + str(rep))
        rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
        os.mkdir(rep_dir)
        state = "speciesTree"
        while True:
            if state == "speciesTree":
                print("----Species Tree")
                sp_tree = simulate_species_tree(
                    n_sp=config["nspecies"],
                    birth_rate=config["birth_rate"],
                    death_rate=0,
                    pop_size_shape=config["pop_size_shape"],
                    pop_size_scale=config["pop_size_scale"],
                    rng=rng)
                state = "geneTree"
                locus = 0
                loci_tried = 0
                gene_trees = dp.TreeList()
            # Simulate a gene tree and an alignment
            elif state == "geneTree":
                print("--------Locus: " + str(locus))
                print("--------Attempt: " + str(loci_tried))
                gene_tree = simulate_gene_tree(
                    sp_tree=sp_tree,
                    n_tips_per_sp=config["ngenomes"],
                    rng=rng)
                align_path = os.path.join(rep_dir, "alignment-{}.phy".format(locus))
                simulate_alignment(
                    path=align_path,
                    rng=rng,
                    gene_tree=gene_tree,
                    locus_length=config["locus_length"])              
                if similarity_thresh > 0.0:
                    state = "checkSimilarity"
                else:
                    state = "saveGene"
            # Check proportion of similarity
            elif state == "checkSimilarity":
                print("------------Check Similarity")
                similar = check_similarity(
                    align_path=align_path,
                    similarity_thresh=similarity_thresh)
                if similar:
                    print("----------------Similar: True")
                    state = "saveGene"
                    loci_tried += 1
                else:
                    print("--------------------Similar: False")
                    if loci_tried < max_loci - 1:
                        loci_tried += 1
                        state = "geneTree"
                    else:
                        state = "speciesTree"
                        discarded_species_trees.append(sp_tree)
                        print("--------------------------------Discard Species Tree")
            # Add gene tree to list
            elif state == "saveGene":
                print("-------------------Save Gene")
                gene_tree.label = str("{}".format(locus))
                gene_trees.append(gene_tree)
                if locus < config["nloci"] - 1:
                    state = "geneTree"
                    locus += 1
                else:
                    state = "finish"
            # Simulations for current replicate complete
            elif state == "finish":
                print("Finish")
                sp_tree.write(
                    path=os.path.join(rep_dir, "species_tree.nex"), 
                    schema="nexus")
                gene_trees.write(
                    path=os.path.join(rep_dir, "gene_trees.nex"), 
                    schema="nexus")
                break
    # Save any discarded species trees
    if len(discarded_species_trees) > 0:
        discarded_species_trees.write(
            path=os.path.join(out_dir, "discarded_species_trees.nex"),
            schema="nexus")


if __name__ == "__main__":
    fire.Fire(gen_data)
