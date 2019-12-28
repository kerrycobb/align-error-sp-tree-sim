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

def simulate_gene_trees(sp_tree, nloci, ngenomes, rng):
    name_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping(
        sp_tree.taxon_namespace,
        ngenomes)
    gene_trees = dp.TreeList()
    for locus in range(0, nloci):
        gene_tree = dp.simulate.treesim.contained_coalescent_tree(
            containing_tree=sp_tree,
            gene_to_containing_taxon_map=name_map,
            edge_pop_size_attr="pop_size",
            rng=rng)
        gene_tree.label = str(locus)
        gene_trees.append(gene_tree)
    return gene_trees

def simulate_alignment(trees, length, rng, dir, concat=False):
    s = seqgen.SeqGen()
    s.char_model = seqgen.SeqGen.GTR
    s.state_freqs = [0, 0, 0.5, 0.5]
    s.general_rates = [0, 0, 0, 0, 0, 1.0]
    s.seq_len = length
    s._rng = rng
    dataset = s.generate(trees)
    if concat:
        path = os.path.join(dir, "alignment.phy")
        align = dp.DnaCharacterMatrix.concatenate(dataset.char_matrices)
        for i in dataset.char_matrices:
            for taxon in align.taxon_namespace:
                taxon.label = taxon.label.replace(' ', '_')
        align.write(path=path, schema="phylip")
    else:
        for i, align in enumerate(dataset.char_matrices):
            for taxon in align.taxon_namespace:
                taxon.label = taxon.label.replace(' ', '_')
            path = os.path.join(dir, "alignment-{}.phy".format(i))
            align.write(path=path, schema="phylip")

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


def gen_data(seed, config_path, eco_path, concat=False):
    """
    concat: Concatenate loci into single alignment
    """
    config = yaml.safe_load(open(config_path))
    eco_config = yaml.safe_load(open(eco_path))
    rng = random.Random(seed)

    # Make output directory
    out_dir = os.path.join(
        "out-sp{sp}-gen{gen}-loc{loc}-len{len}".format(
            sp=config["nspecies"],
            gen=config["ngenomes"],
            loc=config["nloci"],
            len=config["locus_length"]),
        "seed{}-reps{}".format(seed, config["nreps"]),
        "singleton-prob-1.0")
    os.makedirs(out_dir) 

   # Copy configs into output directory
    config["seed"] = rng.randint(0, 1000000000)
    new_config_path = os.path.join(out_dir, "config.yml")
    yaml.dump(config, open(new_config_path, "w"))
    new_eco_config_path = os.path.join(out_dir, "eco-config.yml")
    yaml.dump(eco_config, open(new_eco_config_path, "w"))

    for rep in tqdm(range(0, config["nreps"])):
        rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
        os.mkdir(rep_dir)
    
        # Simulate species tree
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

        # Simulate gene trees
        gene_trees = simulate_gene_trees(
            sp_tree=sp_tree,
            nloci=config["nloci"],
            ngenomes=config["ngenomes"],
            rng=rng)
        gene_trees.write(
            path=os.path.join(rep_dir, "gene_trees.nex"), 
            schema="nexus")

        # Simulate alignments
        simulate_alignment(
            trees=gene_trees, 
            length=config["locus_length"], 
            rng=rng, 
            dir=rep_dir, 
            concat=concat)
    
    print("Simulation Complete")


    # # Simulate data
    # for rep in tqdm(range(0, config["nreps"])):
    #     rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
    #     os.mkdir(rep_dir)
    #     state = "speciesTree"
    #     while True:
    #         if state == "speciesTree":
    #             sp_tree = simulate_species_tree(
    #                 n_sp=config["nspecies"],
    #                 birth_rate=config["birth_rate"],
    #                 death_rate=0,
    #                 pop_size_shape=config["pop_size_shape"],
    #                 pop_size_scale=config["pop_size_scale"],
    #                 rng=rng)
    #             state = "geneTree"
    #             locus = 0
    #             loci_tried = 0
    #             name_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping(
    #                 sp_tree.taxon_namespace,
    #                 config["ngenomes"])
    #             gene_trees = dp.TreeList()
    #         # Simulate a gene tree and an alignment
    #         elif state == "geneTree":
    #             # gene_tree = simulate_gene_tree(
    #             #     sp_tree=sp_tree,
    #             #     n_tips_per_sp=config["ngenomes"],
    #             #     rng=rng)
    #             gene_tree = dp.simulate.treesim.contained_coalescent_tree(
    #                 containing_tree=sp_tree,
    #                 gene_to_containing_taxon_map=name_map,
    #                 edge_pop_size_attr="pop_size",
    #                 rng=rng)
    #             align_path = os.path.join(rep_dir, "alignment-{}.phy".format(locus))
    #             simulate_alignment(
    #                 path=align_path,
    #                 rng=rng,
    #                 gene_tree=gene_tree,
    #                 locus_length=config["locus_length"])              
    #             if similarity_thresh > 0.0:
    #                 state = "checkSimilarity"
    #             else:
    #                 state = "saveGene"
    #         # Check proportion of similarity
    #         elif state == "checkSimilarity":
    #             similar = check_similarity(
    #                 align_path=align_path,
    #                 similarity_thresh=similarity_thresh)
    #             if similar:
    #                 state = "saveGene"
    #                 loci_tried += 1
    #             else:
    #                 if loci_tried < max_loci - 1:
    #                     loci_tried += 1
    #                     state = "geneTree"
    #                 else:
    #                     state = "speciesTree"
    #                     discarded_species_trees.append(sp_tree)
    #                     print("Discarding species tree")
    #         # Add gene tree to list
    #         elif state == "saveGene":
    #             gene_tree.label = str("{}".format(locus))
    #             gene_trees.append(gene_tree)
    #             if locus < config["nloci"] - 1:
    #                 state = "geneTree"
    #                 locus += 1
    #             else:
    #                 state = "finish"
    #         # Simulations for current replicate complete
    #         elif state == "finish":
    #             sp_tree.write(
    #                 path=os.path.join(rep_dir, "species_tree.nex"), 
    #                 schema="nexus")
    #             gene_trees.write(
    #                 path=os.path.join(rep_dir, "gene_trees.nex"), 
    #                 schema="nexus")
    #             break
    # # Save any discarded species trees
    # if len(discarded_species_trees) > 0:
    #     discarded_species_trees.write(
    #         path=os.path.join(out_dir, "discarded_species_trees.nex"),
    #         schema="nexus")


if __name__ == "__main__":
    fire.Fire(gen_data)
