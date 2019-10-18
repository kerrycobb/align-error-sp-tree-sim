#!/usr/bin/env python

import fire
import glob
import yaml
import fire
import subprocess
import os
# import jinja2 as jj
import dendropy as dp
# import math
import random
import numpy as np
# from Bio import AlignIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment
# from Bio.Alphabet import generic_dna
# import pyvolve as pv

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

# def simulate_gene_trees(sp_tree, n_tips_per_sp, n_loci, rng):
#     gene_to_species_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping(
#         sp_tree.taxon_namespace,
#         n_tips_per_sp)
#     gene_trees = dp.TreeList()
#     for i in range(n_loci):
#         gene_tree = dp.simulate.treesim.contained_coalescent_tree(
#                 containing_tree=sp_tree,
#                 gene_to_containing_taxon_map=gene_to_species_map,
#                 edge_pop_size_attr="pop_size",
#                 rng=rng)
#         gene_tree.label = str("{}".format(i))
#         gene_trees.append(gene_tree)
#     return gene_trees

# def simulate_alignments(dir, rng, gene_trees, locus_length):
#     for ix, tree in enumerate(gene_trees):
#         out_path = os.path.join(dir, "alignment-{}.phy".format(ix))
#         newick_string = tree.as_string("newick", suppress_rooting=True)
#         seq_gen_script = [
#             "seq-gen",
#             "-z", str(rng.randint(0, 1000000000)),
#             "-l", str(locus_length),
#             "-m", "gtr",
#             "-f", "0.0, 0.0, 0.5, 0.5",
#             "-r", "0, 0, 0, 0, 0, 1.0", 
#             "-op", 
#             "-q",
#             "<<<", "\"{}\"".format(newick_string),
#             ">", out_path]
#         subprocess.call(" ".join(seq_gen_script), shell=True)


# def gen_data(seed, config_path, eco_path, similarity_thresh=0.0):
#     """
#     similarity_thresh: minimum proportion of sites that must be conserved
#     """
#     config = yaml.safe_load(open(config_path))
#     eco_config = yaml.safe_load(open(eco_path))
#     rng = random.Random(seed)
    
#     # Create output directories
#     out_dir = os.path.join(
#             os.path.abspath(
#                 "out-sp{sp}-gen{gen}-loc{loc}-len{len}".format(
#                     sp=config["nspecies"],
#                     gen=config["ngenomes"],
#                     loc=config["nloci"],
#                     len=config["locus_length"])),
#             "seed{}-reps{}".format(seed, config["nreps"]),
#             "singleton-prob-1.0")
#     os.makedirs(out_dir) 

#    # Copy configs into output directory
#     config["seed"] = rng.randint(0, 1000000000)
#     new_config_path = os.path.join(out_dir, "config.yml")
#     yaml.dump(config, open(new_config_path, "w"))
#     new_eco_config_path = os.path.join(out_dir, "eco-config.yml")
#     yaml.dump(eco_config, open(new_eco_config_path, "w"))

#     # Simulate data
#     for rep in range(0, config["nreps"]):
#         rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
#         os.mkdir(rep_dir)           
      
#         sp_tree = simulate_species_tree(
#             n_sp=config["nspecies"],
#             birth_rate=config["birth_rate"],
#             death_rate=0,
#             pop_size_shape=config["pop_size_shape"],
#             pop_size_scale=config["pop_size_scale"],
#             rng=rng)
 
#         gene_trees = dp.TreeList()
#         for locus in range(0, config["nloci"]):
#             gene_tree = simulate_gene_tree(
#                 sp_tree=sp_tree,
#                 n_tips_per_sp=config["ngenomes"],
#                 rng=rng)
#             gene_tree.label = str("{}".format(locus))
#             gene_trees.append(gene_tree)
#             align_path = os.path.join(rep_dir, "alignment-{}.phy".format(locus))
#             simulate_alignment(
#                 path=align_path,
#                 rng=rng,
#                 gene_tree=gene_tree,
#                 locus_length=config["locus_length"])              
       
#         sp_tree.write(
#                 path=os.path.join(rep_dir, "species_tree.nex"), 
#                 schema="nexus")
#         gene_trees.write(
#                 path=os.path.join(rep_dir, "gene_trees.nex"), 
#                 schema="nexus")

def gen_data(seed, config_path, eco_path, similarity_thresh=0.0):
    """
    similarity_thresh: minimum proportion of sites that must be conserved
    """
    max_loci = 20
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
        print("Rep:" + str(rep))
        rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
        os.mkdir(rep_dir)           

        cnt = 0
        while True:
            print("Species tree attempt:" + str(cnt))
            sp_tree = simulate_species_tree(
                n_sp=config["nspecies"],
                birth_rate=config["birth_rate"],
                death_rate=0,
                pop_size_shape=config["pop_size_shape"],
                pop_size_scale=config["pop_size_scale"],
                rng=rng)


            while True:
                gene_tree = simulate_gene_tree(
                    sp_tree=sp_tree,
                    n_tips_per_sp=config["ngenomes"],
                    rng=rng)

            if cnt < 5:
                cnt += 1
            else:
                break

#TODO Finish this part


        # gene_trees = dp.TreeList()
        # for locus in range(0, config["nloci"]):
        #     gene_tree = simulate_gene_tree(
        #         sp_tree=sp_tree,
        #         n_tips_per_sp=config["ngenomes"],
        #         rng=rng)
        #     gene_tree.label = str("{}".format(locus))
        #     gene_trees.append(gene_tree)
        #     align_path = os.path.join(rep_dir, "alignment-{}.phy".format(locus))
        #     simulate_alignment(
        #         path=align_path,
        #         rng=rng,
        #         gene_tree=gene_tree,
        #         locus_length=config["locus_length"])              
       
        # sp_tree.write(
        #         path=os.path.join(rep_dir, "species_tree.nex"), 
        #         schema="nexus")
        # gene_trees.write(
        #         path=os.path.join(rep_dir, "gene_trees.nex"), 
        #         schema="nexus")


        # print("******************Rep: {}".format(rep))
        # rep_dir = os.path.join(out_dir, "rep-{}".format(rep))
        # os.mkdir(rep_dir)
        # sp_trees_tried = 0
        # while True:
        #     loci_tried = 0
        #     loci = 0
        #     while True:
        #         n


            # print("Species tree: {}".format(rep)) 
            # # Species tree
            # sp_tree = simulate_species_tree(
            #     n_sp=config["nspecies"],
            #     birth_rate=config["birth_rate"],
            #     death_rate=0,
            #     pop_size_shape=config["pop_size_shape"],
            #     pop_size_scale=config["pop_size_scale"],
            #     rng=rng)
                   
            # # Simulate gene trees and alignments
            # locus = 0
            # loci_tried = 0
            # gene_trees = dp.TreeList()
            # while True:
            #     if locus == nloci:
            #         break
            #     print("Locus: {}".format(locus))
            #     align_path = os.path.join(rep_dir, "alignment-{}.phy".format(locus))
            #     while True:
            #         print(loci_tried)
            #         gene_tree = simulate_gene_tree(
            #             sp_tree=sp_tree,
            #             n_tips_per_sp=n_tips_per_sp,
            #             rng=rng)
            #         simulate_alignment(
            #             path=align_path,
            #             rng=rng,
            #             gene_tree=gene_tree,
            #             locus_length=locus_length)
            #         if similarity_thresh > 0.0:
            #             is_similar = check_similarity(
            #                 align_path=align_path,
            #                 similarity_thresh=similarity_thresh) 
            #             if is_similar:
            #                 sp_tree_state = False
            #                 break
            #             else:
            #                 if loci_tried < max_loci: 
            #                     loci_tried += 1
            #                 else:
            #                     print("1,000,000 loci tried, simulating new species tree")
            #                     break
            #         else:
            #             sp_tree_state = False
            #             break
            #     if loci_tried >= max_loci:
            #         print("Breaking locus loop. Locus = {}".format(locus))
            #         break
            #     else:
            #         gene_tree.label = str("{}".format(locus))
            #         gene_trees.append(gene_tree)
            #         locus += 1
            #         loci_tried += 1  
            
            # # Output tree files 
            # sp_tree.write(
            #     path=os.path.join(rep_dir, "species_tree.nex"), 
            #     schema="nexus")
            # gene_trees.write(
            #     path=os.path.join(rep_dir, "gene_trees.nex"), 
            #     schema="nexus")
            # break
            
        # # Gene trees
        # gene_trees = simulate_gene_trees(
        #     sp_tree=sp_tree,
        #     n_tips_per_sp=config["ngenomes"],
        #     n_loci=config["nloci"],
        #     rng=rng)
        # gene_trees.write(
        #     path=os.path.join(rep_dir, "gene_trees.nex"), 
        #     schema="nexus")

        # # Generate alignments    
        # simulate_alignments(
        #     dir=rep_dir, 
        #     rng=rng,
        #     gene_trees=gene_trees,
        #     locus_length=config["locus_length"])

if __name__ == "__main__":
    fire.Fire(gen_data)
