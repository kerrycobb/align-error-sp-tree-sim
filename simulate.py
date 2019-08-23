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

def gen_out_dirs(seed, config, ignore_existing=False):
    contain_dir = os.path.abspath(
        "out-sp{sp}-gen{gen}-loc{loc}-len{len}".format(
            sp=config["nspecies"],
            gen=config["ngenomes"],
            loc=config["nloci"],
            len=config["locus_length"]))
    if not os.path.exists(contain_dir):
        os.mkdir(contain_dir)
    out_dir = os.path.join(
        contain_dir, 
        "seed{}-reps{}".format(seed, config["nreps"]))
    try:
        os.mkdir(out_dir)
    except:
        if ignore_existing:
            pass
    return out_dir

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
    out_dir = os.path.join(dir, "singleton-prob-1.0")
    os.mkdir(out_dir)
    for ix, tree in enumerate(gene_trees):
        out_path = os.path.join(out_dir, "alignment-{}.phy".format(ix))
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

def gen_error(dir, rng, config):
    for prob in config["singleton_sample_prob"]:
        out_dir = os.path.join(dir, "singleton-prob-{}".format(prob))
        os.mkdir(out_dir)
        singletons = 0
        rmvd_singletons = 0       
        for locus in range(0, config["nloci"]):
            align_path = os.path.join(
                dir, 
                "singleton-prob-1.0", 
                "alignment-{}.phy".format(locus))
            dna = AlignIO.read(open(align_path), "phylip")
            ids = []
            seqs = [] 
            for record in dna:
               ids.append(record.id)
               seqs.append(np.array(list(record.seq))) 
            matrix = np.stack(seqs)
            print(matrix)
            for column in range(0, matrix.shape[1]):
                g = np.count_nonzero(matrix[...,column] == "G")
                t = np.count_nonzero(matrix[...,column] == "T")
                if g == 1:
                    singletons += 1
                    if rng.uniform(0, 1) >= prob: 
                        rmvd_singletons += 1
                        for i in range(0, matrix.shape[0]):
                            matrix[i, column] = "T"                      
                elif t == 1:
                    singletons += 1
                    if rng.uniform(0, 1) >= prob:
                        rmvd_singletons += 1
                        for i in range(0, matrix.shape[0]):
                            matrix[i, column] = "G"
            print(matrix)
            new_seqs = ["".join(row) for row in matrix]
            new_seq_records = [] 
            for id, new_seq in zip(ids, new_seqs):
                new_seq_records.append(SeqRecord(Seq(new_seq), id=id))
            new_align = MultipleSeqAlignment(new_seq_records)
            AlignIO.write(
                new_align, 
                os.path.join(out_dir, "alignment-{}.phy".format(locus)), 
                "phylip") 
        # TODO: Store data that could be used to plot histogram of singletons per locus
        stats = dict(singletons=singletons, removed_singletons=rmvd_singletons)
        with open(os.path.join(out_dir, "stats.yaml"), "w") as fh:
            yaml.dump(stats, fh)

def gen_xml(dir, config):
    taxa = ["T{}".format(i) for i in range(1, config["nspecies"]+1)] 
    tips = [str(i) for i in range(1, config["ngenomes"] + 1)]
    env = jj.Environment(loader=jj.FileSystemLoader(config["templates_path"]),
        autoescape=jj.select_autoescape(['html', 'xml']),
        lstrip_blocks=True,
        trim_blocks=True)
    template = env.get_template("starbeast.xml")     
    alignment_dicts = []
    for i in range(0, config["nloci"]):
        align_path = os.path.join(dir, "alignment-{}.phy".format(i))
        matrix = dp.DnaCharacterMatrix.get(path=align_path, schema="phylip")
        seq_dicts = []
        for tax, seq in matrix.items():
            taxon, tip = tax.label.split("_")
            seq_dicts.append(dict(taxon=taxon, tip=tip, seq=str(seq)))
        alignment_dicts.append(dict(id=i, sequences=seq_dicts))
    starbeast_xml = template.render(
        chain_length=config["starbeast_chains"],
        sample_freq=int(config["starbeast_chains"]/10000),
        alignments=alignment_dicts,
        taxa=taxa,
        tips=tips,
        locus_length=config["locus_length"],
        n_loci=config["nloci"],
        pop_size_shape=config["pop_size_shape"],
        pop_size_scale=1/config["pop_size_scale"],
        birth_rate=config["birth_rate"])
    with open(os.path.join(dir, "starbeast.xml"), 'w') as handle:
        handle.write(starbeast_xml)

def gen_ecoevol_alignment(dir, config, gene_trees):
    d = {}
    for i in range(1, config["nspecies"]+1):
        for j in range(1, config["ngenomes"]+1):
            tax = "T{}_{}".format(i, j) 
            d[tax] = [] 
    for i in range(0, config["nreps"]):
        align_path = os.path.join(dir, "alignment-{}.phy".format(i))
        dna = dp.DnaCharacterMatrix.get(path=align_path, schema="phylip")
        for tax, seq in dna.items():
            for j in str(seq):
                if j == "T":
                    d[tax.label].append("0") 
                elif j == "G":
                    d[tax.label].append("1")
                else:
                    quit("Invalid character {} encountered in {}{}".format(j, dir, align_path ))
    for key, value in d.items():
        d[key] = "".join(value)
    matrix = dp.StandardCharacterMatrix.from_dict(d)
    matrix.write(path=os.path.join(dir, "alignment.nex"), schema="nexus")

def gen_data(seed, config_path, eco_path):
    config = yaml.safe_load(open(config_path))
    eco_config = yaml.safe_load(open(eco_path))

    rng = random.Random(seed)

    # Create output directories
    out_dir = gen_out_dirs(seed, config, ignore_existing=True)

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

        # Simulate trees
        sp_tree = sim_species_tree(rng=rng, config=config)
        sp_tree.write(
            path=os.path.join(rep_dir, "sp-tree.tre"), 
            schema="newick", 
            suppress_rooting=True)

        # Generate gene trees
        gene_trees = sim_gene_trees(sp_tree=sp_tree, rng=rng, config=config)
        gene_trees.write(
            path=os.path.join(rep_dir, "gene-trees.tre"), 
            schema="newick", 
            suppress_rooting=True)

        # Generate alignments    
        gen_alignments(
            dir=rep_dir, 
            rng=rng,
            config=config,
            gene_trees=gene_trees)

        # Generate alignments with error
        gen_error(dir=rep_dir, rng=rng, config=config)

        for align_dir in glob.glob(os.path.join(rep_dir, "singleton-prob-*")): 
            # Generate xml 
            gen_xml(dir=align_dir, config=config)

            # Generate alignment for ecoevolity
            gen_ecoevol_alignment(dir=align_dir, config=config, gene_trees=gene_trees)



if __name__ == "__main__":
    fire.Fire(gen_data)
