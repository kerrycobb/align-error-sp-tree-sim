#!/usr/bin/env python

import fire
import yaml
import fire
import subprocess
import os
import jinja2 as jj
import dendropy as dp
import math
import random

def gen_out_dirs(seed, config):
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
    os.mkdir(out_dir)
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
    for ix, tree in enumerate(gene_trees):
        alignment_path = os.path.join(dir, "alignment-{}.phy".format(ix))
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
            ">", alignment_path]
        subprocess.call(" ".join(seq_gen_script), shell=True)

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
    
def gen_data(seed, config_path):
    subprocess.call("module load seq-gen; module load simphy", shell=True)
    config = yaml.safe_load(open(config_path))

    # Create output directories
    out_dir = gen_out_dirs(seed, config)

    # Copy configs into output directory
    new_config_path = os.path.join(out_dir, "config.yml")
    yaml.dump(config, open(new_config_path, "w"))
    # new_eco_config_path = os.path.join(out_dir, "eco-config.yml")
    # yaml.dump(eco_config, open(new_eco_config_path, "w"))

    rng = random.Random(seed)

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
        # Generate gene trees    
        gen_alignments(
            dir=rep_dir, 
            rng=rng,
            config=config,
            gene_trees=gene_trees)

        # Generate xml 
        gen_xml(dir=rep_dir, config=config)
        # Generate alignment for ecoevolity
        gen_ecoevol_alignment(dir=rep_dir, config=config, gene_trees=gene_trees)

if __name__ == "__main__":
    fire.Fire(gen_data)
