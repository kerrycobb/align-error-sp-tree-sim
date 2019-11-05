#!/usr/bin/env python

import fire
import glob
import random
import os
import yaml
import time
import subprocess
import random
import math
from shutil import copyfile, rmtree
import jinja2 as jj
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna
import dendropy as dp

def clean(dir, method):
    for i in glob.glob(os.path.join(dir, "rep-*")):
        if method in ["all", "ecoevolity"]: 
            try:
                os.remove(os.path.join(i, "alignment.nex"))
            except:
                pass
            for j in glob.glob(os.path.join(i, "ecoevo-chain-*")):
                try:
                    rmtree(j)
                except:
                    pass
        if method in ["all", "starbeast"]:
            try:
                os.remove(os.path.join(i, "starbeast.xml"))
            except:
                pass
            for j in glob.glob(os.path.join(i, "starbeast-chain-*")):
                try:
                    rmtree(j)
                except:
                    pass

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
        alignment = AlignIO.read(open(align_path), "phylip") 
        seq_dicts = []
        for record in alignment:
            taxon, tip = record.id.split("_")
            seq_dicts.append(dict(taxon=taxon, tip=tip, seq=record.seq))
        alignment_dicts.append(dict(id=i, sequences=seq_dicts))
    starbeast_xml = template.render(
        chain_length=config["starbeast_chain_length"],
        sample_freq=config["starbeast_sample_freq"],
        alignments=alignment_dicts,
        taxa=taxa,
        tips=tips,
        locus_length=config["locus_length"],
        n_loci=config["nloci"],
        pop_size_shape=config["pop_size_shape"],
        pop_size_scale=1/config["pop_size_scale"],
        birth_rate=config["birth_rate"])
    xml_path = os.path.join(dir, "starbeast.xml")
    with open(xml_path, 'w') as handle:
        handle.write(starbeast_xml)
    return xml_path

def gen_ecoevol_alignment(dir, config):
    d = {}
    for i in range(1, config["nspecies"]+1):
        for j in range(1, config["ngenomes"]+1):
            tax = "T{}_{}".format(i, j) 
            d[tax] = [] 
    for i in range(0, config["nloci"]):
        align_path = os.path.join(dir, "alignment-{}.phy".format(i))
        alignment = AlignIO.read(open(align_path), "phylip") 
        for record in alignment:
            for j in record.seq:
                if j == "T":
                    d[record.id].append("0")
                elif j == "G":
                    d[record.id].append("1")
                else:
                    raise Exception("Invalid character {}".format(j))
    ns = []
    ns.append("#NEXUS\n")
    ns.append("BEGIN TAXA;\n")
    ns.append("DIMENSIONS NTAX={};\n".format(len(d)))
    ns.append("TAXLABELS\n")
    for key in d.keys():
        ns.append("'{}'\n".format(key))
    ns.append(";\nEND;\n")
    ns.append("BEGIN CHARACTERS;\n")
    ns.append("DIMENSIONS NCHAR={};\n".format(config["locus_length"] * config["nloci"]))
    ns.append("FORMAT DATATYPE=STANDARD SYMBOLS=\"01\" MISSING=?;\n")
    ns.append("MATRIX\n")
    for key, value in d.items():
        ns.append("'{}'    {}\n".format(key, "".join(value)))
    ns.append(";\nEND;")

    with open(os.path.join(dir, "alignment.nex"), "w") as fh:
        fh.writelines(ns)
    
def qsub(dir, jobname, walltime, script, ssh=False):
    cmd = ["myqsub", "-N", jobname, "-d", dir, "-t", walltime]
    if ssh:
        cmd.append("-s")
    cmd.append("\"{}\"".format(script)) 
    subprocess.call(" ".join(cmd), shell=True)

def run_starbeast(dir, seed, jobname, xml, rerun=False, ssh=False):
    script = [
        "module load beast \n",
        "module load beagle \n",
        "cd {}\n".format(dir),
        "beast",
        "-beagle",
        "-seed", seed]
    if rerun:
        script.append("-overwrite")
    script.append(xml)
    qsub(dir=dir, jobname=jobname, walltime="200:00:00",
            script=" ".join(script), ssh=ssh)
   
def run_ecoevolity(dir, seed, jobname, eco_config, rerun=False, ssh=False):
    if rerun:
        try:
            os.remove(os.path.join(dir, "ecoevolity-config-operator-run-1.log"))
            os.remove(os.path.join(dir, "ecoevolity-config-state-run-1.log"))
        except:
            pass       
    script = [
        "cd {}\n".format(dir),
        "ecoevolity-ar",
        "--seed", seed,
        "--relax-constant-sites",
        "ecoevolity-config.yml"]
    qsub(dir=dir, jobname=jobname,
            walltime="1:00:00", script=" ".join(script), ssh=ssh)

def run_analyses(dir, ssh=False, overwrite=False, method="all"):
    dir = os.path.abspath(dir)
    basename = os.path.basename(dir)
    if overwrite:
        clean(dir, method)
    config_path = os.path.join(dir, "config.yml")
    config = yaml.safe_load(open(config_path))
    eco_config_path = os.path.join(dir, "eco-config.yml")
    eco_config = yaml.safe_load(open(eco_config_path))
    rng = random.Random(config["seed"])

    for rep in range(0, config["nreps"]):
        rep_dir = os.path.join(dir, "rep-{}".format(rep))

        # Create ecoevolity directories and files and run ecoevolity
        if method in ["all", "ecoevolity"]: 
            gen_ecoevol_alignment(rep_dir, config)
            eco_config["comparisons"][0]["comparison"]["path"] = os.path.join(
                    rep_dir, "alignment.nex")
            for chain in range(1, config["ecoevolity_chains"]+1):
                chain_dir = os.path.join(rep_dir, "ecoevo-chain-{}".format(chain))
                os.mkdir(chain_dir)
                eco_config_path = os.path.join(
                        chain_dir, "ecoevolity-config.yml")
                yaml.dump(eco_config, open(eco_config_path, "w"))
                seed = str(rng.randint(1, 1000000))
                with open(os.path.join(chain_dir, "seed.txt"), "w") as fh:
                    fh.write(seed)
                jobname = "{}-{}-eco-{}-{}".format(
                        basename, config["locus_length"], rep, chain) 
                run_ecoevolity(dir=chain_dir, seed=seed, jobname=jobname, 
                        eco_config=eco_config, ssh=ssh)
        # Creat starbeast directories and run starbeast 
        if method in ["all", "starbeast"]:
            for chain in range(1, config["starbeast_chains"] + 1):
                chain_dir = os.path.join(rep_dir, 
                        "starbeast-chain-{}".format(chain))
                os.mkdir(chain_dir)
                seed = str(rng.randint(1,1000000))
                with open(os.path.join(chain_dir, "seed.txt"), "w") as fh:
                    fh.write(seed)
                xml_path = gen_xml(rep_dir, config)
                jobname = "{}-{}-star-{}-{}".format(
                        basename, config["locus_length"], rep, chain)
                run_starbeast(dir=chain_dir, seed=seed, jobname=jobname, 
                        xml=xml_path, ssh=ssh)

if __name__ == "__main__":
    fire.Fire(run_analyses)
