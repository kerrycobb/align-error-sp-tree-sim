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

def clean(dir):
    for i in glob.glob(os.path.join(dir, "rep-*")):
        try:
            os.remove(os.path.join(i, "alignment.nex"))
        except:
            pass
        try:
            os.remove(os.path.join(i, "starbeast.xml"))
        except:
            pass
        for j in glob.glob(os.path.join(i, "ecoevo-chain-*")):
            try:
                rmtree(j)
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
        sample_freq=int(config["starbeast_chain_length"]/10000),
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
    
    # # seq_records = [] 
    # # for key, value in d.items():
    # #     d[key] = "".join(value)
    # #     seq_records.append(SeqRecord(Seq("".join(value), SingleLetterAlphabet("01")), id=key))
    # # AlignIO.write(
    # #     MultipleSeqAlignment(seq_records),
    # #     os.path.join(dir, "alignment.nex"),
    # #     "nexus")

    # new_alignment = dp.StandardCharacterMatrix.get(path=os.path.join(dir, "alignment.nex"), schema="nexus")
    # new_alignment = dp.StandardCharacterMatrix.from_dict(d)
    # # new_alignment.write(path=os.path.join(dir, "alignment.nex"), schema="nexus")
    # print(new_alignment.as_string(schema="nexus"))

def qsub(dir, jobname, walltime, script, ssh=False):
    qsub = [
        "/cm/shared/apps/torque/6.1.1.1.h3/bin/qsub",
        "-m", "n",
        "-N", jobname,
        "-l", "nodes=1:ppn=1",
        "-l", "mem=1G",
        "-l", "walltime={}".format(walltime),
        "-j", "oe",
        "-d", dir,
        "-o", dir,
        ]
    if time.strptime(walltime, "%H:%M:%S") > time.strptime("1:0:0", "%H:%M:%S"):
        rand = random.random()
        if rand < .5:
            qsub.extend([
            "-q", "general",
            "-W", "group_list=jro0014_lab",
            "-W", "x=FLAGS:ADVRES:jro0014_lab",
            ])
        else:
            qsub.extend([
                "-q", "gen28",
                "-W", "group_list=jro0014_lab",
                "-W", "x=FLAGS:ADVRES:jro0014_s28",
            ])
    ssh_script = ""
    if ssh:
        ssh_script = "ssh kac0070@hopper.auburn.edu "
    cmd = "{ssh}{qsub} <<< \"{script}\"".format(
        ssh=ssh_script, qsub=" ".join(qsub), script=script)
    proc = subprocess.call(cmd, shell=True)

def run_starbeast(dir, rep, config, rng, ssh=False):
    for chain in range(1, config["starbeast_chains"]+1):
        xml = os.path.join(dir, "starbeast.xml")
        chain_dir = os.path.join(dir, "starbeast-chain-{}".format(chain))
        os.mkdir(chain_dir)
        seed = str(rng.randint(1,1000000))
        with open(os.path.join(chain_dir, "seed.txt"), "w") as fh:
            fh.write(seed)
        script = [
            "module load beast \n",
            "module load beagle \n",
            "cd {}\n".format(chain_dir),
            "beast",
            "-beagle",
            "-seed", seed,
            xml]
        qsub(
            dir=chain_dir,
            jobname="star-{}-{}".format(rep, chain),
            walltime="4:00:00",
            script=" ".join(script),
            ssh=ssh)
 
def run_ecoevolity(dir, rep, rng, config, eco_config, ssh=False):
    eco_config["comparisons"][0]["comparison"]["path"] = os.path.join(
        dir, "alignment.nex")
    for chain in range(1, config["ecoevolity_chains"]+1):
        chain_dir = os.path.join(dir, "ecoevo-chain-{}".format(chain))
        if not os.path.exists(chain_dir):
            os.mkdir(chain_dir)
        config_path = os.path.join(chain_dir, "ecoevolity-config.yml")
        yaml.dump(eco_config, open(config_path, "w"))
        seed = str(rng.randint(1, 1000000))
        with open(os.path.join(chain_dir, "seed.txt"), "w") as fh:
            fh.write(seed)
        script = [
            "cd {}\n".format(chain_dir),
            "ecoevolity-ar",
            "--seed", seed,
            "--relax-constant-sites",
            "ecoevolity-config.yml"]
        qsub(
            dir=chain_dir,
            jobname="ecoevo-{}-{}".format(rep, chain),
            walltime="1:00:00",
            script=" ".join(script),
            ssh=ssh)

def run_analyses(dir, ssh=False, overwrite=False):
    dir = os.path.abspath(dir)
    if overwrite:
        clean(dir)
    config_path = os.path.join(dir, "config.yml")
    config = yaml.safe_load(open(config_path))
    eco_config_path = os.path.join(dir, "eco-config.yml")
    eco_config = yaml.safe_load(open(eco_config_path))
    rng = random.Random(config["seed"])
    
    for rep in range(0, config["nreps"]):
        rep_dir = os.path.join(dir, "rep-{}".format(rep))
        gen_xml(rep_dir, config)
        gen_ecoevol_alignment(rep_dir, config)
        run_starbeast(rep_dir, rep, config, rng, ssh=ssh)
        run_ecoevolity(rep_dir, rep, rng, config, eco_config, ssh=ssh)

if __name__ == "__main__":
    fire.Fire(run_analyses)
