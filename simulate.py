#!/usr/bin/env python

import yaml
import fire
import subprocess
import os
import glob
import jinja2 as jj
import dendropy as dp
import math

def generate_xml(
        alignment_path, out_dir, locus_length, total_chars, pop_size_shape,
        pop_size_scale, birth_rate, chain_length, template, template_dir):

    align_base_name = os.path.basename(alignment_path).split(".")[0]
    out_path = os.path.join(out_dir, "starbeast.xml")
    alignment = dp.StandardCharacterMatrix.get(path=alignment_path, schema="nexus")

    # Get sequences and sequence IDs from alignment
    ids = []
    sequences = []
    for tax, seq in alignment.items():
        ids.append(tax.label)
        orig_sequence = seq.symbols_as_string()
        new_sequence = []
        for i in orig_sequence:
            if i == "0":
                new_sequence.append("T")
            elif i == "1":
                new_sequence.append("G")
            else:
                print("Invalid character: {}".format(i))
                quit()
        sequences.append("".join(new_sequence))

    # Variables for XML
    nloci = int(total_chars/locus_length)
    taxa = list({taxon.label.split("-")[0] for taxon in alignment})
    tips = list({taxon.label.split("-")[1] for taxon in alignment})

    # Split sequence into separate loci
    pos = 0
    alignment_dicts = []
    for locus_id in range(1, nloci+1):
        # recoded_alignment_dict = {}
        sequence_dicts = []
        for i, id in enumerate(ids):
            id_split = id.split("-")
            sequence_dicts.append(dict(
                taxon=id_split[0],
                tip=id_split[1],
                seq=sequences[i][pos:pos+locus_length]))
        alignment_dicts.append(dict(
            id=locus_id,
            sequences=sequence_dicts))
        pos += locus_length

    # Build XML file
    env = jj.Environment(loader=jj.FileSystemLoader(template_dir),
        autoescape=jj.select_autoescape(['html', 'xml']),
        lstrip_blocks=True,
        trim_blocks=True)
    template = env.get_template(template)
    xml = template.render(
        chain_length=chain_length,
        sample_freq=int(chain_length / 10000),
        alignments=alignment_dicts,
        taxa=taxa,
        tips=tips,
        locus_length=locus_length,
        n_loci=nloci,
        pop_size_shape=pop_size_shape,
        pop_size_scale=1 / pop_size_scale,
        birth_rate=float(birth_rate))
    with open(out_path, 'w') as handle:
      handle.write(xml)


def gen_data(config_path, eco_config_path, single_snp=False, rmv_cons=False):
    config = yaml.safe_load(open(config_path))
    eco_config = yaml.safe_load(open(eco_config_path))
    seed = config["seed"]
    nreps = config["nreps"]
    singleton_sample_prob = float(config["singleton_sample_prob"])
    nspecies = config["nspecies"]
    ngenomes = config["ngenomes"]
    nloci = config["nloci"]
    locus_length = config["locus_length"]
    ncharacters = nloci * locus_length
    pop_size_shape = eco_config["global_comparison_settings"]["parameters"]\
        ["population_size"]["prior"]["gamma_distribution"]["shape"]
    pop_size_scale = eco_config["global_comparison_settings"]["parameters"]\
        ["population_size"]["prior"]["gamma_distribution"]["scale"]
    wait_time = eco_config["event_time_prior"]["exponential_distribution"]["rate"]
    
    # Create output directories
    snp = ""
    rmvCons = ""
    if single_snp:
        snp = "-snp"
    elif rmv_cons:
        rmvCons = "-rmvcons"
    contain_dir = os.path.abspath(
        "out{snp}{rmvCons}-sp{sp}-gen{gen}-loc{loc}-len{len}-sin{sin}".format(
            snp=snp, rmvCons=rmvCons, sp=nspecies, gen=ngenomes, loc=nloci,
            len=locus_length, sin=singleton_sample_prob))
    if not os.path.exists(contain_dir):
        os.mkdir(contain_dir)
    out_dir = os.path.join(contain_dir, "seed{}-reps{}".format(seed, nreps))
    os.mkdir(out_dir)

    # Copy configs into output directory
    new_config_path = os.path.join(out_dir, "config.yml")
    yaml.dump(config, open(new_config_path, "w"))
    if rmv_cons:
        eco_config["global_comparison_settings"]["constant_sites_removed"] = True
    new_eco_config_path = os.path.join(out_dir, "eco-config.yml")
    yaml.dump(eco_config, open(new_eco_config_path, "w"))

    # Generate dummy alignment
    dummy_alignment_file_name = "dummy-alignment.nex"
    dummy_alignment_path = os.path.join(out_dir, dummy_alignment_file_name)
    dummy_alignment_script = " ".join([
        "./generate-dummy-biallelic-alignment.py",
        "--nspecies", str(nspecies),
        "--ngenomes", str(ngenomes),
        "--ncharacters", str(ncharacters),
        ">", dummy_alignment_path])
    subprocess.call(dummy_alignment_script, shell=True)

    # Simulate data
    simcoevol_script = [
        "simcoevolity-ar",
        "--seed", str(seed),
        "--singleton-sample-prob", str(singleton_sample_prob),
        "--number-of-replicates", str(nreps),
        "-o", out_dir]
    if single_snp:
        simcoevol_script.append("--max-one-variable-site-per-locus")
    if rmv_cons:
        simcoevol_script.append("--relax-constant-sites")
    else:
        simcoevol_script.append("--locus-size " + str(locus_length))
    simcoevol_script.append(new_eco_config_path)
    subprocess.call(" ".join(simcoevol_script), shell=True)

    # Move each replicate into directory, recode alignments, and build xml
    digits = int(math.log10(nreps))+1
    for i in range(0, nreps):
        i = str(i).zfill(digits)
        replicate_dir = os.path.join(out_dir,"replicate-{}".format(i))
        os.mkdir(replicate_dir)

        alignment_path = os.path.join(
            out_dir,
            "simcoevolity-sim-{}-dummy-alignment.nex".format(i))
        new_alignment_path = os.path.join(replicate_dir, "alignment.nex")
        os.rename(alignment_path, new_alignment_path)

        true_value_path = os.path.join(
            out_dir, 
            "simcoevolity-sim-{}-true-values.txt".format(i))
        new_true_value_path = os.path.join(
            replicate_dir, 
            "simulated-true-values.txt")
        os.rename(true_value_path, new_true_value_path)

        rep_config_path = os.path.join(
            out_dir, 
            "simcoevolity-sim-{}-config.yml".format(i))
        new_rep_config_path = os.path.join(
            replicate_dir, 
            "ecoevolity-config.yml")
        os.rename(rep_config_path, new_rep_config_path)
        rep_config = yaml.safe_load(open(new_rep_config_path))
        rep_config["comparisons"][0]["comparison"]["path"] = "../alignment.nex"
        yaml.dump(rep_config, open(new_rep_config_path, "w"))

        generate_xml(
            alignment_path=new_alignment_path,
            out_dir=replicate_dir,
            locus_length=locus_length,
            total_chars=ncharacters,
            pop_size_shape=pop_size_shape,
            pop_size_scale=pop_size_scale,
            birth_rate=wait_time/2, # waiting time = birth rate of 10 * 2 tips
            chain_length=10000000,
            template="starbeast.xml",
            template_dir="templates")

if __name__ == "__main__":
    fire.Fire(gen_data)
