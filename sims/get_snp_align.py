#!/usr/bin/env python

import os
import shutil
import fire
import yaml
import pathlib
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

def align_to_matrix(align):
    """
    Convert phylip alignment to numpy matrix
    """
    ids = []
    seqs = []
    for record in align:
        ids.append(record.id)
        seqs.append(np.array(record.seq)) 
    matrix = np.stack(seqs)
    return ids, matrix

def matrix_to_align(ids, matrix): 
    """
    Output numpy matrix of sequence data to phylip
    """
    new_seqs = ["".join(row) for row in matrix]
    new_seq_records = [] 
    for id, new_seq in zip(ids, new_seqs):
        new_seq_records.append(SeqRecord(Seq(new_seq), id=id))
    new_align = MultipleSeqAlignment(new_seq_records)
    return new_align

def get_snps(matrix, single=True):
    rows = matrix.shape[0]
    polymorphic_sites = 0
    snp_cols = [] 
    for column in range(0, matrix.shape[1]):
        g = np.count_nonzero(matrix[...,column] == "G")
        if g != 0 and g != rows: 
            snp_cols.append(matrix[..., column])
            polymorphic_sites += 1
            if single:
                break 
    return snp_cols, polymorphic_sites

def get_snp_align(dir, single=False, overwrite=False):
    """
    Creat alignment of SNPs only, use --single to only select a single snp 
    """
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    parent_dir = pathlib.Path(dir).parent
    if single:
        new_dir = os.path.join(parent_dir, os.path.basename(dir) + "-snp")
    else:
        new_dir = os.path.join(os.path.dirname(dir), os.path.basename(dir) + "-all-snps")
    if overwrite:
        try:
            shutil.rmtree(new_dir) 
        except:
            pass
    os.mkdir(new_dir) 
    shutil.copy(os.path.join(dir, "config.yml"), new_dir) 
    eco_config = yaml.safe_load(open(os.path.join(dir, "eco-config.yml")))
    eco_config["global_comparison_settings"]["constant_sites_removed"] = True
    yaml.dump(eco_config, open(os.path.join(new_dir, "eco-config.yml"), "w"))
    for rep in range(0, config["nreps"]):
        rep_dir = "rep-{}".format(rep)
        rep_dir_path = os.path.join(dir, rep_dir)
        new_rep_dir_path = os.path.join(new_dir, rep_dir)
        os.mkdir(new_rep_dir_path)
        shutil.copy(os.path.join(rep_dir_path, "species_tree.nex"), 
                new_rep_dir_path)
        new_cols = []
        for locus in range(0, config["nloci"]):
            align_name = "alignment-{}.phy".format(locus) 
            align_path = os.path.join(rep_dir_path, align_name) 
            align = AlignIO.read(open(align_path), "phylip")
            length = align.get_alignment_length()
            ids, align_matrix = align_to_matrix(align)
            snp_cols = get_snps(align_matrix, single=single)[0]
            new_cols.extend(snp_cols)
        concat_matrix = np.column_stack(new_cols) 
        new_align_path = os.path.join(new_rep_dir_path, "snp-alignment.phy")
        new_align = matrix_to_align(ids, concat_matrix)
        AlignIO.write(new_align, open(new_align_path, "w"), "phylip")
    print("Finished getting SNP alignemnt from {}".format(dir))

if __name__ == "__main__":
    fire.Fire(get_snp_align)