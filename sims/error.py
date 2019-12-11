#!/usr/bin/env python

import os
import shutil
import fire
import glob
import yaml
import pathlib
import numpy as np
import random
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna
import itertools as it

def drop_singletons(matrix, prob, rng):
    """
    Change singleton site patterns to major allele with probability of 1 - <prob>
    """
    singletons = 0
    rmvd_singletons = 0
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
    return matrix, (singletons, rmvd_singletons)

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

def drop_het(align, prob, rng):
    lambda_func = lambda x: x.id.split("_")[0]
    align._records.sort(key=lambda_func)
    for key, group in it.groupby(align._records, lambda_func):
        group_records = list(group)
        rng.shuffle(group_records)
        group_records_iter = iter(group_records)
        for record0 in group_records_iter:
            record1 = group_records_iter.__next__()
            if prob < rng.random():
                record0.seq = record1.seq
    return align

def gen_error(dir, prob, single=False, het=False, concat=False):
    """
    Change singleton site patterns to major allele with probability of 1 - <prob>
    Or drop one haplotype for randomly chosen pairs within populations with probability of 1 - <prob>
    """
    if not single and not het:
          quit("Must use either --single or --het argument")
      
    # Get config and make copies in new directories
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    rng = random.Random(config["seed"])
    parent_dir = pathlib.Path(dir).parent
    if single:
        new_dir = os.path.join(parent_dir, "singleton-prob-{}".format(prob))

    elif het:
        new_dir = os.path.join(parent_dir, "het-prob-{}".format(prob))
    os.mkdir(new_dir) 
    shutil.copy(os.path.join(dir, "config.yml"), new_dir) 
    shutil.copy(os.path.join(dir, "eco-config.yml"), new_dir)
    for rep in range(0, config["nreps"]):
        # Create new directory for new alignmetns 
        rep_dir = "rep-{}".format(rep)
        rep_dir_path = os.path.join(dir, rep_dir)
        new_rep_dir_path = os.path.join(new_dir, rep_dir)
        os.mkdir(new_rep_dir_path)
        shutil.copy(os.path.join(rep_dir_path, "species_tree.nex"), 
                new_rep_dir_path)
        # If single concatenated alignment
        if concat:
            align_path = os.path.join(rep_dir_path, "alignment.phy")
            new_align_path = os.path.join(new_rep_dir_path, "alignment.phy")
            align = AlignIO.read(open(align_path), "phylip-relaxed")
            # Drop singletons
            if single:
                ids, seq_matrix = align_to_matrix(align)
                new_seq_matrix = drop_singletons(seq_matrix, prob, rng)[0]
                new_align = matrix_to_align(ids, new_seq_matrix)
            # Drop hets
            elif het:
                new_align = drop_het(align, prob, rng)
            AlignIO.write(new_align, open(new_align_path, "w"), "phylip-relaxed")
        # If separate alignment for each locus
        else:
            for locus in range(0, config["nloci"]):
                # Get alignment name and paths
                align_name = "alignment-{}.phy".format(locus) 
                align_path = os.path.join(rep_dir_path, align_name) 
                new_align_path = os.path.join(new_rep_dir_path, align_name)
                align = AlignIO.read(open(align_path), "phylip-relaxed")
                # Drop singletons
                if single:
                    ids, seq_matrix = align_to_matrix(align)
                    new_seq_matrix = drop_singletons(seq_matrix, prob, rng)[0]
                    new_align = matrix_to_align(ids, new_seq_matrix)
                # Drop hets
                elif het:
                    new_align = drop_het(align, prob, rng)
                AlignIO.write(new_align, open(new_align_path, "w"), "phylip-relaxed")

if __name__ == "__main__":
    fire.Fire(gen_error)