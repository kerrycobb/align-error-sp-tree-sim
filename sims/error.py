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


# def drop_singletons(input, output, prob, rng):
#     dna = AlignIO.read(open(input), "phylip")
#     ids = []
#     seqs = [] 
#     for record in dna:
#        ids.append(record.id)
#        seqs.append(np.array(list(record.seq))) 
#     singletons = 0
#     rmvd_singletons = 0
#     matrix = np.stack(seqs)
#     for column in range(0, matrix.shape[1]):
#         g = np.count_nonzero(matrix[...,column] == "G")
#         t = np.count_nonzero(matrix[...,column] == "T")
#         if g == 1:
#             singletons += 1
#             if rng.uniform(0, 1) >= prob: 
#                 rmvd_singletons += 1
#                 for i in range(0, matrix.shape[0]):
#                     matrix[i, column] = "T"                      
#         elif t == 1:
#             singletons += 1
#             if rng.uniform(0, 1) >= prob:
#                 rmvd_singletons += 1
#                 for i in range(0, matrix.shape[0]):
#                     matrix[i, column] = "G"
#             new_seqs = ["".join(row) for row in matrix]
#             new_seq_records = [] 
#     for id, new_seq in zip(ids, new_seqs):
#         new_seq_records.append(SeqRecord(Seq(new_seq), id=id))
#     new_align = MultipleSeqAlignment(new_seq_records)
#     AlignIO.write(new_align, output, "phylip") 
 
# def gen_error(dir, prob):
#     config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
#     rng = random.Random(config["seed"])
#     parent_dir = pathlib.Path(dir).parent
#     new_dir = os.path.join(parent_dir, "singleton-prob-{}".format(prob))
#     os.mkdir(new_dir) 
#     shutil.copy(os.path.join(dir, "config.yml"), new_dir) 
#     shutil.copy(os.path.join(dir, "eco-config.yml"), new_dir)
#     for rep in range(0, config["nreps"]):
#         rep_dir = "rep-{}".format(rep)
#         rep_dir_path = os.path.join(dir, rep_dir)
#         new_rep_dir_path = os.path.join(new_dir, rep_dir)
#         os.mkdir(new_rep_dir_path)
#         for locus in range(0, config["nloci"]):
#             align_name = "alignment-{}.phy".format(locus) 
#             align_path = os.path.join(rep_dir_path, align_name) 
#             new_align_path = os.path.join(new_rep_dir_path, align_name)
#             drop_singletons(align_path, new_align_path, prob, rng)

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

def phylip_to_matrix(path):
    """
    Convert phylip alignment to numpy matrix
    """
    align = AlignIO.read(open(path), "phylip")
    ids = []
    seqs = [] 
    for record in align:
        ids.append(record.id)
        seqs.append(np.array(record.seq)) 
    matrix = np.stack(seqs)
    return ids, matrix

def matrix_to_phylip(ids, matrix, path): 
    """
    Output numpy matrix of sequence data to phylip
    """
    new_seqs = ["".join(row) for row in matrix]
    new_seq_records = [] 
    for id, new_seq in zip(ids, new_seqs):
        new_seq_records.append(SeqRecord(Seq(new_seq), id=id))
    new_align = MultipleSeqAlignment(new_seq_records)
    AlignIO.write(new_align, path, "phylip")     
 
def gen_error(dir, prob):
    """
    Change singleton site patterns to major allele with probability of 1 - <prob>
    """
    # Get config and make copies in new directories
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    rng = random.Random(config["seed"])
    parent_dir = pathlib.Path(dir).parent
    new_dir = os.path.join(parent_dir, "singleton-prob-{}".format(prob))
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
        for locus in range(0, config["nloci"]):
            # Get alignment name and paths
            align_name = "alignment-{}.phy".format(locus) 
            align_path = os.path.join(rep_dir_path, align_name) 
            new_align_path = os.path.join(new_rep_dir_path, align_name)
            # Read alignment and convert to matrix
            ids, seq_matrix = phylip_to_matrix(align_path)
            # Drop singletons
            new_seq_matrix = drop_singletons(seq_matrix, prob, rng)[0]
            # Output matrix to phylip file 
            matrix_to_phylip(ids, new_seq_matrix, new_align_path)

if __name__ == "__main__":
    fire.Fire(gen_error)