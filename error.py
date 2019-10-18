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


def drop_singletons(input, output, prob, rng):
    dna = AlignIO.read(open(input), "phylip")
    ids = []
    seqs = [] 
    for record in dna:
       ids.append(record.id)
       seqs.append(np.array(list(record.seq))) 
    singletons = 0
    rmvd_singletons = 0
    matrix = np.stack(seqs)
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
            new_seqs = ["".join(row) for row in matrix]
            new_seq_records = [] 
            for id, new_seq in zip(ids, new_seqs):
                new_seq_records.append(SeqRecord(Seq(new_seq), id=id))
            new_align = MultipleSeqAlignment(new_seq_records)
            AlignIO.write(new_align, output, "phylip") 
 

def gen_error(dir, prob):
    config = yaml.safe_load(open(os.path.join(dir, "config.yml")))
    rng = random.Random(config["seed"])
    parent_dir = pathlib.Path(dir).parent
    new_dir = os.path.join(parent_dir, "singleton-prob-{}".format(prob))
    os.mkdir(new_dir) 
    shutil.copy(os.path.join(dir, "config.yml"), new_dir) 
    shutil.copy(os.path.join(dir, "eco-config.yml"), new_dir)
    for rep in range(0, config["nreps"]):
        rep_dir = "rep-{}".format(rep)
        rep_dir_path = os.path.join(dir, rep_dir)
        new_rep_dir_path = os.path.join(new_dir, rep_dir)
        os.mkdir(new_rep_dir_path)
        for locus in range(0, config["nloci"]):
            align_name = "alignment-{}.phy".format(locus) 
            align_path = os.path.join(rep_dir_path, align_name) 
            new_align_path = os.path.join(new_rep_dir_path, align_name)
    #         drop_singletons(align_path, new_align_path, prob, rng)


# def gen_error(dir, rng, config):
    # for prob in config["singleton_sample_prob"]:
    #     out_dir = os.path.join(dir, "singleton-prob-{}".format(prob))
    #     os.mkdir(out_dir)
    #     singletons = 0
    #     rmvd_singletons = 0       
    #     for locus in range(0, config["nloci"]):
    #         align_path = os.path.join(
    #             dir, 
    #             "singleton-prob-1.0", 
    #             "alignment-{}.phy".format(locus))
    #         dna = AlignIO.read(open(align_path), "phylip")
    #         ids = []
    #         seqs = [] 
    #         for record in dna:
    #            ids.append(record.id)
    #            seqs.append(np.array(list(record.seq))) 
    #         matrix = np.stack(seqs)
    #         # print(matrix)
    #         for column in range(0, matrix.shape[1]):
    #             g = np.count_nonzero(matrix[...,column] == "G")
    #             t = np.count_nonzero(matrix[...,column] == "T")
    #             if g == 1:
    #                 singletons += 1
    #                 if rng.uniform(0, 1) >= prob: 
    #                     rmvd_singletons += 1
    #                     for i in range(0, matrix.shape[0]):
    #                         matrix[i, column] = "T"                      
    #             elif t == 1:
    #                 singletons += 1
    #                 if rng.uniform(0, 1) >= prob:
    #                     rmvd_singletons += 1
    #                     for i in range(0, matrix.shape[0]):
    #                         matrix[i, column] = "G"
    #         # print(matrix)
    #         new_seqs = ["".join(row) for row in matrix]
    #         new_seq_records = [] 
    #         for id, new_seq in zip(ids, new_seqs):
    #             new_seq_records.append(SeqRecord(Seq(new_seq), id=id))
    #         new_align = MultipleSeqAlignment(new_seq_records)
    #         AlignIO.write(
    #             new_align, 
    #             os.path.join(out_dir, "alignment-{}.phy".format(locus)), 
    #             "phylip") 
    #     stats = dict(singletons=singletons, removed_singletons=rmvd_singletons)
    #     with open(os.path.join(out_dir, "stats.yaml"), "w") as fh:
    #         yaml.dump(stats, fh)

# def gen_error(dir, prob):
    
    # for prob in config["singleton_sample_prob"]:
    #     out_dir = os.path.join(dir, "singleton-prob-{}".format(prob))
    #     os.mkdir(out_dir)
    #     singletons = 0
    #     rmvd_singletons = 0       
    #     for locus in range(0, config["nloci"]):
    #         align_path = os.path.join(
    #             dir, 
    #             "singleton-prob-1.0", 
    #             "alignment-{}.phy".format(locus))
    #         dna = AlignIO.read(open(align_path), "phylip")
    #         ids = []
    #         seqs = [] 
    #         for record in dna:
    #            ids.append(record.id)
    #            seqs.append(np.array(list(record.seq))) 
    #         matrix = np.stack(seqs)
    #         # print(matrix)
    #         for column in range(0, matrix.shape[1]):
    #             g = np.count_nonzero(matrix[...,column] == "G")
    #             t = np.count_nonzero(matrix[...,column] == "T")
    #             if g == 1:
    #                 singletons += 1
    #                 if rng.uniform(0, 1) >= prob: 
    #                     rmvd_singletons += 1
    #                     for i in range(0, matrix.shape[0]):
    #                         matrix[i, column] = "T"                      
    #             elif t == 1:
    #                 singletons += 1
    #                 if rng.uniform(0, 1) >= prob:
    #                     rmvd_singletons += 1
    #                     for i in range(0, matrix.shape[0]):
    #                         matrix[i, column] = "G"
    #         # print(matrix)
    #         new_seqs = ["".join(row) for row in matrix]
    #         new_seq_records = [] 
    #         for id, new_seq in zip(ids, new_seqs):
    #             new_seq_records.append(SeqRecord(Seq(new_seq), id=id))
    #         new_align = MultipleSeqAlignment(new_seq_records)
    #         AlignIO.write(
    #             new_align, 
    #             os.path.join(out_dir, "alignment-{}.phy".format(locus)), 
    #             "phylip") 
    #     stats = dict(singletons=singletons, removed_singletons=rmvd_singletons)
    #     with open(os.path.join(out_dir, "stats.yaml"), "w") as fh:
    #         yaml.dump(stats, fh)

if __name__ == "__main__":
    fire.Fire(gen_error)