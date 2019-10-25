#! /usr/bin/env python

"""
CLI program for generating a dummy alignment of biallelic characters. 
"""

import sys
import os
import argparse
import logging
import random
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

logging.basicConfig(level=logging.INFO)
_LOG = logging.getLogger(os.path.basename(__file__))
_RNG = random.Random()

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': '0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2016 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}


def arg_is_positive_int(i):
    """
    Returns int if argument is a positive integer; returns argparse error
    otherwise.

    >>> arg_is_positive_int(1) == 1
    True
    """

    try:
        if int(i) < 1:
            raise
    except:
        msg = '{0!r} is not a positive integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)

def check_prefix_and_delimiter(prefix, delimiter):
    if len(delimiter) != 1:
        raise Exception("Delimiter must consist of a one character.")
    if prefix.count(delimiter) > 0:
        raise Exception("Prefix cannot contain the delimiter.")

def write_dummy_biallelic_alignment(
        nspecies = 2,
        ngenomes = 10,
        ncharacters = 1000,
        prefix = "sp",
        delimiter = "-",
        out = sys.stdout,
        prob_missing = 0.0,
        rng = _RNG):
    check_prefix_and_delimiter(prefix, delimiter)
    nspecies_padding = len(str(nspecies))
    ngenomes_padding = len(str(ngenomes))
    label_padding = len(prefix) + nspecies_padding + 1 + ngenomes_padding + 4
    row_idx = 0
    for sp_idx in range(nspecies):
        sp_label = prefix + "{n:0{padding}d}".format(
                n = sp_idx + 1,
                padding = nspecies_padding)
        for g_idx in range(ngenomes):
            row_label = "\'" + sp_label + delimiter + "{n:0{padding}d}".format(
                    n = g_idx + 1,
                    padding = ngenomes_padding) + "\'"
            out.write("{s:{padding}}".format(
                    s = row_label,
                    padding = label_padding))
            if (row_idx % 2 == 0):
                out.write("{0}\n".format("".join((str(0) if rng.random() >= prob_missing else "?" for i in range(ncharacters)))))
            else:
                out.write("{0}\n".format("".join((str(1) if rng.random() >= prob_missing else "?" for i in range(ncharacters)))))
            row_idx += 1

def write_dummy_biallelic_nexus_alignment(
        nspecies = 2,
        ngenomes = 10,
        ncharacters = 1000,
        locus_length = 1,
        prefix = "sp",
        delimiter = "-",
        out = sys.stdout,
        prob_missing = 0.0,
        rng = _RNG):
    check_prefix_and_delimiter(prefix, delimiter)
    out.write("#NEXUS\n")
    out.write("Begin data;\n")
    out.write("    Dimensions ntax={0} nchar={1};\n".format(nspecies * ngenomes, ncharacters))
    out.write("    Format datatype=standard symbols=\"01\" missing=? gap=-;\n")
    out.write("    Matrix\n")
    write_dummy_biallelic_alignment(
        nspecies = nspecies,
        ngenomes = ngenomes,
        ncharacters = ncharacters,
        prefix = prefix,
        delimiter = delimiter,
        out = out,
        prob_missing = prob_missing,
        rng = rng)
    out.write("    ;\n")
    out.write("End;\n")
    if locus_length > 1:
        out.write("\n")
        write_charset_block(
                ncharacters = ncharacters,
                locus_length = locus_length,
                out = out)

def write_charset_block(ncharacters, locus_length,
        out = sys.stdout):
    out.write("Begin sets;\n")
    locus_index = 1
    next_site_index = 1
    while True:
        locus_end_site_index = next_site_index + locus_length - 1
        if locus_end_site_index >= ncharacters:
            out.write("    charset locus{locus_idx}={next_site_idx}-{nchars};\n".format(
                    locus_idx = locus_index,
                    next_site_idx = next_site_index,
                    nchars = ncharacters))
            break
        out.write("    charset locus{locus_idx}={next_site_idx}-{locus_end_idx};\n".format(
                locus_idx = locus_index,
                next_site_idx = next_site_index,
                locus_end_idx = locus_end_site_index))
        next_site_index = locus_end_site_index + 1
        locus_index += 1
    out.write("End;\n")

def main_cli(argv = sys.argv):
    description = '{name} Version {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument('-s', '--nspecies',
            action = 'store',
            type = arg_is_positive_int,
            default = 2,
            help = ('The number of populations.'))
    parser.add_argument('-g', '--ngenomes',
            action = 'store',
            type = arg_is_positive_int,
            default = 10,
            help = ('The number of genomes sampled per population.'))
    parser.add_argument('-c', '--ncharacters',
            action = 'store',
            type = arg_is_positive_int,
            default = 1000,
            help = ('The number of biallelic characters.'))
    parser.add_argument('-l', '--locus-length',
            action = 'store',
            type = arg_is_positive_int,
            default = 1,
            help = ('The length of loci. If greater than 1, a charset block is '
                    'written.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "sp", 
            help = ('Prefix for species labels.'))
    parser.add_argument('-d', '--delimiter',
            action = 'store',
            type = str,
            default = "-", 
            help = ('Symbol used to delimit species label from genome label.'))
    parser.add_argument('-m', '--missing-probability',
            action = 'store',
            type = float,
            default = 0.0,
            help = ('Probability that an allele is missing.'))
    parser.add_argument('--seed',
            action = 'store',
            type = int,
            help = 'Random number seed. Used for random missing data.')

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    _RNG.seed(args.seed)

    write_dummy_biallelic_nexus_alignment(
            nspecies = args.nspecies,
            ngenomes = args.ngenomes,
            ncharacters = args.ncharacters,
            locus_length = args.locus_length,
            prefix = args.prefix,
            delimiter = args.delimiter,
            out = sys.stdout,
            prob_missing = args.missing_probability,
            rng = _RNG)
    

if __name__ == "__main__":
    main_cli()
