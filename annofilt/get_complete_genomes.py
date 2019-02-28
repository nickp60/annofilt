#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Nick Waters
"""

import os
import sys
import glob
import re
import datetime
import subprocess
import argparse
import logging
import shutil
import pandas as pd
import multiprocessing
import copy
import urllib.request


from . import shared_methods as sm

from argparse import Namespace
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_args(): #pragma nocover
    parser = argparse.ArgumentParser(
        description="Given a genus and species, create a pangenome using Roary from complete genomes of the organism",
        add_help=False)
    parser.add_argument(
        "-g", "--genus",
        help="genus")
    parser.add_argument(
        "-s", "--species",
        help="genus")
    parser.add_argument(
        "-n", "--number_of_strains",
        type=int,
        help="genus")
    parser.add_argument(
        "-o",
        "--output",
        help="output dir", required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        "--assembly_summary",
        dest="assembly_summary",
        help="Path to assembly_summary.txt from NCBI; if not present, " +
        "will be downloaded")
    optional.add_argument(
        "-v",
        "--verbosity",
        dest='verbosity', action="store",
        default=2, type=int,
        help="1 = debug(), 2 = info(), 3 = warning(), " +
        "4 = error() and 5 = critical(); " +
        "default: %(default)s")
    optional.add_argument(
        "-h", "--help",
        action="help", default=argparse.SUPPRESS,
        help="Displays this help message")

    args = parser.parse_args()
    return args


def check_exes():
    for exe in ["wget"]:
        if shutil.which(exe) is None:
            raise ValueError("%s executable not found" % exe)

def get_or_check_assembly_metadata(args, logger):
    if args.assembly_summary is not None:
        if not os.path.isfile(args.assembly_summary):
            raise ValueError("Assembly file %s invalid" % args.assembly_summary)
        else:
            return args.ssembly_summary
    else:
        new_path = os.path.join(".", "assembly_summary.txt")
        if not os.path.isfile(new_path):
            logger.info("downloading assembly_summary.txt")
            urllib.request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt", new_path)
        return new_path

def filter_assembly_metadata(args, logger):
    names = "assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material".split("\t")
    d = pd.read_csv(args.assembly_summary, sep="\t", skiprows=(0,1),header=(0),names=names)
    print(d.shape)
    d = d[d.genome_rep=="Full"]
    print(d.shape)
    d = d[(d.assembly_level == "Contig") | (d.assembly_level == "Scaffold")]
    print(d.shape)
    d = d[(d.excluded_from_refseq.notnull())]
    print(d.shape)
    qname = args.genus + " " + args.species
    print(qname)
    d = d[d.organism_name.str.startswith(qname)]
    print(d.shape)
    print(d.genome_rep.head(3))
    return d


def get_n_genomes(d, genomes_dir, n, logger, seed=12345):
    d = d.sample(n, random_state=seed)
    print(str(d.head(1).ftp_path))
    for path in d.ftp_path:
        outpath = os.path.join(genomes_dir, os.path.basename(path) + ".fna.gz")
        if not os.path.isfile(outpath):
            cmd = str("wget {0}/{1}_genomic.fna.gz -O {2}").format(
                path, os.path.basename(path), outpath)
            print(cmd)
            subprocess.run(cmd, shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, check=True)
    unzip_cmd = "gunzip %s*gz" % genomes_dir
    subprocess.run(unzip_cmd, shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, check=True)


def run_roary():
    pass


def main(args=None, logger=None):
    # get args
    if args is None:
        args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    if not os.path.isdir(output_root):
        sys.stderr.write("creating output directory %s\n" % output_root)
        os.makedirs(output_root)
    else:
        sys.stderr.write("Output Directory already exists!\n")
        sys.exit(1)
    check_exes()
    if logger is None:
        logger = sm.set_up_logging(
            outfile=os.path.join(output_root, "make_annofilt_pangenome.log"),
            name="annofilt",
            verbosity=args.verbosity)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    args.assembly_summary = get_or_check_assembly_metadata(args, logger)
    species_df = filter_assembly_metadata(args, logger)
    genomes_dir = os.path.join(output_root, "")
    get_n_genomes(
        d=species_df, genomes_dir=genomes_dir,
        n=args.number_of_strains, logger=logger, seed=12345)
    logger.debug("Done!")



if __name__ == "__main__":
    main()
