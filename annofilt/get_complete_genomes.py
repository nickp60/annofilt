#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Nick Waters
"""

import os
import sys
import subprocess
import argparse
import shutil
import pandas as pd
import urllib.request


from . import shared_methods as sm


def get_args():  # pragma nocover
    parser = argparse.ArgumentParser(
        description="Given a genus and species, download complete " +
        "genomes from NCBI",
        add_help=False)
    parser.add_argument(
        "-g", "--genus",
        help="genus", required=True)
    parser.add_argument(
        "-s", "--species",
        help="genus", required=True)
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
        "--max_contigs",
        dest="max_contigs",
        default=10,
        help="Maximum number of contigs tolerated. Many fragmented " +
        " genome assemblies are categorized as 'complete' by ncbi; " +
        "so if you expect at most 1 chromosome and 10 plasmids, " +
        "set --max_contigs to 11. default: %(default)s")
    optional.add_argument(
        "--max_strains",
        dest="max_strains",
        help="how many genomes to try to download to " +
        " reach the --number_of_strains")
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
            raise ValueError("Metadata file %s invalid" %
                             args.assembly_summary)
        else:
            return args.assembly_summary
    else:
        new_path = os.path.join(".", "assembly_summary.txt")
        if not os.path.isfile(new_path):
            logger.info("downloading assembly_summary.txt")
            urllib.request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/" +
                                       "genomes/genbank/bacteria/" +
                                       "assembly_summary.txt", new_path)
        return new_path


def filter_assembly_metadata(args, logger):
    names = "assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material".split("\t")
    d = pd.read_csv(args.assembly_summary, sep="\t",
                    skiprows=(0, 1), header=(0), names=names)
    # print(d.shape)
    d = d[d.genome_rep == "Full"]
    # print(d.shape)
    d = d[(d.assembly_level == "Complete Genome") |
          (d.assembly_level == "Chromosome")]
    # print(d.shape)
    d = d[(d.excluded_from_refseq.notnull())]
    # print(d.shape)
    qname = args.genus + " " + args.species
    # print(qname)
    d = d[d.organism_name.str.startswith(qname)]
    # print(d.shape)
    # print(d.genome_rep.head(3))
    return d


def get_genome(path, genomes_dir, logger):
    outpath = os.path.join(genomes_dir, os.path.basename(path) + ".fna.gz")
    if not os.path.isfile(outpath):
        cmd = str("wget {0}/{1}_genomic.fna.gz -O {2}").format(
            path, os.path.basename(path), outpath)
        logger.debug("Getting assembly " + path)
        subprocess.run(cmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
    return outpath


def main(args=None, logger=None):
    # get args
    if args is None:
        args = get_args()
    if args.max_strains is None:
        args.max_strains = args.number_of_strains + 10
    else:
        if args.max_strains < args.number_of_strains:
            raise ValueError("Maximum number of strains to try cannot be " +
                             "less than the number of genomes")
    output_root = os.path.abspath(os.path.expanduser(args.output))
    if not os.path.isdir(output_root):
        sys.stderr.write("creating output directory %s\n" % output_root)
        os.makedirs(output_root)
    else:
        sys.stderr.write("WARNING! Using existing output directory!\n")
    check_exes()
    if logger is None:
        logger = sm.set_up_logging(
            outfile=os.path.join(output_root, "make_annofilt_pangenome.log"),
            name="annofilt",
            verbosity=args.verbosity)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    logger.debug("sorting out the assembly_metadata.txt file")
    args.assembly_summary = get_or_check_assembly_metadata(args, logger)
    logger.debug("filtering assembly metadata")
    species_df = filter_assembly_metadata(args, logger)
    genomes_dir = os.path.join(output_root, "")
    not_enough_strains = True
    not_reached_max_tries = True
    counter = 0
    n_successful = 0
    # shuffle the order of the strains we have left after filtering
    species_df = species_df.sample(
        frac=1, random_state=12345).reset_index(drop=True)
    logger.debug("Attempting to download strains")
    if species_df.shape[0] < args.number_of_strains:
        logger.warning("Only %i elegible strains available!" %
                       species_df.shape[0])
        args.number_of_strains = species_df.shape[0]
    while not_enough_strains and not_reached_max_tries:
        logger.debug("getting strain %i of %i" %
                     (n_successful+1, args.number_of_strains))
        this_path = get_genome(
            path=species_df.ftp_path[counter], genomes_dir=genomes_dir,
            logger=logger)
        logger.debug("    unzipping %s" % this_path)
        # -f force overwrite (on osx, -o overwrites like with `unzip`, but not linux?)
        unzip_cmd = "gunzip -f %s" % this_path
        this_path = this_path.replace(".gz", "")

        subprocess.run(unzip_cmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)

        logger.debug("    checking number of contigs")
        grep_cmd = "grep '>' %s" % this_path
        result = subprocess.run(grep_cmd, shell=sys.platform != "win32",
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, check=True)
        n_contigs = len(result.stdout.decode("utf-8").split("\n"))
        if n_contigs > args.max_contigs:
            logger.info("rejecting %s: %i contigs" % (this_path, n_contigs))
            os.remove(this_path)
        else:
            n_successful = n_successful + 1
        # status update
        if n_successful == args.number_of_strains:
            not_enough_strains = False
        if counter == args.max_strains:
            not_reached_max_tries = False
        counter = counter + 1
    logger.debug("Saved {} strains!".format(n_successful))
    logger.debug("Done!")


if __name__ == "__main__":
    main()
