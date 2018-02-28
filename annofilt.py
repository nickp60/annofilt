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

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline


logger = logging.getLogger('root')

def get_args():
    parser = argparse.ArgumentParser(
        description="Blast assembly against core genome to find and " +
        "eliminate truncated genes due to bad assembly, " +
        "returning a assembly  with genes that meet a set of criteria. ",
        add_help=False)
    parser.add_argument(
        "ref_gb",
        help="path to .gb or .gbk Genbank genome file")
    parser.add_argument(
        "prokka_dir",
        help="output dir from prokka ")
    parser.add_argument(
        "-o",
        "--output",
        help="output dir ")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        "-r"
        "--reciprocal", dest="reciprocal",
        action="store_true",
        help="reciprocal blast for stringent checking")
    optional.add_argument(
        "-p"
        "--min_percent_id", dest="min_percent_id",
        help="minimum percentage id", type=float, default=.9)
    optional.add_argument(
        "-l"
        "--min_length", dest="min_length",
        help="minimum percentage length", type=float, default=.9)
    optional.add_argument(
        "-e"
        "--min_evalue",
        dest="min_evalue",
        default=.00005,
        help="minimum e value")
    optional.add_argument(
        "-t"
        "--threads",
        dest="threads",
        default=4,
        help="number of threads to use")
    optional.add_argument("-v", "--verbosity",
                        dest='verbosity', action="store",
                        default=2, type=int,
                        help="1 = debug(), 2 = info(), 3 = warning(), " +
                        "4 = error() and 5 = critical(); " +
                        "default: %(default)s")
    args = parser.parse_args()
    return args


def make_prot_prot_blast_cmds(
        query_list, date, evalue, output, threads=1,
        reciprocal=False,subject_file=None, logger=None):
    """given a file, make a blast cmd, and return path to output csv
    """
    assert logger is not None, "must use logging"
    assert isinstance(threads, int), "threads must be integer"
    logger.info("Creating protein BLAST database")
    db_dir = os.path.join(output,
                          os.path.splitext(os.path.basename(subject_file))[0])
    os.makedirs(db_dir, exist_ok=True)
    protdb = os.path.join(db_dir,
                          os.path.splitext(os.path.basename(subject_file))[0])

    setup_blast_db(input_file=subject_file,
                   input_type="fasta",
                   dbtype="prot",
                   out=protdb, logger=logger)
    blast_cmds = []
    blast_outputs = []
    recip_blast_outputs = []
    print(query_list)
    for f in query_list:
        # run forward, nuc aganst prot, blast
        output_path_tab = str(
            os.path.join(output, date) + "_simpleOrtho_results_" +
            os.path.basename(f) + "_vs_protdb.tab")
        blast_cline = NcbiblastpCommandline(query=f,
                                            db=protdb, evalue=evalue,
                                            outfmt=6, out=output_path_tab)
        add_params = " -num_threads {} -num_alignments 20".format(threads)
        blast_command = str(str(blast_cline) + add_params)
        blast_cmds.append(blast_command)
        blast_outputs.append(output_path_tab)
        # run reverse, prot against nuc, blast
        recip_output_path_tab = str(
            os.path.join(output, date) + "_simpleOrtho_results_" +
            "prot_vs_" + os.path.basename(f) + ".tab")
        recip_blast_cline = NcbiblastpCommandline(query=subject_file,
                                                   subject=f,
                                                   evalue=evalue,
                                                   outfmt=6, out=recip_output_path_tab)
        recip_blast_command = str(str(recip_blast_cline) + add_params)
        if reciprocal:
            blast_cmds.append(recip_blast_command)
            recip_blast_outputs.append(recip_output_path_tab)
    return(blast_cmds, blast_outputs, recip_blast_outputs)



def filter_BLAST_df(df1, df2, min_length_percent, min_percent, reciprocal, logger=None):
    """ results from pd.read_csv with default BLAST output 6 columns
    df1 must be genomes against genes, and df2 must be genes against genomes,
    because we have to split the names so all all the contigs are recognized
    as coming from one genome.  returns a df
    """
    assert logger is not None, "must use a logger"
    logger.debug("shape of blast results")
    df1['genome'] = df1.query_id.str.split('_').str.get(0)
    logger.debug(df1.shape)

    if args.reciprocal:
        logger.debug("shape of recip blast results")
        logger.debug(df2.shape)
        df2['genome'] = df2.subject_id.str.split('_').str.get(0)


    if not reciprocal:
        filtered = pd.DataFrame(columns=df1.columns)
        unq_subject = df1.subject_id.unique()
        unq_query = df1.genome.unique()
        recip_hits = []
        nonrecip_hits = []  # should be renamed to nogood_hits

        for gene in unq_subject:
            for genome in unq_query:
                logger.debug("Checking %s in %s for reciprocity" % (gene, genome))
                tempdf1 = df1.loc[(df1["subject_id"] == gene) &
                              (df1["genome"] == genome), ]
                # here we get the best hit (the one maximizing e value)
                subset1 = tempdf1.loc[
                    (tempdf1["identity_perc"] > min_percent) &
                    (tempdf1["bit_score"] == tempdf1["bit_score"].max())]
                logger.debug("grouped df shape: ")
                logger.debug(tempdf1.shape)
                if subset1.empty:
                    logger.info("No full hits for %s in %s", gene, genome)
                    logger.debug(tempdf1)
                    nonrecip_hits.append([gene, genome])
                else:
                    if all(subset1["alignment_length"] > subset1["subject_length"] * min_length_percent):
                        filtered = filtered.append(subset1)
            # logger.debug(subset.shape)
        return(filtered)


    # recip structure
    filtered = pd.DataFrame(columns=df1.columns)
    unq_subject = df1.subject_id.unique()
    unq_query = df1.genome.unique()
    recip_hits = []
    nonrecip_hits = []
    for gene in unq_subject:
        for genome in unq_query:
            logger.debug("Checking %s in %s for reciprocity" % (gene, genome))
            tempdf1 = df1.loc[(df1["subject_id"] == gene) &
                              (df1["genome"] == genome), ]
            tempdf2 = df2.loc[(df2["query_id"] == gene) &
                              (df2["genome"] == genome), ]
            if tempdf1.empty or tempdf2.empty:
                logger.info("skipping %s in %s; no match found", gene, genome)
            else:
                subset1 = tempdf1.loc[
                    (tempdf1["identity_perc"] > min_percent) &
                    (tempdf1["bit_score"] == tempdf1["bit_score"].max())]
                # (tempdf1["alignement_l"] == tempdf1["bit_score"].max())]
                subset2 = tempdf2.loc[
                    (tempdf2["identity_perc"] > min_percent) &
                    (tempdf2["bit_score"] == tempdf2["bit_score"].max())]
                logger.debug("grouped df shape: ")
                logger.debug(tempdf1.shape)
                logger.debug("grouped df2 shape: " )
                logger.debug(tempdf2.shape)
                if subset1.empty or subset2.empty:
                    logger.info("No reciprocol hits for %s in %s", gene, genome)
                    logger.debug(tempdf1)
                    logger.debug(tempdf2)
                    nonrecip_hits.append([gene, genome])
                else:
                    # logger.debug(tempdf1)
                    # logger.debug("tempdf2")
                    # logger.debug(tempdf2)
                    # logger.debug("subset1")
                    # logger.debug(subset1)
                    # logger.debug("subset2")
                    # logger.debug(subset2)
                    if subset1.iloc[0]["query_id"] == subset2.iloc[0]["subject_id"]:
                        recip_hits.append([gene, genome])
                        filtered = filtered.append(subset1)
                        logger.info("Reciprocol hits for %s in %s!", gene, genome)
                    else:
                        nonrecip_hits.append([gene, genome])
                        logger.info("No reciprocol hits for %s in %s", gene, genome)

            # logger.debug(subset.shape)
    logger.debug("Non-reciprocal genes:")
    logger.debug(nonrecip_hits)
    logger.debug("Reciprocal genes:")
    logger.debug(recip_hits)
    logger.debug("filtered shape:")
    logger.debug(filtered.shape)
    return(filtered)


def write_pipe_extract_cmds(df, outfile, logger=None):
    #% parse output
    assert logger is not None, "must use a logger"
    logger.debug("cleaning up the csv output")
    with open(outfile, "a") as outf:
        for index, row in df.iterrows():
            if row['q_start'] > row['q_end']:
                logger.debug("hit is on the (-) strand")
                line = "{0}-RC@{1} :{2}:{3}".format(
                    row['subject_id'],
                    row['query_id'],
                    int(row['q_end']),
                    int(row['q_start']))
            else:
                line = "{0}@{1} :{2}:{3}".format(
                    row['subject_id'],
                    row['query_id'],
                    int(row['q_start']),
                    int(row['q_end']))
            sys.stdout.write(line + "\n")
            outf.write(line + "\n")


def set_up_logging(verbosity, outfile, name):
    """
    Set up logging a la pyani, with
    a little help from:
    https://aykutakin.wordpress.com/2013/08/06/
        logging-to-console-and-file-in-python/
    requires logging, os, sys, time
    logs debug level to file, and [verbosity] level to stderr
    return a logger object
    """
    if (verbosity * 10) not in range(10, 60, 10):
        raise ValueError('Invalid log level: %s' % verbosity)
    # logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to given verbosity
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(level=(verbosity * 10))
    console_err_format = logging.Formatter(str("%(asctime)s - " +
                                               "%(levelname)s - %(message)s"),
                                           "%Y-%m-%d %H:%M:%S")
    console_err.setFormatter(console_err_format)
    logger.addHandler(console_err)
    # create debug file handler and set level to debug
    try:
        logfile_handler = logging.FileHandler(outfile, "w")
        logfile_handler.setLevel(logging.DEBUG)
        logfile_handler_formatter = \
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        logfile_handler.setFormatter(logfile_handler_formatter)
        logger.addHandler(logfile_handler)
    except:
        logger.error("Could not open {0} for logging".format(outfile))
        sys.exit(1)
    logger.debug("Initializing logger")
    logger.debug("logging at level {0}".format(verbosity))
    return logger


def setup_blast_db(input_file, input_type="fasta", dbtype="prot",
                   title="blastdb", out="blastdb",
                   makeblastdb_exe='', logger=None):
    """
    This runs make blast db with the given parameters
    requires logging, os, subprocess, shutil
    """
    if makeblastdb_exe == '':
        makeblastdb_exe = shutil.which("makeblastdb")
        if logger:
            logger.info("makeblastdb executable: %s", makeblastdb_exe)
    makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
                    "-out {4}").format(makeblastdb_exe,
                                       input_file,
                                       input_type,
                                       dbtype,
                                       out)
    if logger:
        logger.info("Making blast db: {0}".format(makedbcmd))
    try:
        subprocess.run(makedbcmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        logging.debug("BLAST database '{0}' created here: {1}".format(
            title, out))
        return 0
    except:
        if logger:
            logging.error("Something bad happened when trying to make " +
                          "a blast database")
        sys.exit(1)


def merge_outfiles(filelist, outfile_name):
    """
    """
    # only grab .tab files, ie, the blast output
    filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
    if len(filelist) == 1:
        # print("only one file found! no merging needed")
        return(filelist[0])
    else:
        # print("merging all the blast results to %s" % outfile_name)
        nfiles = len(filelist)
        fout = open(outfile_name, "a")
        # first file:
        for line in open(filelist[0]):
            fout.write(line)
        #  now the rest:
        for num in range(1, nfiles):
            f = open(filelist[num])
            for line in f:
                fout.write(line)
            f.close()  # not really needed
        fout.close()
    return(outfile_name)


def BLAST_tab_to_df(path):
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"]
    raw_csv_results = pd.read_csv(
        open(path), comment="#", sep="\t", names=colnames)
    return raw_csv_results


def make_new_genbank(genbank, new_genbank, approved_accessions, logger):
    newrecs = []
    with open(genbank, "r") as inf:
        for rec in SeqIO.parse(inf, "genbank"):
            newrec = copy.deepcopy(rec)
            newrec.features = []
            for idx, feat in enumerate(rec.features):
                if feat.type != "CDS":
                    newrec.features.append(feat)
                    continue
                if feat.qualifiers.get("locus_tag")[0] in approved_accessions:
                    newrec.features.append(feat)
                else:
                    logger.info(feat.qualifiers.get("locus_tag")[0] +
                                " was rejected")
                newrecs.append(new_rec)
    with open(new_genbank, "w") as outf:
        for nr in newrecs:
            SeqIO.write(nr, outf, "genbank")

if __name__ == "__main__":
    # get args
    args = get_args()
    output_root = os.path.abspath(os.path.expanduser(args.output))
    if not os.path.isdir(output_root):
        sys.stderr.write("creating output directory %s\n" % output_root)
        os.makedirs(output_root)
    else:
        sys.stderr.write("Output Directory already exists!\n")
        sys.exit(1)
    logger = set_up_logging(outfile=os.path.join(output_root, "log.log"),
                            name="annofilt", verbosity=args.verbosity)
    print(logger)
    # verify imputes
    # build temp blast DB from
    # Read in Genbank
    # for each gene:
    #     blast agianst the database
    #     if no match:
    #         if args.keep_nohits:
    #             keep Annotation
    #         else:
    #             ditch that annotation
    #     if the best match passes thresholds:
    #         keep Annotation
    #     else:
    #         ditch that annotation
    # logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
    #                     format='%(name)s (%(levelname)s): %(message)s')
    # logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    date = str(datetime.datetime.now().strftime('%Y%m%d'))

    query_gb = glob.glob(args.prokka_dir + "*.faa")
    if len(query_gb) is None:
        raise ValueError("no .faa file found in prokka dir!")
    if len(query_gb) > 1:
        raise ValueError("Multiple AA fastas files found in prokka ouput!")
    query_gb = query_gb[0]
    commands, paths_to_outputs, paths_to_recip_outputs = \
        make_prot_prot_blast_cmds(
            query_list=query_gb,
            evalue=args.min_evalue,
            reciprocal=args.reciprocal,
            subject_file=args.ref_gb,
            threads=args.threads,
            output=output_root, date=date,
            logger=logger)
    logger.debug("cmds:" )
    logger.debug(commands)
    logger.debug("paths to outputs:" )
    logger.debug(paths_to_outputs )
    logger.debug("paths to recip outputs:" )
    logger.debug(paths_to_recip_outputs )
    # check for existing blast results

    if not all([os.path.isfile(x) for x in paths_to_outputs]):
        pool = multiprocessing.Pool()
        logger.debug("Running the following commands in parallel " +
                     "(this could take a while):")
        logger.debug("\n" + "\n".join([x for x in commands]))
        results = [
            pool.apply_async(subprocess.run,
                             (cmd,),
                             {"shell": sys.platform != "win32",
                              "stdout": subprocess.PIPE,
                              "stderr": subprocess.PIPE,
                              "check": True})
            for cmd in commands]
        pool.close()
        pool.join()
        reslist = []
        reslist.append([r.get() for r in results])
    else:
        pass
    merged_tab = os.path.join(output_root,
                              "merged_results.tab")
    recip_merged_tab = os.path.join(output_root,
                                    "recip_merged_results.tab")
    post_merged_tab = merge_outfiles(filelist=paths_to_outputs,
                                            outfile_name=merged_tab)
    resultsdf = BLAST_tab_to_df(post_merged_tab)
    if args.reciprocal:
        recip_resultsdf = BLAST_tab_to_df(post_recip_merged_tab)
        post_recip_merged_tab = merge_outfiles(filelist=paths_to_recip_outputs,
                                               outfile_name=recip_merged_tab)
    else:
        recip_resultsdf = None
    # go through the subject file to aretrieve the ideal lengths.
    with open(args.ref_gb, "r") as inf:
        resultsdf["subject_length"] = 0
        for rec in SeqIO.parse(inf, "fasta"):
            resultsdf.ix[resultsdf.subject_id==rec.id, 'subject_length'] = len(rec.seq)
    # sanity check that this worked
    if all(resultsdf["subject_length"] == 0):
        raise ValueError(
            "error matching blast gene names to those parsed by " +
            "BioPython!  Until this is sorted, we wont be able to " +
            "match up the total expected gene length. Exiting")

    filtered_hits = filter_BLAST_df(
        df1=resultsdf,
        df2=recip_resultsdf,
        reciprocal=args.reciprocal,
        min_percent=args.min_percent_id,
        min_length_percent=args.min_length,
        logger=logger)
    print(filtered_hits)
    make_new_genbank(
        genbank=query_gb,
        new_genbank=os.path.splitext(query_gb)[0] + "_filtered.gbk",
        approved_accessions=filtered_hits["query_id"],
        logger=logger)
    filtered_hits.to_csv(
        os.path.join(output_root, "simpleOrtho_filtered_hits.csv"))
    write_pipe_extract_cmds(
        outfile=os.path.join(output_root,
                             "simpleOrtho_regions.txt"),
        df=filtered_hits, logger=logger)
