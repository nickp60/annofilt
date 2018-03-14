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

from argparse import Namespace
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
        "reference",
        help="path to an fasta (either protein or nucleotide) " +
        "containing al the genes in the  genome;  for instance, " +
        "the pan_genome_reference.fa from roary")
    parser.add_argument(
        "prokka_dir",
        help="output dir from prokka ")
    parser.add_argument(
        "-o",
        "--output",
        help="output dir", required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        "--full", dest="full",
        action="store_true",
        help="check ALL genes for completeness, not just on ends of contigs." +
        " This is much slower.")
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
    optional.add_argument(
        "--genes", dest="just_genes",
        action="store_true",
        help="only deal with the genes; default is to deal with " +
        "anything that has a locus tag")
    optional.add_argument(
        "-v", "--verbosity",
        dest='verbosity', action="store",
        default=2, type=int,
        help="1 = debug(), 2 = info(), 3 = warning(), " +
        "4 = error() and 5 = critical(); " +
        "default: %(default)s")

    args = parser.parse_args()
    return args


def make_prokka_files_object(prokka_dir):
    """
    """
    checked = []
    for ext in ["faa", "gbk", "gff"]:
        files = glob.glob(prokka_dir + "*" + ext)
        if len(files) is None:
            raise ValueError("No " + ext + " file found in prokka dir!")
        if len(files) > 1:
            raise ValueError(
                "Multiple " + ext + " files found in prokka ouput!")
        checked.append(files[0])
    return Namespace(faa=checked[0],
                     gbk=checked[1],
                     gff=checked[2])


def return_list_of_locus_tags(gbk=None, faa=None):
    lt_list = []
    if gbk is not None:
        with open(gbk, "r") as inf:
            for rec in SeqIO.parse(inf, "genbank"):
                for feat in rec.features:
                    lt = feat.qualifiers.get("locus_tag")
                    if lt is not None:
                        # I dont know why locus tags are stored as lists;
                        # never seen one with more than one
                        lt_list.append(lt[0])
    else:
        assert faa is not None, "must submit either a genbank or AA fasta"
        with open(faa, "r") as inf:
            for rec in SeqIO.parse(inf, "fasta"):
                lt_list.append(rec.id)
    # we return a list after getting the unique entries, as there is
    # usually a CDS and gene/tRNA/rRNA entry for each thing
    return list(set(lt_list))


def make_prot_prot_blast_cmds(
        query_file, date, evalue, output, threads=1,
        reciprocal=False,subject_file=None, logger=None):
    """given a file, make a blast cmd, and return path to output csv
    This should handle both protein and nucleotide references.
    """
    assert logger is not None, "must use logging"
    assert isinstance(threads, int), "threads must be integer"
    # if dir exists, we already have the db created
    #  (as we make it in the output dir)
    db_dir = os.path.join(output,
                          os.path.splitext(os.path.basename(subject_file))[0])
    subjectdb = os.path.join(db_dir,
                             os.path.splitext(os.path.basename(subject_file))[0])
    first_record = SeqIO.parse(subject_file, "fasta").__next__()
    subject_is_protein = False

    if not os.path.isdir(db_dir):
        logger.info("Creating protein BLAST database")
        os.makedirs(db_dir, exist_ok=True)

        if first_record.seq.alphabet == IUPAC.IUPACProtein():
            subject_is_protein = True
            setup_blast_db(input_file=subject_file,
                           input_type="fasta",
                           dbtype="prot",
                           out=subjectdb, logger=logger)
        else:
            setup_blast_db(input_file=subject_file,
                           input_type="fasta",
                           dbtype="nucl",
                           out=subjectdb, logger=logger)
    blast_cmds = []
    blast_outputs = []
    recip_blast_outputs = []
    for f in [query_file]:
        qname = os.path.splitext(os.path.basename(f))[0]
        # run forward, nuc aganst prot, blast
        output_path_tab = os.path.join(
            output,
            qname + "_vs_protdb.tab")
        if subject_is_protein:
            blast_cline = NcbiblastpCommandline(
                query=f,
                db=subjectdb, evalue=evalue, out=output_path_tab)
        else:
            blast_cline = NcbitblastnCommandline(
                query=f,
                db=subjectdb, evalue=evalue, out=output_path_tab)
        add_params = str(
            " -num_threads {} -num_alignments 20 " +
            "-outfmt '6 qaccver saccver pident length mismatch " +
            "gapopen qstart qend sstart send evalue bitscore slen'"
        ).format(threads)
        blast_command = str(str(blast_cline) + add_params)
        blast_cmds.append(blast_command)
        blast_outputs.append(output_path_tab)
        # run reverse, prot against nuc, blast
        recip_output_path_tab = os.path.join(
            output,
            "protdb_vs_" + qname + ".tab")
        if subject_is_protein:
            recip_blast_cline = NcbiblastpCommandline(
                query=subject_file,
                subject=f,
                evalue=evalue,
                out=recip_output_path_tab)
        else:
            recip_blast_cline = NcbiblastxCommandline(
                query=subject_file,
                subject=f,
                evalue=evalue,
                out=recip_output_path_tab)
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

    # here we store the loci that fail filtering, so we can get the
    # subest of loci to include later
    bad_loci = []

    if reciprocal:
        logger.debug("shape of recip blast results")
        logger.debug(df2.shape)
        df2['genome'] = df2.subject_id.str.split('_').str.get(0)

    if not reciprocal:
        filtered = pd.DataFrame(columns=df1.columns)
        unq_subject = df1.subject_id.unique()
        unq_query = df1.query_id.unique()
        recip_hits = []
        nonrecip_hits = []  # should be renamed to nogood_hits
        for gene in unq_query:
            logger.debug("scoring %s" % (gene))
            tempdf1 = df1.loc[(df1["query_id"] == gene), ]
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
                bad_loci.append(gene)
            elif all(subset1["alignment_length"] > subset1["subject_length"] * min_length_percent):
                filtered = filtered.append(subset1)
            else:
                bad_loci.append(gene)
        return(filtered, bad_loci)


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
                    bad_loci.append(gene)
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
                        bad_loci.append(gene)

            # logger.debug(subset.shape)
    logger.debug("Non-reciprocal genes:")
    logger.debug(nonrecip_hits)
    logger.debug("Reciprocal genes:")
    logger.debug(recip_hits)
    logger.debug("filtered shape:")
    logger.debug(filtered.shape)
    return(filtered, bad_loci)


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
    """  merge the output from multiple BLAST runs of type 6 output (no headers)
    """
    # only grab .tab files, ie, the blast output
    with open(outfile_name, "a") as outf:
        for idx, f in enumerate(filelist):
            with open(f, "r") as inf:
                for line in inf:
                    outf.write(line)
    return outfile_name


def BLAST_tab_to_df(path):
    """ parse blast results when outputed as type 6 defaults with sequence length added
    """
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score", "subject_length"]
    with open(path) as ph:
        raw_csv_results = pd.read_csv(
            ph, comment="#", sep="\t", names=colnames)
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
                newrecs.append(newrec)
    with open(new_genbank, "w") as outf:
        for nr in newrecs:
            SeqIO.write(nr, outf, "genbank")


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
    if logger is None:
        logger = set_up_logging(outfile=os.path.join(output_root, "log.log"),
                            name="annofilt", verbosity=args.verbosity)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    date = str(datetime.datetime.now().strftime('%Y%m%d'))

    prokka_files = make_prokka_files_object(args.prokka_dir)

    # if args.full:
    #     all_loci = return_list_of_locus_tags(faa=prokka_files.faa)
    #     # check all the genes for completeness
    #     commands, paths_to_outputs, paths_to_recip_outputs = \
    #         make_prot_prot_blast_cmds(
    #             query_file=prokka_files.faa,
    #             evalue=args.min_evalue,
    #             reciprocal=args.reciprocal,
    #             subject_file=args.reference,
    #             threads=args.threads,
    #             output=output_root, date=date,
    #             logger=logger)
    # else:
    # check only selected the genes for completeness
    all_loci = return_list_of_locus_tags(gbk=prokka_files.gbk)
    genes_dirpath = os.path.join(output_root, "query_genes")
    gene_queries = []
    os.makedirs(genes_dirpath)
    # for each record, get the first and the last feature that isnt a "source"
    with open(prokka_files.gbk, "r") as inf:
        for rec in SeqIO.parse(inf, "genbank"):
            nfeats = len(rec.features)
            FIRST = True
            for idx, feat in enumerate(rec.features):
                # if the first or last records
                if feat.type == "source":
                    continue
                if (args.full or FIRST or idx == nfeats - 1):
                    FIRST = False
                    # print(type(feat.qualifiers.get('translation')))
                    gene = feat.extract(rec)
                    gene.seq = gene.seq.translate()
                    gene.id = feat.qualifiers.get("locus_tag")[0]
                    genepath = os.path.join(
                        genes_dirpath, gene.id + ".faa" )
                    SeqIO.write(gene, genepath, "fasta")
                    gene_queries.append(genepath)
    commands, paths_to_outputs, paths_to_recip_outputs = [], [], []
    for query in gene_queries:
        pcommands, ppaths_to_outputs, ppaths_to_recip_outputs = \
            make_prot_prot_blast_cmds(
                query_file=query,
                evalue=args.min_evalue,
                reciprocal=args.reciprocal,
                subject_file=args.reference,
                threads=args.threads,
                output=output_root, date=date,
                logger=logger)
        commands.extend(pcommands),
        paths_to_outputs.extend(ppaths_to_outputs)
        paths_to_recip_outputs.extend(ppaths_to_recip_outputs)
    logger.debug("cmds:" )
    logger.debug(commands)
    logger.debug("paths to outputs:" )
    logger.debug(paths_to_outputs )
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
    logger.debug(resultsdf)
    if args.reciprocal:
        logger.debug("paths to recip outputs:")
        logger.debug(paths_to_recip_outputs)
        post_recip_merged_tab = merge_outfiles(filelist=paths_to_recip_outputs,
                                               outfile_name=recip_merged_tab)
        recip_resultsdf = BLAST_tab_to_df(post_recip_merged_tab)
    else:
        recip_resultsdf = None


    filtered_hits, bad_loci = filter_BLAST_df(
        df1=resultsdf,
        df2=recip_resultsdf,
        reciprocal=args.reciprocal,
        min_percent=args.min_percent_id,
        min_length_percent=args.min_length,
        logger=logger)

    good_loci = [x for x in all_loci if x not in bad_loci]
    make_new_genbank(
        genbank=prokka_files.gbk,
        new_genbank=os.path.join(
            output_root,
            os.path.basename(os.path.splitext(prokka_files.faa)[0]) + "_filtered.gbk"),
        approved_accessions=good_loci,
        logger=logger)
    sys.stdout.write("Total\tKept\tLost\n")
    sys.stdout.write("{0}\t{1}\t{2}\n".format(
        len(all_loci), len(good_loci), len(bad_loci)))
    filtered_hits.to_csv(os.path.join(output_root, "filtered_hits.csv"))
    with open(os.path.join(output_root, "bad_loci.txt"), "w") as outf:
        for loci in bad_loci:
            outf.write(loci + "\n")
    with open(os.path.join(output_root, "good_loci.txt"), "w") as outf:
        for loci in good_loci:
            outf.write(loci + "\n")



if __name__ == "__main__":
    main()
