#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Nick Waters
"""

import os
import sys
import glob
import subprocess
import argparse
import logging
import shutil
import pandas as pd
import multiprocessing
import copy
from . import __version__
from . import shared_methods as sm

from argparse import Namespace
from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
# from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline


def get_args():  # pragma nocover
    parser = argparse.ArgumentParser(
        description="Blast assembly against core genome to find and " +
        "eliminate truncated genes due to bad assembly, " +
        "returning a assembly  with genes that meet a set of criteria. ",
        add_help=False)
    parser.add_argument(
        "reference",
        help="path to a nucleotide fasta" +
        "containing al the genes in the  genome; for instance, " +
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
        "--full",
        dest="full",
        action="store_true",
        help="check ALL genes for completeness, not just on ends of contigs." +
        " This is much slower.")
    optional.add_argument(
        "--local_quick",
        dest="local_quick",
        action="store_true",
        help="blast using prokkas' protein fasta;  can speed up preformance " +
        "by avoiding writing out genes to separate files, but this is " +
        "less efficient to paralellize")
    optional.add_argument(
        "-r",
        "--reciprocal",
        dest="reciprocal",
        action="store_true",
        help="reciprocal blast for stringent checking")
    optional.add_argument(
        "-p",
        "--min_id_percent",
        dest="min_id_percent",
        help="minimum percentage id (1-100)",
        type=float,
        default=90)
    optional.add_argument(
        "-l",
        "--min_length",
        dest="min_length_frac",
        help="minimum percentage length (0-1)",
        type=float,
        default=.9)
    optional.add_argument(
        "-e",
        "--min_evalue",
        dest="min_evalue",
        default=1,
        help="minimum e value")
    optional.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=4,
        help="number of threads to use; default to 4")
    optional.add_argument(
        "-b",
        "--blast_algorithm",
        dest="blast_algorithm",
        choices={"blastn", "blastp", "blastx", "tblastn", "tblastx"},
        default="blastn",
        help="Which BLAST algorithm to use. See BLAST documentation for " +
        "choosing. Most common options are blastn and tblastn, as Roary " +
        "provies a nucleotide pangenome.")
    optional.add_argument(
        "-s"
        "--sge",
        dest="sge",
        action="store_true",
        help="use Sun Grid Engine for job distribution; " +
        "default is to use multiprocessing")
    optional.add_argument(
        "--genes",
        dest="just_genes",
        action="store_true",
        help="only deal with the genes; default is to deal with " +
        "anything that has a locus tag")
    optional.add_argument(
        "--stringent",
        dest="stringent",
        action="store_true",
        help="reject genes with no blast hits")
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
    optional.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return args


def make_prokka_files_object(prokka_dir):
    """ givne a prokka output directory, find the AA fasta, genbank, and gff
    returns a Namespace object
    """
    checked = []
    for ext in ["faa", "gbk", "gff"]:
        files = glob.glob(prokka_dir + os.path.sep + "*" + ext)
        if len(files) == 0:
            raise ValueError("No %s file found in %s!" % (ext, prokka_dir))
        if len(files) > 1:
            raise ValueError(
                "Multiple %s files found in %s!" % (ext, prokka_dir))
        checked.append(files[0])
    # its prokka output; we can assume they all have the same prefix
    prefix = os.path.basename(os.path.splitext(checked[0])[0])
    return Namespace(faa=checked[0],
                     gbk=checked[1],
                     gff=checked[2],
                     prefix=prefix)


def return_list_of_locus_tags(gbk=None, faa=None, cds_only=False):
    lt_list = []
    nrecs = 0
    if gbk is not None:
        with open(gbk, "r") as inf:
            for rec in SeqIO.parse(inf, "genbank"):
                nrecs = nrecs + 1
                for feat in rec.features:
                    if feat.type != "CDS" and cds_only:
                        continue
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
    return (list(set(lt_list)), nrecs)


def make_blast_params(algo):
    # ensure figure out which type of fasta each seach scheme needs
    PROTEIN_SUBJECT = True
    PROTEIN_QUERY = True
    if algo in ["tblastn", "blastn", "tblastx"]:
        PROTEIN_SUBJECT = False
    if algo in ["tblastx", "blastx", "blastn"]:
        PROTEIN_QUERY = False
    return (PROTEIN_SUBJECT, PROTEIN_QUERY)


def make_blast_cmds(
        query_file, evalue, output, subject_file, threads=1, algo=None,
        protein_subject=False, makedb=False,
        reciprocal=False, logger=None):
    """given a file, make a blast cmd, and return path to output csv
    This should handle both protein and nucleotide references.
    """
    assert logger is not None, "must use logging"
    db_dir = os.path.join(
        output,
        os.path.splitext(os.path.basename(subject_file))[0])
    subjectdb = os.path.join(
        db_dir,
        os.path.splitext(os.path.basename(subject_file))[0])
    # first_record = SeqIO.parse(subject_file, "fasta").__next__()
    DBNAME = "protdb" if protein_subject else "nucdb"
    # logger.debug("Protein Subject: %s", PROTEIN_SUBJECT)
    # logger.debug("Protein Query: %s", PROTEIN_QUERY)
    if makedb:
        logger.debug("Creating  BLAST database")
        os.makedirs(db_dir, exist_ok=False)
        if protein_subject:
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
            qname + "_vs_" + DBNAME + ".tab")
        # if PROTEIN_QUERY and PROTEIN_SUBJECT:
        #     blast_cline = NcbiblastpCommandline(
        #         query=f,
        #         db=subjectdb, evalue=evalue, out=output_path_tab)
        # elif not PROTEIN_QUERY and not PROTEIN_SUBJECT and algo == "tblastx":
        #     blast_cline = NcbitblastxCommandline(
        #         query=f,
        #         db=subjectdb, evalue=evalue, out=output_path_tab)
        # elif not PROTEIN_QUERY and not PROTEIN_SUBJECT and algo == "blastn":
        #     blast_cline = NcbiblastnCommandline(
        #         query=f,
        #         db=subjectdb, evalue=evalue, out=output_path_tab)
        # elif not PROTEIN_QUERY and PROTEIN_SUBJECT:
        #     blast_cline = NcbiblastxCommandline(
        #         query=f,
        #         db=subjectdb, evalue=evalue, out=output_path_tab)
        # elif PROTEIN_QUERY and not PROTEIN_SUBJECT:
        #     blast_cline = NcbiblastxCommandline(
        #         query=f,
        #         db=subjectdb, evalue=evalue, out=output_path_tab)
        if algo == "blastn":
            blast_cline = NcbiblastnCommandline(
                query=f,
                db=subjectdb, evalue=evalue, out=output_path_tab)
        elif algo == "tblastn":
            blast_cline = NcbitblastnCommandline(
                query=f,
                db=subjectdb, evalue=evalue, out=output_path_tab)
        elif algo == "blastp":
            blast_cline = NcbiblastpCommandline(
                query=f,
                db=subjectdb, evalue=evalue, out=output_path_tab)
        elif algo == "blastx":
            blast_cline = NcbiblastxCommandline(
                query=f,
                db=subjectdb, evalue=evalue, out=output_path_tab)
        else:
            assert algo == "tblastx", "error parsing algrithm!"
            blast_cline = NcbitblastxCommandline(
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
            DBNAME + "_vs_" + qname + ".tab")
        if protein_subject:
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


def filter_BLAST_df(df1, df2, algo, min_evalue, min_length_frac,
                    min_id_percent, reciprocal, logger=None):
    """ results from pd.read_csv with default BLAST output 6 columns
    df1 must be genomes against genes, and df2 must be genes against genomes,
    because we have to split the names so all all the contigs are recognized
    as coming from one genome.  returns a df
    """
    assert logger is not None, "must use a logger"
    logger.debug("shape of blast results")
    df1['genome'] = df1.query_id.str.split('_').str.get(0)
    logger.debug(df1.shape)
    # here we
    if algo in ["tblastn", "tblastx", "blastp"]:
        CODON_ADJ = 3
    else:
        CODON_ADJ = 1
    # here we store the loci that fail filtering, so we can get the
    # subest of loci to exclude later
    bad_loci = []
    # subest of loci to lacking blast hits
    nohit_loci = []

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
            if tempdf1.empty:
                # no hit was found in the pangenome; lets be
                # conservative and ignore it
                nohit_loci.append(gene)
                continue
            subset1 = tempdf1.loc[
                (tempdf1["identity_perc"] > min_id_percent) &
                (tempdf1["bit_score"] == tempdf1["bit_score"].max()) &
                (tempdf1["evalue"] < min_evalue)]
            if subset1.empty:
                logger.debug("No full hits for %s", gene)
                logger.debug(tempdf1)
                bad_loci.append(gene)
            elif all(
                    (subset1["alignment_length"] * CODON_ADJ) >=
                    (subset1["subject_length"] * min_length_frac)):
                filtered = filtered.append(subset1)
            else:
                bad_loci.append(gene)
        return(filtered, bad_loci, nohit_loci)

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
                logger.debug("skipping %s in %s; no match found", gene, genome)
            else:
                subset1 = tempdf1.loc[
                    (tempdf1["identity_perc"] > min_id_percent) &
                    (tempdf1["bit_score"] == tempdf1["bit_score"].max())]
                # (tempdf1["alignement_l"] == tempdf1["bit_score"].max())]
                subset2 = tempdf2.loc[
                    (tempdf2["identity_perc"] > min_id_percent) &
                    (tempdf2["bit_score"] == tempdf2["bit_score"].max())]
                logger.debug("grouped df shape: ")
                logger.debug(tempdf1.shape)
                logger.debug("grouped df2 shape: ")
                logger.debug(tempdf2.shape)
                if subset1.empty or subset2.empty:
                    logger.debug("No reciprocol hits for %s in %s",
                                 gene, genome)
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
                    if subset1.iloc[0]["query_id"] == \
                       subset2.iloc[0]["subject_id"]:
                        recip_hits.append([gene, genome])
                        filtered = filtered.append(subset1)
                        logger.debug("Reciprocol hits for %s in %s!",
                                     gene, genome)
                    else:
                        nonrecip_hits.append([gene, genome])
                        logger.debug("No reciprocol hits for %s in %s",
                                     gene, genome)
                        bad_loci.append(gene)

            # logger.debug(subset.shape)
    logger.debug("Non-reciprocal genes:")
    logger.debug(nonrecip_hits)
    logger.debug("Reciprocal genes:")
    logger.debug(recip_hits)
    logger.debug("filtered shape:")
    logger.debug(filtered.shape)
    return(filtered, bad_loci)


def setup_blast_db(input_file, input_type="fasta", dbtype="prot",
                   title="blastdb", out="blastdb",
                   makeblastdb_exe='', logger=None):
    """
    This runs make blast db with the given parameters
    requires logging, os, subprocess, shutil
    """
    if makeblastdb_exe == '':
        makeblastdb_exe = shutil.which("makeblastdb")
        logger.debug("makeblastdb executable: %s", makeblastdb_exe)
    makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
                    "-out {4} 2> {4}.log").format(
                        makeblastdb_exe,
                        input_file,
                        input_type,
                        dbtype,
                        out)
    logger.debug("Making blast db: %s", makedbcmd)
    try:
        subprocess.run(makedbcmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        logging.debug("BLAST database '%s' created here: %s",
                      title, out)
        return 0
    except Exception as e:
        print(e)
        logging.error("Something bad happened when trying to make " +
                      "a blast database")
        sys.exit(1)


def merge_outfiles(filelist, outfile_name):
    """  merge the output from multiple BLAST runs of type 6 output
    (no headers)
    """
    # only grab .tab files, ie, the blast output
    with open(outfile_name, "a") as outf:
        for idx, f in enumerate(filelist):
            with open(f, "r") as inf:
                for line in inf:
                    outf.write(line)
    return outfile_name


def BLAST_tab_to_df(path):
    """ parse blast results when outputed as type 6
    defaults with sequence length added
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
                try:
                    # if feat.type in ["source", "misc_feature",
                    #                  "repeat_region", "assembly_gap"]:
                    #     continue
                    if not feat.type:
                        continue
                    if feat.type != "CDS":
                        continue
                    if feat.qualifiers.get("locus_tag")[0] in \
                       approved_accessions:
                        newrec.features.append(feat)
                    else:
                        logger.debug(feat.qualifiers.get("locus_tag")[0] +
                                     " was rejected")
                except TypeError as e:
                    logger.error(feat)
                    raise(e)
            newrecs.append(newrec)
    with open(new_genbank, "w") as outf:
        for nr in newrecs:
            SeqIO.write(nr, outf, "genbank")


def make_filter_gff_cmd(gff, baddies, newgff):
    """ given a gff file and a file of unwanted locus tags, run inverse grep

    Note 2019-04-25 this is a ticking time bomb
    """
    # -f means get pattern from file
    # -v means return inverse match
    return "grep {0} -f {1} -v > {2}".format(gff, baddies, newgff)


def get_genewise_blast_cmds(output_root, prokka_files, args,
                            debug=False, logger=None):
    genes_dirpath = os.path.join(output_root, "query_genes")
    gene_queries = []
    os.makedirs(genes_dirpath)
    TRANSLATE = args.blast_algorithm in ["tblastn", "blastp"]
    ext = ".faa" if TRANSLATE else ".fna"
    # for each record, get the first and the last feature that isnt a "source"
    logger.debug("writing out genes for blasting")
    with open(prokka_files.gbk, "r") as inf:
        for rec in SeqIO.parse(inf, "genbank"):
            nfeats = len(rec.features)
            logger.debug("found %d features on %s", nfeats, rec.id)
            FIRST = True
            for idx, feat in enumerate(rec.features):
                # note we dont disciminate against misc_feature because some
                # can technically have locus tags, and are therefor fair
                # game for downstreaem analysis
                # If a weird feature comes up that doesnt have a locus
                # tag, it is logged and ignored.
                if feat.type in ["source", "repeat_region",
                                 "assembly_gap"]:
                    continue
                if (args.full or FIRST or idx == nfeats - 1):
                    FIRST = False
                    gene = feat.extract(rec)
                    if TRANSLATE:
                        try:
                            gene.seq = gene.seq.translate()
                        except BiopythonWarning as bpw:
                            logger.warning(bpw)
                    try:
                        gene.id = feat.qualifiers.get("locus_tag")[0]
                    except TypeError:
                        logger.warning("Could not get a locus tag for the " +
                                       "following feature: %s", feat)
                        continue
                    genepath = os.path.join(
                        genes_dirpath, gene.id + ext)
                    SeqIO.write(gene, genepath, "fasta")
                    gene_queries.append(genepath)
    commands, paths_to_outputs, paths_to_recip_outputs = [], [], []
    logger.debug("making blast commands")
    blast_params = make_blast_params(algo=args.blast_algorithm)
    for idx, query in enumerate(gene_queries):
        pcommands, ppaths_to_outputs, ppaths_to_recip_outputs = \
            make_blast_cmds(
                query_file=query,
                evalue=1,
                reciprocal=args.reciprocal,
                subject_file=args.reference,
                protein_subject=blast_params[0],
                threads=1,
                makedb=idx == 0 and not debug,
                algo=args.blast_algorithm,
                output=output_root,
                logger=logger)
        commands.extend(pcommands),
        paths_to_outputs.extend(ppaths_to_outputs)
        paths_to_recip_outputs.extend(ppaths_to_recip_outputs)
    return (commands, paths_to_outputs, paths_to_recip_outputs)


def get_genome_blast_cmds(output_root, prokka_files, args, logger=None):
    commands, paths_to_outputs, paths_to_recip_outputs = \
        make_blast_cmds(
            query_file=prokka_files.faa,
            evalue=1,
            reciprocal=args.reciprocal,
            subject_file=args.reference,
            protein_subject=True,
            threads=args.threads,
            output=output_root,
            makedb=True,
            algo=args.blast_algorithm,
            logger=logger)
    return (commands, paths_to_outputs, paths_to_recip_outputs)


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
        logger = sm.set_up_logging(
            outfile=os.path.join(output_root, "annofilt.log"),
            name="annofilt", verbosity=args.verbosity)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    logger.debug("processing prokka output")
    prokka_files = make_prokka_files_object(args.prokka_dir)
    # check only selected the genes for completeness
    logger.info("Making list of all locus tags in assembly")
    all_loci, nrecs = return_list_of_locus_tags(gbk=prokka_files.gbk)
    logger.info("Generating BLAST commands")
    if not args.local_quick and args.full:
        # build blast command using prokka faa as query
        commands, paths_to_outputs, paths_to_recip_outputs = \
            get_genewise_blast_cmds(
                output_root=output_root,
                prokka_files=prokka_files,
                args=args, logger=logger)
    else:
        # wtrite out genes to use as queries (
        # if full, all the genes, else just the 1rst and last)
        commands, paths_to_outputs, paths_to_recip_outputs = \
            get_genewise_blast_cmds(
                output_root=os.path.join(output_root, "BLAST_results"),
                prokka_files=prokka_files,
                args=args, logger=logger)

    logger.debug("writing out blast commands for reference")
    with open(os.path.join(output_root, "blast_cmds.txt"), "w") as outf:
        for c in commands:
            outf.write(c + "\n")
    logger.debug("paths to outputs:")
    logger.debug(paths_to_outputs)
    logger.info("Running BLAST commands")
    if not args.sge:
        pool = multiprocessing.Pool(processes=args.threads)
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
        raise ValueError("SGE scheduling no yet implemented!")
    logger.info("Consolidating BLAST outputs")
    merged_tab = merge_outfiles(filelist=paths_to_outputs,
                                outfile_name=os.path.join(
                                    output_root, "merged_results.tab"))
    resultsdf = BLAST_tab_to_df(merged_tab)

    if args.reciprocal:
        recip_merged_tab = merge_outfiles(
            filelist=paths_to_recip_outputs,
            outfile_name=os.path.join(output_root, "recip_merged_results.tab"))
        recip_resultsdf = BLAST_tab_to_df(recip_merged_tab)
    else:
        recip_resultsdf = None

    filtered_hits, bad_loci, nohit_loci = filter_BLAST_df(
        df1=resultsdf,
        df2=recip_resultsdf,
        algo=args.blast_algorithm,
        reciprocal=args.reciprocal,
        min_evalue=args.min_evalue,
        min_id_percent=args.min_id_percent,
        min_length_frac=args.min_length_frac,
        logger=logger)
    if args.stringent:
        good_loci = [x for x in all_loci if
                     x not in bad_loci.extend(nohit_loci)]
    else:
        good_loci = [x for x in all_loci if x not in bad_loci]
    # write out a .gbk file that lacks the genes deemed "bad" (truncated, etc)
    logger.info("Writing out filtered Genbank file")
    make_new_genbank(
        genbank=prokka_files.gbk,
        new_genbank=os.path.join(
            output_root,
            prokka_files.prefix + "_filtered.gbk"),
        approved_accessions=good_loci,
        logger=logger)
    filtered_hits.to_csv(os.path.join(output_root, "filtered_hits.csv"))

    # write out locus tags that will be used with grep to filter the gff file
    baddies_file = os.path.join(output_root, "bad_loci.txt")
    with open(baddies_file, "w") as outf:
        for loci in bad_loci:
            outf.write(loci + "\n")
    with open(os.path.join(output_root, "good_loci.txt"), "w") as outf:
        for loci in good_loci:
            outf.write(loci + "\n")
    with open(os.path.join(output_root, "nohit_loci.txt"), "w") as outf:
        for loci in nohit_loci:
            outf.write(loci + "\n")

    newgff = os.path.join(
        output_root,
        prokka_files.prefix + "_filtered.gff")
    logger.debug("writing filtered gff to %s", newgff)
    if len(bad_loci) == 0:
        filter_cmd = "cp {0} {1}".format(prokka_files.gff, newgff)
    else:
        filter_cmd = make_filter_gff_cmd(
            gff=prokka_files.gff,
            baddies=baddies_file,
            newgff=newgff)
    logger.info("Creating filtered GFF files")
    subprocess.run(filter_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, check=True)
    if args.stringent:
        logger.info("%i genes lacked any BLAST hits, and will be rejected",
                    len(nohit_loci))
        logger.info("Total from %i contigs: %i\tKept: %i\tLost: %i",
                    nrecs, len(all_loci), len(good_loci),
                    len(bad_loci) + len(nohit_loci))
        sys.stdout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
            nrecs, prokka_files.prefix, len(all_loci), len(good_loci),
            len(bad_loci) + len(nohit_loci)))
    else:
        logger.info("%s genes lacked any BLAST hits, and were retained",
                    len(nohit_loci))
        logger.info("Total from %i contigs: %i\tKept: %i\tLost: %i",
                    nrecs, len(all_loci), len(good_loci), len(bad_loci))
        sys.stdout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
            nrecs, prokka_files.prefix, len(all_loci),
            len(good_loci), len(bad_loci)))

    logger.debug("Done!")


if __name__ == "__main__":
    main()
