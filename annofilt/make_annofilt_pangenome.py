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
import shutil


from . import shared_methods as sm


def get_args():  # pragma nocover
    parser = argparse.ArgumentParser(
        description="Given a directory of genomes, create a " +
        "pangenome using prokka and Roary from complete genomes " +
        "of the organism",
        add_help=False)
    parser.add_argument(
        "-g", "--genomes",
        help="genomes")
    parser.add_argument(
        "-o",
        "--output",
        help="output dir", required=True)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        "--add_roary_args",
        default="-e -n",
        help="Args to pass on to Roary for pangenome construction. " +
        "Default args to give to roary. " +
        "-e -n means a fast core multifast is generated with MAFFT")
    optional.add_argument(
        "-t",
        "--threads",
        dest='threads', action="store",
        default=2, type=int,
        help="number of cores to use " +
        "default: %(default)s")
    optional.add_argument(
        "-e",
        "--experiment_name",
        dest='experiment_name', action="store",
        default="annofilt",
        help="name default: %(default)s")
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
    for exe in ["prokka", "roary"]:
        if shutil.which(exe) is None:
            raise ValueError("%s executable not found" % exe)


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
    check_exes()
    if logger is None:
        logger = sm.set_up_logging(
            outfile=os.path.join(output_root, "make_annofilt_pangenome.log"),
            name="annofilt",
            verbosity=args.verbosity)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    genomes = glob.glob(args.genomes + "*.fna")
    if len(genomes) < 2:
        raise ValueError("Prokka needs a minimum of 2 genomes to run! " +
                         "check the contents of your --genomes dir.  " +
                         "Genomes need to end in .fna")
    logger.info("Running Prokka commands")
    for i, genome in enumerate(genomes):
        thisname = os.path.basename(os.path.splitext(genome)[0])
        outdir = os.path.join(output_root,
                              os.path.basename(os.path.splitext(genome)[0]))
        if not os.path.exists(outdir):
            prokka_cmd = str(
                "prokka --outdir {outdir} --prefix " +
                "{args.experiment_name}-{i} --compliant --genus Genus " +
                "--species species --cpus {args.threads} " +
                "{genome}").format(**locals())
            logger.debug(prokka_cmd)
            subprocess.run(prokka_cmd, shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, check=True)
    roary_out = os.path.join(output_root, args.experiment_name)
    roary_cmd = str(
        "roary -p {args.threads} -f {roary_out} {args.add_roary_args} " +
        "-v {output_root}/*/*.gff >> {output_root}/roary.log 2>&1").format(
            **locals())
    logger.info("Preparing Roary cmd")
    logger.debug(roary_cmd)
    subprocess.run(roary_cmd, shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, check=True)
    logger.debug("Done!")


if __name__ == "__main__":
    main()
