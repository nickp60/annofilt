import os
import unittest
import glob
import sys
import time
import copy
import logging
import glob
import shutil
from argparse import Namespace
from Bio import SeqIO

from annofilt import annofilt as af

logger = logging

@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class annofilt(unittest.TestCase):
    """ tests for annofilt.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_annofilt_tests")
        self.data_dir = os.path.join(os.path.dirname(__file__),
                                     "testdata")
        self.merged_tab = os.path.join(self.data_dir,
                                       "merged_results.tab")
        self.tofilter_tab = os.path.join(self.data_dir,
                                       "tofilter.tab")
        self.ref_pangenome = os.path.join(os.path.dirname(__file__),
                                          "testdata", "miniref.fna")
        self.ref_gb = os.path.join(os.path.dirname(__file__),
                                   "testdata", "PROKKA.gbk")
        self.ref_faa = os.path.join(os.path.dirname(__file__),
                                    "testdata", "PROKKA.faa")
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.startTime = time.time()

    def test_loop_through_genbank(self):
        """ test creation of new gb after filtering.
        Not the greatest test
        """
        newf = os.path.join(self.test_dir, "new.gbk")
        af.make_new_genbank(
            genbank=self.ref_gb,
            new_genbank=newf,
            approved_accessions=["IPDMBIDP_00043"], logger=logger)
        nlines = 0
        with open(newf, "r") as inf:
            for line in inf:
                nlines = nlines + 1
        self.assertEqual(nlines, 1516)

    def test_BLAST_tab_to_df(self):
        self.assertEqual(
            af.BLAST_tab_to_df(self.merged_tab).shape,
            (26, 13))

    def test_merge_outfiles(self):
        blastfiles = glob.glob(
            os.path.join(self.data_dir, "blasthits", "") + "*.tab")
        dest = os.path.join(self.test_dir, "merged.tab")
        af.merge_outfiles(filelist=blastfiles,
                          outfile_name=dest)
        self.assertEqual(
            af.BLAST_tab_to_df(self.merged_tab).shape,
            af.BLAST_tab_to_df(dest).shape
            )
        self.to_be_removed.append(dest)

    def test_make_prokka_files_object(self):
        file_ob = af.make_prokka_files_object(self.data_dir)
        self.assertEqual(file_ob.faa, self.ref_faa)

    def test_make_prokka_files_object_bad(self):
        with self.assertRaises(ValueError):
            badfile_ob = af.make_prokka_files_object(self.test_dir)

    def test_return_list_of_locus_tags_gbk(self):
        self.assertEqual(
            len(af.return_list_of_locus_tags(gbk=self.ref_gb)[0]),
            89)
        self.assertEqual(
            len(af.return_list_of_locus_tags(gbk=self.ref_gb, cds_only=True)[0]),
            65)

    def test_return_list_of_locus_tags_faa(self):
        self.assertEqual(
            len(af.return_list_of_locus_tags(faa=self.ref_faa)[0]),
            65)

    def test_return_list_of_locus_tags_both(self):
        """ we know this fails;
        the faa file only has proteins, while the gbk file has all annotations
        """
        gbkl, nrecsg = af.return_list_of_locus_tags(gbk=self.ref_gb, cds_only=True)
        faal, nrecsf = af.return_list_of_locus_tags(faa=self.ref_faa)
        for i in range(len(gbkl)):
            self.assertEqual(gbkl[i], faal[i])

    def test_make_filter_gff_cmd(self):
        self.assertEqual(
            af.make_filter_gff_cmd(gff="gff", baddies="yuk", newgff="ohboy"),
            "grep gff -f yuk -v > ohboy"
        )

    def test_blast_cmds(self):
        self.maxDiff =None
        cmds, opaths, ropaths = af.make_blast_cmds(
            query_file = "query.fa", subject_file = self.ref_faa, makedb=False,
            evalue=1, output=self.test_dir, threads=3, algo="tblastn",
            reciprocal=False, protein_subject=False, logger=logger)
        self.assertEqual(
            cmds[0],
            "tblastn -out " +
            self.test_dir +
            "/query_vs_nucdb.tab -query query.fa -db " +
            self.test_dir +
            "/PROKKA/PROKKA -evalue 1 -num_threads 3 -num_alignments " +
            "20 -outfmt '6 qaccver saccver pident length mismatch " +
            "gapopen qstart qend sstart send evalue bitscore slen'")
        self.to_be_removed.extend(glob.glob(os.path.join(
            self.test_dir, "PROKKA", "*")))
        # self.to_be_removed.append(
        #     os.path.join(self.test_dir, "PROKKA"))
        # self.to_be_removed.append(os.path.join(self.test_dir, "output"))

    # @unittest.skipIf(shutil.which("makeblastdb") is None,
    #                  "blast executables not found; skipping test")
    # if we point this to an existing dir, no need for make blast db
    def test_get_genewise_blast_cmds(self):
        file_ob = af.make_prokka_files_object(self.data_dir)
        args=Namespace(min_evalue=1, reciprocal=False,
                       blast_algorithm="tblastn",
                       threads=1, reference=self.ref_pangenome, full=False)
        commands, paths_to_outputs, paths_to_recip_outputs = \
            af.get_genewise_blast_cmds(output_root=self.test_dir, prokka_files=file_ob,
                                       args=args, debug=True, logger=logger)
        self.to_be_removed.extend(glob.glob(os.path.join(
            self.test_dir, "query_genes", "*")))
        self.to_be_removed.append(
            os.path.join(self.test_dir, "query_genes"))

    def test_filter_BLAST_df_nofilt(self):
        nofilt = af.filter_BLAST_df(
            df1=af.BLAST_tab_to_df(self.tofilter_tab),
            df2=None,
            algo="blastn",
            min_length_frac=.1, min_id_percent=0, min_evalue=1, reciprocal=False, logger=logger)
        self.assertEqual(len(nofilt[1]), 0)

    def test_filter_BLAST_df_low_id(self):
        filt_id = af.filter_BLAST_df(
            df1=af.BLAST_tab_to_df(self.tofilter_tab),
            df2=None,
            algo="blastn",
            min_length_frac=.1, min_id_percent=100, min_evalue=1, reciprocal=False, logger=logger)
        self.assertEqual(len(filt_id[1]), 1)

    def test_filter_BLAST_df_low_length(self):
        filt_length = af.filter_BLAST_df(
            df1=af.BLAST_tab_to_df(self.tofilter_tab),
            df2=None,
            algo="blastn",
            min_length_frac=.99, min_id_percent=0, min_evalue=1, reciprocal=False, logger=logger)
        self.assertEqual(len(filt_length[1]), 1)

    def test_filter_BLAST_df_low_evalue(self):
        filt_evalue = af.filter_BLAST_df(
            df1=af.BLAST_tab_to_df(self.tofilter_tab),
            df2=None,
            algo="blastn",
            min_length_frac=.1, min_id_percent=0, min_evalue=.0001, reciprocal=False, logger=logger)
        self.assertEqual(filt_evalue[0].shape[0], 2)
        self.assertEqual(len(filt_evalue[1]), 0)

    def test_make_blast_params_blastn(self):
        self.assertEqual(af.make_blast_params("blastn"),
                         #subject, query
                         (False, False))

    def test_make_blast_params_tblastn(self):
        self.assertEqual(af.make_blast_params("tblastn"),
                         #subject, query
                         (False, True))

    def test_make_blast_params_blastp(self):
        self.assertEqual(af.make_blast_params("blastp"),
                         #subject, query
                         (True, True))

    def test_make_blast_params_blastx(self):
        self.assertEqual(af.make_blast_params("blastx"),
                         #subject, query
                         (True, False))

    def test_make_blast_params_tblastx(self):
        self.assertEqual(af.make_blast_params("blastn"),
                         #subject, query
                         (False, False))

    def tearDown(self):
        """ delete temp files if no errors
        """
        for f in self.to_be_removed:
            if os.path.isfile(f):
                os.unlink(f)
            else:
                os.removedirs(f)


if __name__ == '__main__':
    unittest.main()
