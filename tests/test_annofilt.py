import os
import unittest
import glob
import sys
import time
import copy
import logging
import glob
from Bio import SeqIO

from annofilt import annofilt as af

logger = logging
@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class annofilt(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_annofilt_tests")
        self.data_dir = os.path.join(os.path.dirname(__file__),
                                     "testdata")
        self.merged_tab = os.path.join(self.data_dir,
                                       "merged_results.tab")
        self.ref_gb = os.path.join(os.path.dirname(__file__),
                                   "testdata", "PROKKA.gbk")
        self.ref_faa = os.path.join(os.path.dirname(__file__),
                                    "testdata", "PROKKA.faa")
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.startTime = time.time()

    def tet_loop_through_genbank(self):
        """ construct spades command that check for file presense
        this is useful when multiprocessing and unable to check before
        sending the command out
        """
        af.make_new_genbank(genbank=self.ref_gb,
                         new_genbank=os.path.join(self.test_dir, "new.gb"),
                         approved_accessions=["ECUMN_0005"], logger=logger)

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

    def test_filter_BLAST_df(self):
        filtered_hits, bad_loci = af.filter_BLAST_df(
            df1=af.BLAST_tab_to_df(self.merged_tab),
            df2="notafile",
            reciprocal=False,
            min_percent=.75,
            min_length_percent=.75,
            logger=logger)

    def test_make_prokka_files_object(self):
        file_ob = af.make_prokka_files_object(self.data_dir)
        self.assertEqual(file_ob.faa, self.ref_faa)

    def test_make_prokka_files_object_bad(self):
        with self.assertRaises(ValueError):
            badfile_ob = af.make_prokka_files_object(self.test_dir)

    def test_return_list_of_locus_tags_gbk(self):
        self.assertEqual(
            len(af.return_list_of_locus_tags(gbk=self.ref_gb)),
            89)
        self.assertEqual(
            len(af.return_list_of_locus_tags(gbk=self.ref_gb, cds_only=True)),
            65)

    def test_return_list_of_locus_tags_faa(self):
        self.assertEqual(
            len(af.return_list_of_locus_tags(faa=self.ref_faa)),
            65)

    def test_return_list_of_locus_tags_both(self):
        """ we know this fails;
        the faa file only has proteins, while the gbk file has all annotations
        """
        gbkl = af.return_list_of_locus_tags(gbk=self.ref_gb, cds_only=True)
        faal = af.return_list_of_locus_tags(faa=self.ref_faa)
        for i in range(len(gbkl)):
            self.assertEqual(gbkl[i], faal[i])

    def test_make_filter_gff_cmd(self):
        self.assertEqual(
            af.make_filter_gff_cmd(gff="gff", baddies="yuk", newgff="ohboy"),
            "grep gff -f yuk -v > ohboy"
        )

    def test_blast_cmds(self):
        cmds, opaths, ropaths = af.make_prot_prot_blast_cmds(
            query_file = "query.fa", subject_file = self.ref_faa,
            evalue=1, output=self.test_dir + "/output/", threads=3,
            reciprocal=False, protein_subject=True, logger=logger)
        self.assertEqual(
            cmds[0],
            "blastp -out " +
            self.test_dir +
            "/output/query_vs_protdb.tab -query query.fa -db " +
            self.test_dir +
            "/output/PROKKA/PROKKA -evalue 1 -num_threads 3 -num_alignments " +
            "20 -outfmt '6 qaccver saccver pident length mismatch " +
            "gapopen qstart qend sstart send evalue bitscore slen'")
        self.to_be_removed.extend(glob.glob(os.path.join(
            self.test_dir, "output", "PROKKA", "*")))
        self.to_be_removed.append(
            os.path.join(self.test_dir, "output", "PROKKA"))
        # self.to_be_removed.append(os.path.join(self.test_dir, "output"))

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
