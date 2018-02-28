import os
import unittest
import sys
import time
import copy
from Bio import SeqIO

from annofilt import make_new_genbank

@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class annofilt(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "output_annofilt_tests")
        self.ref_gb = os.path.join(os.path.dirname(__file__),
                                   "test_data", "NC_011751.1.gb")
        self.to_be_removed = []
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir, exist_ok=True)
        self.startTime = time.time()

    def test_loop_through_genbank(self):
        """ construct spades command that check for file presense
        this is useful when multiprocessing and unable to check before
        sending the command out
        """
        make_new_genbank(genbank=self.ref_gb,
                         new_genbank=os.path.join(self.test_dir, "new.gb")
                         approved_accessions=["ECUMN_0005"])

    def tearDown(self):
        """ delete temp files if no errors
        """
        for filename in self.to_be_removed:
            os.unlink(filename)
        pass

if __name__ == '__main__':
    unittest.main()
