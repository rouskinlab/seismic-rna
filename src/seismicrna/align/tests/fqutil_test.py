"""

Tests for Alignment FASTQ Utilities Module

========================================================================

"""

import unittest as ut

from ...core.path import IN_TEST_DIR


class TestRunFASTQC(ut.TestCase):
    """ Test function `run_fastqc`. """

    temp_fq = "test-sample.fq"
    temp_fq1 = "test-sample_R1.fq"
    temp_fq2 = "test-sample_R2.fq"

    def setUp(self):
        """ Make temporary FASTQ file(s) on which to run FASTQC. """
        IN_TEST_DIR.mkdir(parents=True, exist_ok=False)

    def tearDown(self):
        """ Delete the test files and directories. """
