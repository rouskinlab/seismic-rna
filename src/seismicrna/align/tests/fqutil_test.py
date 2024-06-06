"""

Tests for Alignment FASTQ Utilities Module

========================================================================

"""

import unittest as ut

from ...core.path import INP_TEST_DIR


class TestRunFASTQC(ut.TestCase):
    """ Test function `run_fastqc`. """

    test_fq = "test-sample.fq"
    test_fq1 = "test-sample_R1.fq"
    test_fq2 = "test-sample_R2.fq"

    def setUp(self):
        """ Make temporary FASTQ file(s) on which to run FASTQC. """
        INP_TEST_DIR.mkdir(parents=True, exist_ok=False)

    def tearDown(self):
        """ Delete the test files and directories. """

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
