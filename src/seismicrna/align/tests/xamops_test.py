import unittest as ut
from pathlib import Path

from seismicrna.align.xamops import flags_cmds


class TestFlagsCmds(ut.TestCase):
    XAM_INP = Path("a/b/c.bam")
    XAM_OUT = Path("x/y/z.bam")
    TMP_PFX = XAM_OUT.with_suffix("")

    def test_0_flags_bam(self):
        result = flags_cmds(self.XAM_INP, self.XAM_OUT)
        expect = [f"samtools view -@ 0 -o {self.XAM_OUT} {self.XAM_INP}"]
        self.assertEqual(result, expect)

    def test_0_flags_none(self):
        result = flags_cmds(self.XAM_INP, None)
        expect = [f"samtools view -@ 0 -u {self.XAM_INP}"]
        self.assertEqual(result, expect)

    def test_1_flag_bam(self):
        result = flags_cmds(self.XAM_INP,
                            self.XAM_OUT,
                            flags_req=16,
                            flags_exc=128)
        expect = [
            f"samtools view -@ 0 -f 16 -F 128 -o {self.XAM_OUT} {self.XAM_INP}"
        ]
        self.assertEqual(result, expect)

    def test_1_flag_none(self):
        result = flags_cmds(self.XAM_INP,
                            None,
                            flags_req=16,
                            flags_exc=128)
        expect = [f"samtools view -@ 0 -f 16 -F 128 -u {self.XAM_INP}"]
        self.assertEqual(result, expect)

    def test_2_flags_bam(self):
        result = flags_cmds(self.XAM_INP,
                            self.XAM_OUT,
                            flags_req=[16, 32],
                            flags_exc=[128, 256])
        expect = [
            f"( samtools view -@ 0 -f 16 -F 128 -h {self.XAM_INP} "
            f"; samtools view -@ 0 -f 32 -F 256 {self.XAM_INP} )",
            f"samtools collate -@ 0 -f -T {self.TMP_PFX} -o {self.XAM_OUT} -"
        ]
        self.assertEqual(result, expect)

    def test_2_flags_none(self):
        result = flags_cmds(self.XAM_INP,
                            None,
                            flags_req=[16, 32],
                            flags_exc=[128, 256])
        expect = [
            f"( samtools view -@ 0 -f 16 -F 128 -h {self.XAM_INP} "
            f"; samtools view -@ 0 -f 32 -F 256 {self.XAM_INP} )",
            "samtools collate -@ 0 -f -O -u -"
        ]
        self.assertEqual(result, expect)


if __name__ == "__main__":
    ut.main(verbosity=2)
