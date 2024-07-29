import unittest as ut
from itertools import chain, product, combinations
from string import ascii_letters

from seismicrna.core.rel.pattern import HalfRelPattern, RelPattern


def powerset(items: list | set):
    """ Set of all subsets of items, including empty and all items. """
    return chain.from_iterable(combinations(items, n)
                               for n in range(len(items) + 1))


def iter_codes(all_letters: bool = False):
    if all_letters:
        ab = product(ascii_letters, repeat=2)
    else:
        ab = product("acgtACGT", "acgtdiACGTDI")
    for a, b in ab:
        plain = HalfRelPattern.fmt_plain.format(a, b)
        fancy = HalfRelPattern.fmt_fancy.format(a, b)
        yield a, b, plain, fancy


class TestHalfRelPattern(ut.TestCase):

    def test_bases_muts(self):
        self.assertEqual(HalfRelPattern.ref_bases, "ACGT")
        self.assertEqual(HalfRelPattern.read_bases, "ACGTDI")
        self.assertEqual(HalfRelPattern.mut_bits,
                         bytes([16, 32, 64, 128, 2, 12]))

    def test_re_patterns(self):
        self.assertEqual(HalfRelPattern.fmt_plain, "{}{}")
        self.assertEqual(HalfRelPattern.fmt_fancy, "{} -> {}")
        for a, b, plain, fancy in iter_codes(all_letters=True):
            with self.subTest(a=a, b=b):
                if a in "acgt" and b in "acgtdi":
                    self.assertRegex(plain, HalfRelPattern.ptrn_plain)
                else:
                    self.assertNotRegex(plain, HalfRelPattern.ptrn_plain)
                if a in "ACGT" and b in "ACGTDI":
                    self.assertRegex(fancy, HalfRelPattern.ptrn_fancy)
                else:
                    self.assertNotRegex(fancy, HalfRelPattern.ptrn_fancy)

    def test_as_match(self):
        for a, b, plain, fancy in iter_codes(all_letters=True):
            with self.subTest(a=a, b=b):
                if a in "acgtACGT" and b in "acgtdiACGTDI":
                    self.assertEqual(HalfRelPattern.as_match(plain).group(),
                                     plain.lower())
                    self.assertEqual(HalfRelPattern.as_match(fancy).group(),
                                     fancy.upper())
                else:
                    self.assertRaisesRegex(ValueError,
                                           "Failed to match code",
                                           HalfRelPattern.as_match,
                                           plain)
                    self.assertRaisesRegex(ValueError,
                                           "Failed to match code",
                                           HalfRelPattern.as_match,
                                           fancy)

    def test_as_plain(self):
        for _, _, plain, fancy in iter_codes():
            with self.subTest(fancy=fancy):
                self.assertEqual(HalfRelPattern.as_plain(plain), plain.lower())
                self.assertEqual(HalfRelPattern.as_plain(fancy), plain.lower())
                self.assertEqual(HalfRelPattern.as_fancy(plain), fancy.upper())
                self.assertEqual(HalfRelPattern.as_fancy(fancy), fancy.upper())

    def test_compile_example(self):
        codes = ["A -> C", "A -> T", "C -> D", "G -> A",
                 "G -> G", "G -> I", "T -> A", "T -> T"]
        expect = {"A": 160, "C": 2, "G": 29, "T": 17}
        self.assertDictEqual(HalfRelPattern.compile(codes), expect)

    def test_decompile_example(self):
        patterns = {"A": 160, "C": 2, "G": 29, "T": 17}
        expect = ["A -> C", "A -> T", "C -> D", "G -> A",
                  "G -> G", "G -> I", "T -> A", "T -> T"]
        self.assertListEqual(sorted(HalfRelPattern.decompile(patterns)), expect)

    def test_compile_decompile(self):
        for _, _, plain, fancy in iter_codes():
            with self.subTest(fancy=fancy):
                compiled = HalfRelPattern.compile([fancy])
                self.assertDictEqual(compiled, HalfRelPattern.compile([plain]))
                self.assertListEqual(list(HalfRelPattern.decompile(compiled)),
                                     [fancy.upper()])

    def test_from_report_format(self):
        for _, _, _, fancy in iter_codes():
            with self.subTest(fancy=fancy):
                pattern = HalfRelPattern.from_report_format([fancy])
                self.assertIsInstance(pattern, HalfRelPattern)
                self.assertListEqual(pattern.codes, [fancy.upper()])
                expect = HalfRelPattern.compile([fancy])
                self.assertEqual(pattern.a, expect["A"])
                self.assertEqual(pattern.c, expect["C"])
                self.assertEqual(pattern.g, expect["G"])
                self.assertEqual(pattern.t, expect["T"])

    def test_from_counts(self):
        pattern = HalfRelPattern.from_counts()
        self.assertIsInstance(pattern, HalfRelPattern)
        self.assertDictEqual(pattern.patterns,
                             {"A": 0, "C": 0, "G": 0, "T": 0})
        pattern = HalfRelPattern.from_counts(count_ref=True)
        self.assertDictEqual(pattern.patterns,
                             {"A": 1, "C": 1, "G": 1, "T": 1})
        pattern = HalfRelPattern.from_counts(count_sub=True)
        self.assertDictEqual(pattern.patterns,
                             {"A": 224, "C": 208, "G": 176, "T": 112})
        pattern = HalfRelPattern.from_counts(count_del=True)
        self.assertDictEqual(pattern.patterns,
                             {"A": 2, "C": 2, "G": 2, "T": 2})
        pattern = HalfRelPattern.from_counts(count_ins=True)
        self.assertDictEqual(pattern.patterns,
                             {"A": 12, "C": 12, "G": 12, "T": 12})
        pattern = HalfRelPattern.from_counts(count_ref=True,
                                             discount=["A -> G", "cc", "ta"])
        self.assertDictEqual(pattern.patterns,
                             {"A": 1, "C": 0, "G": 1, "T": 1})
        pattern = HalfRelPattern.from_counts(count_sub=True,
                                             discount=["A -> G", "cc", "ta"])
        self.assertDictEqual(pattern.patterns,
                             {"A": 160, "C": 208, "G": 176, "T": 96})
        pattern = HalfRelPattern.from_counts(count_ref=True,
                                             count_sub=True,
                                             count_del=True,
                                             count_ins=True,
                                             discount=["ag", "G -> C"])
        self.assertDictEqual(pattern.patterns,
                             {"A": 175, "C": 223, "G": 159, "T": 127})

    def test_allc(self):
        pattern = HalfRelPattern.allc()
        self.assertIsInstance(pattern, HalfRelPattern)
        self.assertDictEqual(pattern.patterns,
                             {"A": 239, "C": 223, "G": 191, "T": 127})

    def test_muts(self):
        pattern = HalfRelPattern.muts()
        self.assertIsInstance(pattern, HalfRelPattern)
        self.assertDictEqual(pattern.patterns,
                             {"A": 238, "C": 222, "G": 190, "T": 126})

    def test_refs(self):
        pattern = HalfRelPattern.refs()
        self.assertIsInstance(pattern, HalfRelPattern)
        self.assertDictEqual(pattern.patterns,
                             {"A": 1, "C": 1, "G": 1, "T": 1})

    def test_none(self):
        pattern = HalfRelPattern.none()
        self.assertIsInstance(pattern, HalfRelPattern)
        self.assertDictEqual(pattern.patterns,
                             {"A": 0, "C": 0, "G": 0, "T": 0})

    def test_fits_refs(self):
        pattern = HalfRelPattern.refs()
        for base in "ACGT":
            for rel in range(256):
                with self.subTest(base=base, rel=rel):
                    self.assertEqual(pattern.fits(base, rel),
                                     rel == 0 or rel == 1)

    def test_fits_subs(self):
        pattern = HalfRelPattern.from_counts(count_sub=True)
        subs = {"A": [32, 64, 128],
                "C": [16, 64, 128],
                "G": [16, 32, 128],
                "T": [16, 32, 64]}
        fits = {base: {sum(subset) for subset in powerset(base_subs)}
                for base, base_subs in subs.items()}
        for base, base_fits in fits.items():
            for rel in range(256):
                with self.subTest(base=base, rel=rel):
                    self.assertEqual(pattern.fits(base, rel), rel in base_fits)

    def test_fits_muts(self):
        pattern = HalfRelPattern.muts()
        muts = {"A": [2, 4, 8, 32, 64, 128],
                "C": [2, 4, 8, 16, 64, 128],
                "G": [2, 4, 8, 16, 32, 128],
                "T": [2, 4, 8, 16, 32, 64]}
        fits = {base: {sum(subset) for subset in powerset(base_muts)}
                for base, base_muts in muts.items()}
        for base, base_fits in fits.items():
            for rel in range(256):
                with self.subTest(base=base, rel=rel):
                    self.assertEqual(pattern.fits(base, rel), rel in base_fits)

    def test_fits_allc(self):
        pattern = HalfRelPattern.allc()
        allc = {"A": [1, 2, 4, 8, 32, 64, 128],
                "C": [1, 2, 4, 8, 16, 64, 128],
                "G": [1, 2, 4, 8, 16, 32, 128],
                "T": [1, 2, 4, 8, 16, 32, 64]}
        fits = {base: {sum(subset) for subset in powerset(base_allc)}
                for base, base_allc in allc.items()}
        for base, base_fits in fits.items():
            for rel in range(256):
                with self.subTest(base=base, rel=rel):
                    self.assertEqual(pattern.fits(base, rel), rel in base_fits)

    def test_to_report_format(self):
        for _, _, plain, fancy in iter_codes():
            with self.subTest(fancy=fancy):
                pattern = HalfRelPattern(fancy)
                self.assertListEqual(pattern.to_report_format(),
                                     [fancy.upper()])

    def test_intersect(self):
        p1 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        p2 = HalfRelPattern("ac", "at", "ca", "cc", "ga", "gg", "gi", "td")
        expect = HalfRelPattern("ac", "ca", "cc", "ga", "gi", "td")
        intersect12 = p1.intersect(p2)
        self.assertIsInstance(intersect12, HalfRelPattern)
        intersect21 = p2.intersect(p1)
        self.assertIsInstance(intersect21, HalfRelPattern)
        self.assertEqual(intersect12, intersect21)
        self.assertEqual(intersect12, expect)

    def test_equal(self):
        p1 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        p2 = HalfRelPattern("ac", "at", "ca", "cc", "ga", "gg", "gi", "td")
        p3 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        self.assertNotEqual(p1, p2)
        self.assertEqual(p1, p3)
        self.assertNotEqual(p2, p3)


class TestRelPattern(ut.TestCase):

    def test_from_counts(self):
        for _, _, _, fancy in iter_codes():
            for count_del in [False, True]:
                for count_ins in [False, True]:
                    with self.subTest(fancy=fancy,
                                      count_del=count_del,
                                      count_ins=count_ins):
                        pattern = RelPattern.from_counts(count_del=count_del,
                                                         count_ins=count_ins,
                                                         discount=[fancy])
                        self.assertIsInstance(pattern, RelPattern)
                        self.assertEqual(
                            pattern.yes,
                            HalfRelPattern.from_counts(count_sub=True,
                                                       count_del=count_del,
                                                       count_ins=count_ins,
                                                       discount=[fancy])
                        )
                        self.assertEqual(
                            pattern.nos,
                            HalfRelPattern.from_counts(count_ref=True,
                                                       discount=[fancy])
                        )

    def test_allc(self):
        pattern = RelPattern.allc()
        self.assertEqual(pattern.yes, HalfRelPattern.allc())
        self.assertEqual(pattern.nos, HalfRelPattern.none())

    def test_muts(self):
        pattern = RelPattern.muts()
        self.assertEqual(pattern.yes, HalfRelPattern.muts())
        self.assertEqual(pattern.nos, HalfRelPattern.refs())

    def test_yes_nos(self):
        hrp1 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        hrp2 = HalfRelPattern("ac", "at", "ca", "cc", "ga", "gg", "gi", "td")
        pattern = RelPattern(hrp1, hrp2)
        self.assertIsInstance(pattern.yes, HalfRelPattern)
        self.assertEqual(pattern.yes, hrp1)
        self.assertNotEqual(pattern.yes, hrp2)
        self.assertIsInstance(pattern.nos, HalfRelPattern)
        self.assertEqual(pattern.nos, hrp2)
        self.assertNotEqual(pattern.nos, hrp1)

    def test_equal(self):
        hrp1 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        hrp2 = HalfRelPattern("ac", "at", "ca", "cc", "ga", "gg", "gi", "td")
        pattern1 = RelPattern(hrp1, hrp2)
        pattern2 = RelPattern(hrp2, hrp1)
        pattern3 = RelPattern(hrp1, hrp2)
        pattern4 = RelPattern(hrp1, hrp1)
        self.assertNotEqual(pattern1, pattern2)
        self.assertEqual(pattern1, pattern3)
        self.assertNotEqual(pattern2, pattern3)
        self.assertNotEqual(pattern1, pattern4)
        self.assertNotEqual(pattern2, pattern4)

    def test_fits(self):
        for hrp1, hrp2 in product([HalfRelPattern.none(),
                                   HalfRelPattern.refs(),
                                   HalfRelPattern.muts(),
                                   HalfRelPattern.allc()],
                                  repeat=2):
            pattern = RelPattern(hrp1, hrp2)
            for base in "ACGT":
                for rel in range(256):
                    with self.subTest(base=base, rel=rel):
                        info, fits = pattern.fits(base, rel)
                        if hrp1 is hrp2:
                            self.assertFalse(info)
                        else:
                            self.assertEqual(info, (hrp1.fits(base, rel)
                                                    != hrp2.fits(base, rel)))
                        self.assertEqual(fits, hrp1.fits(base, rel))

    def test_intersect_none(self):
        hrp1 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        hrp2 = HalfRelPattern("ac", "at", "ca", "cc", "ga", "gg", "gi", "td")
        self.assertEqual(RelPattern(hrp1, hrp2).intersect(None),
                         RelPattern(hrp1, hrp2))
        self.assertEqual(RelPattern(hrp1, hrp2).intersect(None, invert=True),
                         RelPattern(hrp2, hrp1))

    def test_intersect_other(self):
        hrp1 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        hrp2 = HalfRelPattern("ac", "at", "ca", "cc", "ga", "gg", "gi", "td")
        hrp12 = HalfRelPattern("ac", "ca", "cc", "ga", "gi", "td")
        hrp3 = HalfRelPattern("ag", "ca", "cc", "cg", "ct", "ga", "td", "ti")
        hrp4 = HalfRelPattern("at", "cg", "ga", "tc", "td")
        hrp34 = HalfRelPattern("cg", "ga", "td")
        self.assertEqual(
            RelPattern(hrp1, hrp3).intersect(RelPattern(hrp2, hrp4)),
            RelPattern(hrp12, hrp34)
        )
        self.assertEqual(
            RelPattern(hrp1, hrp3).intersect(RelPattern(hrp2, hrp4),
                                             invert=True),
            RelPattern(hrp34, hrp12)
        )

    def test_invert(self):
        hrp1 = HalfRelPattern("aa", "ac", "ca", "cc", "ga", "gt", "gi", "td")
        hrp2 = HalfRelPattern("ac", "at", "ca", "cc", "ga", "gg", "gi", "td")
        self.assertEqual(RelPattern(hrp1, hrp2).invert(),
                         RelPattern(hrp2, hrp1))


if __name__ == "__main__":
    ut.main()
