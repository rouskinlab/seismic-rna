import unittest as ut
from itertools import combinations, product

import numpy as np
import pandas as pd

from seismicrna.core.header import (AVERAGE_PREFIX,
                                    CLUSTER_PREFIX,
                                    CLUST_NAME,
                                    NUM_CLUSTS_NAME,
                                    REL_NAME,
                                    Header,
                                    RelHeader,
                                    ClustHeader,
                                    RelClustHeader,
                                    format_clust_name,
                                    format_clust_names,
                                    validate_k_clust,
                                    validate_ks,
                                    deduplicate_rels,
                                    list_clusts,
                                    list_k_clusts,
                                    list_ks_clusts,
                                    make_header,
                                    parse_header)


class TestConstants(ut.TestCase):

    def test_average_prefix(self):
        self.assertEqual(AVERAGE_PREFIX, "average")

    def test_cluster_prefix(self):
        self.assertEqual(CLUSTER_PREFIX, "cluster")

    def test_clust_name(self):
        self.assertEqual(CLUST_NAME, "Cluster")

    def test_k_name(self):
        self.assertEqual(NUM_CLUSTS_NAME, "K")

    def test_rel_name(self):
        self.assertEqual(REL_NAME, "Relationship")


class TestValidateKClust(ut.TestCase):

    def test_float_k(self):
        self.assertRaisesRegex(TypeError,
                               "k must be int",
                               validate_k_clust,
                               k=1.,
                               clust=1)

    def test_float_clust(self):
        self.assertRaisesRegex(TypeError,
                               "clust must be int",
                               validate_k_clust,
                               k=1,
                               clust=1.)

    def test_negative_zero(self):
        self.assertRaisesRegex(ValueError,
                               "k must be ≥ 0",
                               validate_k_clust,
                               k=-1,
                               clust=0)

    def test_zero_negative(self):
        self.assertRaisesRegex(ValueError,
                               "clust must be ≥ 0",
                               validate_k_clust,
                               k=0,
                               clust=-1)

    def test_zero(self):
        self.assertIsNone(validate_k_clust(0, 0))

    def test_zero_one_allowed(self):
        self.assertRaisesRegex(ValueError,
                               "clust must be ≤ k",
                               validate_k_clust,
                               k=0,
                               clust=1)

    def test_one_zero_allowed(self):
        self.assertRaisesRegex(ValueError,
                               "clust must be ≥ 1",
                               validate_k_clust,
                               k=1,
                               clust=0)

    def test_positive_positive(self):
        for k in range(1, 11):
            for clust in range(1, 11):
                if clust <= k:
                    self.assertIsNone(validate_k_clust(k, clust))
                else:
                    self.assertRaisesRegex(ValueError,
                                           "clust must be ≤ k",
                                           validate_k_clust,
                                           k=k,
                                           clust=clust)


class TestValidateKs(ut.TestCase):

    def test_empty(self):
        self.assertListEqual(validate_ks([]), [])

    def test_zero(self):
        self.assertRaisesRegex(ValueError,
                               "k must be ≥ 1, but got 0",
                               validate_ks,
                               [0])

    def test_valid(self):
        for min_k in range(1, 5):
            for max_k in range(min_k, 7):
                ks = list(range(min_k, max_k + 1))
                self.assertListEqual(validate_ks(ks), ks)
                self.assertListEqual(validate_ks(reversed(ks)), ks)

    def test_duplicated(self):
        self.assertRaisesRegex(ValueError,
                               "Duplicate k: 2",
                               validate_ks,
                               [1, 2, 3, 2, 1])


class TestDeduplicateRels(ut.TestCase):

    def test_empty(self):
        self.assertListEqual(deduplicate_rels([]), [])

    def test_no_duplicates(self):
        rels = "qwerty"
        self.assertListEqual(deduplicate_rels(rels), list(rels))

    def test_duplicates(self):
        rels = "qwqeerwtytyr"
        self.assertListEqual(deduplicate_rels(rels), list("qwerty"))


class TestFormatClustName(ut.TestCase):

    def test_zero(self):
        self.assertEqual(format_clust_name(0, 0),
                         AVERAGE_PREFIX)

    def test_positive(self):
        self.assertEqual(format_clust_name(1, 1), "cluster 1-1")
        self.assertEqual(format_clust_name(2, 1), "cluster 2-1")
        self.assertEqual(format_clust_name(2, 2), "cluster 2-2")


class TestFormatClustNames(ut.TestCase):

    def test_zero(self):
        self.assertListEqual(format_clust_names([(0, 0)]),
                             [AVERAGE_PREFIX])

    def test_positive_no_dups(self):
        self.assertListEqual(format_clust_names([(1, 1), (2, 1), (2, 2)]),
                             ["cluster 1-1", "cluster 2-1", "cluster 2-2"])

    def test_positive_valid_dups(self):
        self.assertListEqual(format_clust_names([(1, 1), (2, 1), (1, 1)],
                                                allow_duplicates=True),
                             ["cluster 1-1", "cluster 2-1", "cluster 1-1"])

    def test_positive_invalid_dups(self):
        self.assertRaisesRegex(ValueError,
                               "Duplicate clusters",
                               format_clust_names,
                               [(1, 1), (2, 1), (1, 1)])


class TestListClusts(ut.TestCase):

    def test_positive(self):
        self.assertListEqual(list_clusts(1), [1])
        self.assertListEqual(list_clusts(2), [1, 2])
        self.assertListEqual(list_clusts(3), [1, 2, 3])
        self.assertListEqual(list_clusts(4), [1, 2, 3, 4])

    def test_zero(self):
        self.assertRaisesRegex(ValueError,
                               "k must be ≥ 1, but got 0",
                               list_clusts,
                               0)


class TestListKClusts(ut.TestCase):

    def test_zero(self):
        self.assertRaisesRegex(ValueError,
                               "k must be ≥ 1, but got 0",
                               list_k_clusts,
                               0)

    def test_positive(self):
        self.assertListEqual(list_k_clusts(1),
                             [(1, 1)])
        self.assertListEqual(list_k_clusts(2),
                             [(2, 1), (2, 2)])
        self.assertListEqual(list_k_clusts(3),
                             [(3, 1), (3, 2), (3, 3)])
        self.assertListEqual(list_k_clusts(4),
                             [(4, 1), (4, 2), (4, 3), (4, 4)])


class TestListKsClusts(ut.TestCase):

    def test_empty(self):
        self.assertListEqual(list_ks_clusts([]), [])

    def test_zero(self):
        self.assertRaisesRegex(ValueError,
                               "k must be ≥ 1, but got 0",
                               list_ks_clusts,
                               [0])

    def test_one_k(self):
        for k in range(1, 10):
            self.assertListEqual(list_ks_clusts([k]),
                                 list_k_clusts(k))

    def test_ks(self):
        self.assertListEqual(list_ks_clusts([1, 2]),
                             [(1, 1), (2, 1), (2, 2)])
        self.assertListEqual(list_ks_clusts([1, 3]),
                             [(1, 1), (3, 1), (3, 2), (3, 3)])
        self.assertListEqual(list_ks_clusts([1, 2, 3]),
                             [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])


class TestHeader(ut.TestCase):

    def test_abstract(self):
        self.assertRaisesRegex(TypeError,
                               "Can't instantiate abstract class Header",
                               Header)


class TestRelHeader(ut.TestCase):

    def test_clustered(self):
        self.assertFalse(RelHeader.clustered())

    def test_levels(self):
        self.assertEqual(RelHeader.levels(), dict(rel=REL_NAME))

    def test_num_levels(self):
        self.assertEqual(RelHeader.num_levels(), 1)

    def test_level_keys(self):
        self.assertListEqual(RelHeader.level_keys(), ["rel"])

    def test_level_names(self):
        self.assertListEqual(RelHeader.level_names(), [REL_NAME])

    def test_rels_normal(self):
        rels = RelHeader(rels=list("qwerty")).rels
        self.assertListEqual(rels, list("qwerty"))

    def test_rels_duplicated(self):
        rels = RelHeader(rels=list("banana")).rels
        self.assertListEqual(rels, list("ban"))

    def test_rels_empty(self):
        rels = RelHeader(rels=[]).rels
        self.assertListEqual(rels, [])

    def test_ks(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertListEqual(header.ks, [0])

    def test_clusts(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertListEqual(header.clusts, [(0, 0)])

    def test_names(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertListEqual(header.names, [AVERAGE_PREFIX])

    def test_signature(self):
        rels = list("qwerty")
        sig = RelHeader(rels=rels).signature
        self.assertListEqual(list(sig.keys()), ["rels"])
        sig_rels = sig["rels"]
        self.assertIsInstance(sig_rels, list)
        self.assertListEqual(sig_rels, rels)

    def test_index(self):
        index = RelHeader(rels=list("qwerty")).index
        self.assertIsInstance(index, pd.Index)
        self.assertNotIsInstance(index, pd.MultiIndex)
        self.assertEqual(index.name, REL_NAME)
        self.assertListEqual(list(index.names), [REL_NAME])
        self.assertListEqual(index.to_list(), list("qwerty"))

    def test_iter_clust_indexes(self):
        header = RelHeader(rels=list("qwerty"))
        clust_indexes = list(header.iter_clust_indexes())
        self.assertEqual(len(clust_indexes), 1)
        self.assertEqual(len(clust_indexes), len(header.clusts))
        self.assertIsInstance(clust_indexes[0], pd.Index)
        self.assertTrue(clust_indexes[0].equals(header.index))

    def test_size(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.size, len("qwerty"))

    def test_select_none(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select()
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertTrue(selection.equals(header.index))

    def test_select_rel(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select(rel="w")
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(), ["w"])

    def test_select_one_rels(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select(rel=["w"])
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(), ["w"])

    def test_select_two_rels(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select(rel=["t", "w"])
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(), ["w", "t"])

    def test_select_invalid(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertRaisesRegex(ValueError,
                               "Expected rel to be one of",
                               header.select,
                               rel="v")

    def test_select_extra(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertRaisesRegex(TypeError,
                               "Unexpected keyword arguments",
                               header.select,
                               k=1)

    def test_select_extra_zero(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select(k=0)
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertTrue(selection.equals(header.index))

    def test_modified_none(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.modified(), header)

    def test_modified_rels(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.modified(rels=list("uiop")),
                         make_header(rels=list("uiop")))

    def test_modified_rels_empty(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.modified(rels=[]),
                         make_header(rels=[]))

    def test_modified_rels_none(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertRaisesRegex(TypeError,
                               "Must give rels, ks, or both, but got neither",
                               header.modified,
                               rels=None)

    def test_modified_ks(self):
        header = RelHeader(rels=list("qwerty"))
        for max_k in range(4):
            ks = list(range(1, max_k + 1))
            modified = header.modified(ks=ks)
            self.assertIsInstance(modified, RelClustHeader)
            self.assertEqual(modified, make_header(rels=list("qwerty"), ks=ks))


class TestClustHeader(ut.TestCase):

    def test_clustered(self):
        self.assertTrue(ClustHeader.clustered())

    def test_levels(self):
        self.assertEqual(ClustHeader.levels(),
                         dict(k=NUM_CLUSTS_NAME, clust=CLUST_NAME))

    def test_num_levels(self):
        self.assertEqual(ClustHeader.num_levels(), 2)

    def test_level_keys(self):
        self.assertListEqual(ClustHeader.level_keys(),
                             ["k", "clust"])

    def test_level_names(self):
        self.assertListEqual(ClustHeader.level_names(),
                             [NUM_CLUSTS_NAME, CLUST_NAME])

    def test_ks_valid(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = ClustHeader(ks=ks)
                self.assertListEqual(header.ks, ks)
                header = ClustHeader(ks=reversed(ks))
                self.assertListEqual(header.ks, ks)

    def test_ks_invalid(self):
        self.assertRaises(ValueError,
                          ClustHeader,
                          ks=[0])
        self.assertRaises(ValueError,
                          ClustHeader,
                          ks=[1, 1])

    def test_clusts(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = ClustHeader(ks=ks)
                self.assertListEqual(header.clusts,
                                     list_ks_clusts(ks))

    def test_names(self):
        header = ClustHeader(ks=[4, 3])
        self.assertListEqual(header.names,
                             ["cluster 3-1",
                              "cluster 3-2",
                              "cluster 3-3",
                              "cluster 4-1",
                              "cluster 4-2",
                              "cluster 4-3",
                              "cluster 4-4"])

    def test_signature(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                sig = ClustHeader(ks=ks).signature
                self.assertListEqual(list(sig.keys()), ["ks"])
                self.assertListEqual(sig["ks"], ks)

    def test_index(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                header = ClustHeader(ks=range(min_k, max_k + 1))
                index = header.index
                self.assertIsInstance(index, pd.MultiIndex)
                self.assertListEqual(list(index.names),
                                     [NUM_CLUSTS_NAME, CLUST_NAME])
                self.assertListEqual(index.tolist(), header.clusts)

    def test_iter_clust_indexes(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                header = ClustHeader(ks=range(min_k, max_k + 1))
                clust_indexes = list(header.iter_clust_indexes())
                self.assertEqual(len(clust_indexes), len(header.clusts))
                for index, clust in zip(clust_indexes,
                                        header.clusts,
                                        strict=True):
                    self.assertIsInstance(index, pd.MultiIndex)
                    self.assertListEqual(index.to_list(), [clust])

    def test_select_none(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select()
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertTrue(selection.equals(header.index))

    def test_select_k(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select(k=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [(3, 1), (3, 2), (3, 3)])

    def test_select_ks(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select(k=[3, 2])
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])

    def test_select_clust(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select(clust=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [(3, 3), (4, 3)])

    def test_select_clusts(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select(clust=[3, 2])
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [(2, 2), (3, 2), (3, 3), (4, 2), (4, 3)])

    def test_select_k_clust_exist(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select(k=3, clust=1)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [(3, 1)])

    def test_select_ks_clusts_exist(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select(k=[3, 2], clust=[1, 3])
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [(2, 1), (3, 1), (3, 3)])

    def test_select_k_clust_empty(self):
        header = ClustHeader(ks=range(1, 5))
        selection = header.select(k=1, clust=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [])

    def test_select_invalid_k(self):
        header = ClustHeader(ks=range(1, 5))
        self.assertRaisesRegex(ValueError,
                               "Expected k to be one of",
                               header.select,
                               k=5)

    def test_select_invalid_clust(self):
        header = ClustHeader(ks=range(1, 5))
        self.assertRaisesRegex(ValueError,
                               "Expected clust to be one of",
                               header.select,
                               clust=5)

    def test_select_extra(self):
        header = ClustHeader(ks=range(1, 5))
        self.assertRaisesRegex(TypeError,
                               "Unexpected keyword arguments",
                               header.select,
                               rel="w")

    def test_modified_none(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                header = ClustHeader(ks=range(min_k, max_k + 1))
                self.assertEqual(header.modified(), header)

    def test_modified_rels(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = ClustHeader(ks=ks)
                modified = header.modified(rels=list("qwerty"))
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified,
                                 make_header(rels=list("qwerty"), ks=ks))

    def test_modified_ks(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = ClustHeader(ks=ks)
                for new_min_k in range(1, 4):
                    for new_max_k in range(new_min_k, 6):
                        new_ks = list(range(new_min_k, new_max_k + 1))
                        modified = header.modified(ks=new_ks)
                        self.assertIsInstance(modified, ClustHeader)
                        self.assertEqual(modified, make_header(ks=new_ks))
                modified = header.modified(ks=[])
                self.assertIsInstance(modified, ClustHeader)
                self.assertEqual(modified, make_header(ks=[]))
                self.assertRaisesRegex(
                    TypeError,
                    "Must give rels, ks, or both, but got neither",
                    header.modified,
                    ks=None
                )


class TestRelClustHeader(ut.TestCase):

    def test_clustered(self):
        self.assertTrue(RelClustHeader.clustered())

    def test_levels(self):
        self.assertEqual(RelClustHeader.levels(),
                         dict(rel=REL_NAME,
                              k=NUM_CLUSTS_NAME,
                              clust=CLUST_NAME))

    def test_num_levels(self):
        self.assertEqual(RelClustHeader.num_levels(), 3)

    def test_level_keys(self):
        self.assertListEqual(RelClustHeader.level_keys(),
                             ["rel", "k", "clust"])

    def test_level_names(self):
        self.assertListEqual(RelClustHeader.level_names(),
                             [REL_NAME, NUM_CLUSTS_NAME, CLUST_NAME])

    def test_signature(self):
        rels = list("qwerty")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                sig = RelClustHeader(rels=rels, ks=ks).signature
                self.assertListEqual(sorted(sig.keys()), ["ks", "rels"])
                self.assertListEqual(sig["ks"], ks)
                self.assertListEqual(sig["rels"], rels)

    def test_ks(self):
        rels = list("qwerty")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=rels, ks=ks)
                self.assertListEqual(header.ks, ks)

    def test_clusts(self):
        rels = list("qwerty")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=rels, ks=ks)
                self.assertListEqual(header.clusts, list_ks_clusts(ks))

    def test_index(self):
        index = RelClustHeader(rels=["a", "b"], ks=[2, 3]).index
        self.assertIsInstance(index, pd.MultiIndex)
        self.assertListEqual(list(index.names),
                             [REL_NAME, NUM_CLUSTS_NAME, CLUST_NAME])
        self.assertTrue(np.array_equal(index.get_level_values(REL_NAME),
                                       list("aaaaabbbbb")))
        self.assertTrue(np.array_equal(index.get_level_values(NUM_CLUSTS_NAME),
                                       [2, 2, 3, 3, 3, 2, 2, 3, 3, 3]))
        self.assertTrue(np.array_equal(index.get_level_values(CLUST_NAME),
                                       [1, 2, 1, 2, 3, 1, 2, 1, 2, 3]))

    def test_iter_clust_indexes(self):
        rels = list("qwerty")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=rels, ks=ks)
                clust_indexes = list(header.iter_clust_indexes())
                self.assertEqual(len(clust_indexes), len(header.clusts))
                for index, clust in zip(clust_indexes,
                                        header.clusts,
                                        strict=True):
                    self.assertIsInstance(index, pd.MultiIndex)
                    self.assertEqual(index.size, len(rels))
                    self.assertListEqual(index.to_list(),
                                         [(rel, *clust) for rel in rels])

    def test_select_none(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select()
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertTrue(selection.equals(header.index))

    def test_select_rel(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(rel="b")
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [("b", 2, 1),
                              ("b", 2, 2),
                              ("b", 3, 1),
                              ("b", 3, 2),
                              ("b", 3, 3)])

    def test_select_ks(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(k=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [("a", 3, 1),
                              ("a", 3, 2),
                              ("a", 3, 3),
                              ("b", 3, 1),
                              ("b", 3, 2),
                              ("b", 3, 3)])

    def test_select_clust(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(clust=2)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [("a", 2, 2),
                              ("a", 3, 2),
                              ("b", 2, 2),
                              ("b", 3, 2)])

    def test_select_k_clust_exist(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(rel="a", k=2, clust=1)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [("a", 2, 1)])

    def test_select_k_clust_empty(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(rel="a", k=2, clust=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [])

    def test_select_invalid_rel(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        self.assertRaisesRegex(ValueError,
                               "Expected rel to be one of",
                               header.select,
                               rel="c")

    def test_select_invalid_k(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        self.assertRaisesRegex(ValueError,
                               "Expected k to be one of",
                               header.select,
                               k=1)

    def test_select_invalid_clust(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        self.assertRaisesRegex(ValueError,
                               "Expected clust to be one of",
                               header.select,
                               clust=4)

    def test_select_extra(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        self.assertRaisesRegex(TypeError,
                               "Unexpected keyword arguments",
                               header.select,
                               extra="x")

    def test_select_extra_none(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(rel="a", k=2, clust=1, extra=None)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [("a", 2, 1)])

    def test_select_extra_zero(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(rel="a", k=2, clust=1, extra=0)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [("a", 2, 1)])

    def test_select_extra_emptystr(self):
        header = RelClustHeader(rels=["a", "b"], ks=[2, 3])
        selection = header.select(rel="a", k=2, clust=1, extra="")
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertListEqual(selection.to_list(),
                             [("a", 2, 1)])

    def test_modified_none(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                header = RelClustHeader(rels=list("qwerty"),
                                        ks=range(min_k, max_k + 1))
                self.assertEqual(header.modified(), header)

    def test_modified_rels(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=list("qwerty"), ks=ks)
                modified = header.modified(rels=list("uiop"))
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified,
                                 make_header(rels=list("uiop"), ks=ks))

    def test_modified_rels_empty(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=list("qwerty"), ks=ks)
                modified = header.modified(rels=[])
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified, make_header(rels=[], ks=ks))

    def test_modified_rels_none(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=list("qwerty"), ks=ks)
                modified = header.modified(rels=None)
                self.assertIsInstance(modified, ClustHeader)
                self.assertEqual(modified, make_header(ks=ks))

    def test_modified_ks(self):
        rels = list("qwerty")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=rels, ks=ks)
                for new_min_k in range(1, 4):
                    for new_max_k in range(new_min_k, 6):
                        new_ks = list(range(new_min_k, new_max_k + 1))
                        modified = header.modified(ks=new_ks)
                        self.assertIsInstance(modified, RelClustHeader)
                        self.assertEqual(modified,
                                         make_header(rels=rels, ks=new_ks))

    def test_modified_ks_empty(self):
        rels = list("qwerty")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=rels, ks=ks)
                modified = header.modified(ks=[])
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified,
                                 make_header(rels=rels, ks=[]))

    def test_modified_ks_none(self):
        rels = list("qwerty")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                ks = list(range(min_k, max_k + 1))
                header = RelClustHeader(rels=rels, ks=ks)
                modified = header.modified(ks=None)
                self.assertIsInstance(modified, RelHeader)
                self.assertEqual(modified, make_header(rels=rels))

    def test_modified_all(self):
        new_rels = list("uiop")
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                header = RelClustHeader(rels=list("qwerty"),
                                        ks=range(min_k, max_k + 1))
                for new_min_k in range(1, 4):
                    for new_max_k in range(new_min_k, 6):
                        new_ks = list(range(new_min_k, new_max_k + 1))
                        modified = header.modified(rels=new_rels, ks=new_ks)
                        self.assertIsInstance(modified, ClustHeader)
                        self.assertEqual(modified,
                                         make_header(rels=new_rels, ks=new_ks))

    def test_modified_nullified(self):
        for min_k in range(1, 4):
            for max_k in range(min_k, 6):
                header = RelClustHeader(rels=list("qwerty"),
                                        ks=range(min_k, max_k + 1))
                self.assertRaisesRegex(
                    TypeError,
                    "Must give rels, ks, or both, but got neither",
                    header.modified,
                    rels=None,
                    ks=None
                )


class TestEqualHeaders(ut.TestCase):

    def test_relheaders(self):
        for rels1, rels2 in product(["qwerty", "uiop", "asdf", "ghjkl"],
                                    repeat=2):
            header1 = RelHeader(rels=list(rels1))
            header2 = RelHeader(rels=list(rels2))
            if rels1 == rels2:
                self.assertEqual(header1, header2)
            else:
                self.assertNotEqual(header1, header2)

    def test_clustheaders(self):
        for min_k1 in range(1, 4):
            for min_k2 in range(1, 4):
                for max_k1 in range(min_k1, 4):
                    for max_k2 in range(min_k2, 4):
                        header1 = ClustHeader(ks=range(min_k1, max_k1 + 1))
                        header2 = ClustHeader(ks=range(min_k2, max_k2 + 1))
                        if min_k1 == min_k2 and max_k1 == max_k2:
                            self.assertEqual(header1, header2)
                        else:
                            self.assertNotEqual(header1, header2)

    def test_relclustheaders(self):
        for min_k1 in range(1, 4):
            for min_k2 in range(1, 4):
                for max_k1 in range(min_k1, 4):
                    for max_k2 in range(min_k2, 4):
                        for rels1, rels2 in product(["qwerty",
                                                     "uiop",
                                                     "asdf",
                                                     "ghjkl"],
                                                    repeat=2):
                            header1 = RelClustHeader(rels=list(rels1),
                                                     ks=range(min_k1,
                                                              max_k1 + 1))
                            header2 = RelClustHeader(rels=list(rels2),
                                                     ks=range(min_k2,
                                                              max_k2 + 1))
                            if (rels1 == rels2
                                    and min_k1 == min_k2
                                    and max_k1 == max_k2):
                                self.assertEqual(header1, header2)
                            else:
                                self.assertNotEqual(header1, header2)

    def test_different_types(self):
        for rels in ["qwerty", "uiop", "asdf", "ghjkl"]:
            for min_k in range(1, 4):
                for max_k in range(min_k, 6):
                    ks = list(range(min_k, max_k + 1))
                    rh = RelHeader(rels=list(rels))
                    ch = ClustHeader(ks=ks)
                    rch = RelClustHeader(rels=list(rels), ks=ks)
                    for header1, header2 in combinations([rh, ch, rch], 2):
                        self.assertNotEqual(header1, header2)


class TestMakeHeader(ut.TestCase):

    def test_none(self):
        self.assertRaisesRegex(TypeError,
                               "Must give rels, ks, or both, but got neither",
                               make_header)

    def test_rels(self):
        rels = ["a", "b"]
        header = make_header(rels=rels)
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.tolist(), rels)

    def test_ks(self):
        header = make_header(ks=[3, 2])
        self.assertIsInstance(header, ClustHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.tolist(),
                             [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])

    def test_all(self):
        header = make_header(rels=["a", "b"], ks=[3, 2])
        self.assertIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.tolist(),
                             [("a", 2, 1),
                              ("a", 2, 2),
                              ("a", 3, 1),
                              ("a", 3, 2),
                              ("a", 3, 3),
                              ("b", 2, 1),
                              ("b", 2, 2),
                              ("b", 3, 1),
                              ("b", 3, 2),
                              ("b", 3, 3)])


class TestParseHeader(ut.TestCase):

    def test_empty(self):
        header = parse_header(pd.Index([]))
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.to_list(), [])

    def test_rel_index(self):
        header = parse_header(pd.Index(["a", "b"]))
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.to_list(), ["a", "b"])

    def test_rel_index_valid_name(self):
        header = parse_header(pd.Index(["a", "b"], name=REL_NAME))
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.to_list(), ["a", "b"])

    def test_rel_index_invalid_name(self):
        self.assertRaisesRegex(ValueError,
                               f"Expected index named {repr(REL_NAME)}",
                               parse_header,
                               pd.Index(["a", "b"], name=NUM_CLUSTS_NAME))

    def test_rel_multiindex(self):
        header = parse_header(pd.MultiIndex.from_arrays([["a", "b"]],
                                                        names=[REL_NAME]))
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.to_list(), ["a", "b"])

    def test_clust(self):
        header = parse_header(pd.MultiIndex.from_tuples([("1", "1"),
                                                         ("2", "1"),
                                                         ("2", "2")],
                                                        names=[NUM_CLUSTS_NAME,
                                                               CLUST_NAME]))
        self.assertIsInstance(header, ClustHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.to_list(),
                             [(1, 1), (2, 1), (2, 2)])

    def test_relclust(self):
        header = parse_header(pd.MultiIndex.from_tuples([("b", "1", "1"),
                                                         ("b", "2", "1"),
                                                         ("b", "2", "2"),
                                                         ("a", "1", "1"),
                                                         ("a", "2", "1"),
                                                         ("a", "2", "2")],
                                                        names=[REL_NAME,
                                                               NUM_CLUSTS_NAME,
                                                               CLUST_NAME]))
        self.assertIsInstance(header, RelClustHeader)
        self.assertListEqual(header.index.to_list(),
                             [("b", 1, 1),
                              ("b", 2, 1),
                              ("b", 2, 2),
                              ("a", 1, 1),
                              ("a", 2, 1),
                              ("a", 2, 2)])

    def test_rel_index_repeated(self):
        self.assertRaisesRegex(ValueError,
                               "Indexes do not match",
                               parse_header,
                               pd.Index(["a", "b", "a"]))

    def test_missing_index_names(self):
        self.assertRaisesRegex(TypeError,
                               "Must give rels, ks, or both, but got neither",
                               parse_header,
                               pd.MultiIndex.from_arrays([["a", "b"]],
                                                         names=["random"]))

    def test_extra_index_names(self):
        self.assertRaisesRegex(ValueError,
                               "Expected index names",
                               parse_header,
                               pd.MultiIndex.from_arrays([["a", "b"],
                                                          [0, 0]],
                                                         names=[REL_NAME,
                                                                "random"]))

    def test_missing_values(self):
        self.assertRaisesRegex(ValueError,
                               "Invalid index level",
                               parse_header,
                               pd.MultiIndex.from_tuples([("1", "1"),
                                                          ("2", "2")],
                                                         names=[NUM_CLUSTS_NAME,
                                                                CLUST_NAME]))

    def test_extra_values(self):
        self.assertRaisesRegex(ValueError,
                               "Invalid index level",
                               parse_header,
                               pd.MultiIndex.from_tuples([("1", "1"),
                                                          ("2", "1"),
                                                          ("2", "1"),
                                                          ("2", "2")],
                                                         names=[NUM_CLUSTS_NAME,
                                                                CLUST_NAME]))

    def test_nonnumeric(self):
        self.assertRaisesRegex(ValueError,
                               "invalid literal for int",
                               parse_header,
                               pd.MultiIndex.from_tuples([("a", "x", "x"),
                                                          ("a", "y", "x"),
                                                          ("a", "y", "y"),
                                                          ("b", "x", "x"),
                                                          ("b", "y", "x"),
                                                          ("b", "y", "y")],
                                                         names=[REL_NAME,
                                                                NUM_CLUSTS_NAME,
                                                                CLUST_NAME]))

    def test_invalid_numeric(self):
        self.assertRaisesRegex(ValueError,
                               "k must be ≥ 1, but got 0",
                               parse_header,
                               pd.MultiIndex.from_tuples([("a", "0", "0"),
                                                          ("b", "0", "0")],
                                                         names=[REL_NAME,
                                                                NUM_CLUSTS_NAME,
                                                                CLUST_NAME]))


if __name__ == "__main__":
    ut.main()

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
