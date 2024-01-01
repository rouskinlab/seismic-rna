import unittest as ut
from itertools import combinations, product

import numpy as np
import pandas as pd

from seismicrna.core.header import (AVERAGE_PREFIX,
                                    CLUSTER_PREFIX,
                                    CLUST_NAME,
                                    ORDER_NAME,
                                    REL_NAME,
                                    Header,
                                    RelHeader,
                                    ClustHeader,
                                    RelClustHeader,
                                    format_clust_name,
                                    format_clust_names,
                                    validate_order_clust,
                                    list_clusts,
                                    list_orders,
                                    list_order_clusts,
                                    list_orders_clusts,
                                    index_clusts,
                                    index_orders,
                                    index_order_clusts,
                                    index_orders_clusts,
                                    make_header,
                                    parse_header)


class TestConstants(ut.TestCase):

    def test_average_prefix(self):
        self.assertEqual(AVERAGE_PREFIX, "average")

    def test_cluster_prefix(self):
        self.assertEqual(CLUSTER_PREFIX, "cluster")

    def test_clust_name(self):
        self.assertEqual(CLUST_NAME, "Cluster")

    def test_order_name(self):
        self.assertEqual(ORDER_NAME, "Order")

    def test_rel_name(self):
        self.assertEqual(REL_NAME, "Relationship")


class TestValidateOrderClust(ut.TestCase):

    def test_float_order(self):
        self.assertRaisesRegex(TypeError,
                               "order must be an int",
                               validate_order_clust,
                               order=1.,
                               clust=1)

    def test_float_clust(self):
        self.assertRaisesRegex(TypeError,
                               "clust must be an int",
                               validate_order_clust,
                               order=1,
                               clust=1.)

    def test_negative_zero_allowed(self):
        self.assertRaisesRegex(ValueError,
                               "order must be ≥ 0",
                               validate_order_clust,
                               order=-1,
                               clust=0,
                               allow_zero=True)

    def test_zero_negative_allowed(self):
        self.assertRaisesRegex(ValueError,
                               "clust must be ≥ 0",
                               validate_order_clust,
                               order=0,
                               clust=-1,
                               allow_zero=True)

    def test_zero_zero_allowed(self):
        self.assertIsNone(validate_order_clust(0, 0, allow_zero=True))

    def test_zero_zero_unallowed(self):
        self.assertRaisesRegex(ValueError,
                               "order must be ≥ 1",
                               validate_order_clust,
                               order=0,
                               clust=0,
                               allow_zero=False)

    def test_zero_one_allowed(self):
        self.assertRaisesRegex(ValueError,
                               "clust must be ≤ order",
                               validate_order_clust,
                               order=0,
                               clust=1,
                               allow_zero=True)

    def test_zero_one_unallowed(self):
        self.assertRaisesRegex(ValueError,
                               "order must be ≥ 1",
                               validate_order_clust,
                               order=0,
                               clust=1,
                               allow_zero=False)

    def test_one_zero_allowed(self):
        self.assertRaisesRegex(ValueError,
                               "clust must be ≥ 1",
                               validate_order_clust,
                               order=1,
                               clust=0,
                               allow_zero=True)

    def test_one_zero_unallowed(self):
        self.assertRaisesRegex(ValueError,
                               "clust must be ≥ 1",
                               validate_order_clust,
                               order=1,
                               clust=0,
                               allow_zero=False)

    def test_positive_positive(self):
        for order in range(1, 11):
            for clust in range(1, 11):
                if clust <= order:
                    self.assertIsNone(validate_order_clust(order, clust))
                else:
                    self.assertRaisesRegex(ValueError,
                                           "clust must be ≤ order",
                                           validate_order_clust,
                                           order=order,
                                           clust=clust)


class TestFormatClustName(ut.TestCase):

    def test_zero_zero_allowed(self):
        self.assertEqual(format_clust_name(0, 0, allow_zero=True),
                         AVERAGE_PREFIX)

    def test_positive(self):
        self.assertEqual(format_clust_name(1, 1), "cluster 1-1")
        self.assertEqual(format_clust_name(2, 1), "cluster 2-1")
        self.assertEqual(format_clust_name(2, 2), "cluster 2-2")


class TestFormatClustNames(ut.TestCase):

    def test_zero_zero_allowed(self):
        self.assertEqual(format_clust_names([(0, 0)],
                                            allow_zero=True),
                         [AVERAGE_PREFIX])

    def test_positive_no_dups(self):
        self.assertEqual(format_clust_names([(1, 1), (2, 1), (2, 2)]),
                         ["cluster 1-1", "cluster 2-1", "cluster 2-2"])

    def test_positive_valid_dups(self):
        self.assertEqual(format_clust_names([(1, 1), (2, 1), (1, 1)],
                                            allow_duplicates=True),
                         ["cluster 1-1", "cluster 2-1", "cluster 1-1"])

    def test_positive_invalid_dups(self):
        self.assertRaisesRegex(ValueError,
                               "Duplicate clusters",
                               format_clust_names,
                               [(1, 1), (2, 1), (1, 1)])


class TestListClusts(ut.TestCase):

    def test_zero(self):
        self.assertEqual(list_clusts(0), [0])

    def test_positive(self):
        self.assertEqual(list_clusts(1), [1])
        self.assertEqual(list_clusts(2), [1, 2])
        self.assertEqual(list_clusts(3), [1, 2, 3])
        self.assertEqual(list_clusts(4), [1, 2, 3, 4])

    def test_negative(self):
        self.assertRaisesRegex(ValueError,
                               "order must be ≥ 0",
                               list_clusts,
                               -1)


class TestListOrders(ut.TestCase):

    def test_zero_none(self):
        self.assertEqual(list_orders(0), [0])

    def test_zero_zero(self):
        self.assertRaisesRegex(ValueError,
                               "min_order must be ≥ 1",
                               list_orders,
                               0,
                               0)

    def test_zero_one(self):
        self.assertEqual(list_orders(0, 1),
                         list_orders(0))

    def test_zero_two(self):
        self.assertRaisesRegex(ValueError,
                               "If max_order = 0, then min_order must be 1",
                               list_orders,
                               0,
                               2)

    def test_positive_none(self):
        self.assertEqual(list_orders(1), [1])
        self.assertEqual(list_orders(2), [1, 2])
        self.assertEqual(list_orders(3), [1, 2, 3])
        self.assertEqual(list_orders(4), [1, 2, 3, 4])

    def test_positive_zero(self):
        self.assertRaisesRegex(ValueError,
                               "min_order must be ≥ 1",
                               list_orders,
                               1,
                               0)

    def test_positive_positive(self):
        for max_order in range(1, 11):
            for min_order in range(1, 11):
                if max_order >= min_order:
                    self.assertEqual(list_orders(max_order, min_order),
                                     list(range(min_order, max_order + 1)))
                else:
                    self.assertRaisesRegex(ValueError,
                                           "If not 0, max_order must be ≥",
                                           list_orders,
                                           max_order,
                                           min_order)

    def test_negative_none(self):
        self.assertRaisesRegex(ValueError,
                               "max_order must be ≥ 0",
                               list_orders,
                               -1)

    def test_negative_zero(self):
        self.assertRaisesRegex(ValueError,
                               "min_order must be ≥ 1",
                               list_orders,
                               -1,
                               0)

    def test_negative_positive(self):
        for min_order in range(1, 11):
            self.assertRaisesRegex(ValueError,
                                   "max_order must be ≥ 0",
                                   list_orders,
                                   -1,
                                   min_order)


class TestIndexClusts(ut.TestCase):

    def test_valid(self):
        for order in range(11):
            index = index_clusts(order)
            self.assertEqual(index.to_list(),
                             list_clusts(order))
            self.assertEqual(index.names,
                             [CLUST_NAME])


class TestIndexOrders(ut.TestCase):

    def test_valid(self):
        for max_order in range(11):
            for min_order in range(1, max_order + 1):
                index = index_orders(max_order, min_order)
                self.assertEqual(index.to_list(),
                                 list_orders(max_order, min_order))
                self.assertEqual(index.names,
                                 [ORDER_NAME])


class TestListOrderClusts(ut.TestCase):

    def test_zero(self):
        self.assertEqual(list_order_clusts(0),
                         [(0, 0)])

    def test_positive(self):
        self.assertEqual(list_order_clusts(1),
                         [(1, 1)])
        self.assertEqual(list_order_clusts(2),
                         [(2, 1), (2, 2)])
        self.assertEqual(list_order_clusts(3),
                         [(3, 1), (3, 2), (3, 3)])
        self.assertEqual(list_order_clusts(4),
                         [(4, 1), (4, 2), (4, 3), (4, 4)])


class TestListOrdersClusts(ut.TestCase):

    def test_zero_none(self):
        self.assertEqual(list_orders_clusts(0),
                         [(0, 0)])

    def test_zero_one(self):
        self.assertEqual(list_orders_clusts(0, 1),
                         list_orders_clusts(0))

    def test_positive_none(self):
        self.assertEqual(list_orders_clusts(1),
                         [(1, 1)])
        self.assertEqual(list_orders_clusts(2),
                         [(1, 1), (2, 1), (2, 2)])
        self.assertEqual(list_orders_clusts(3),
                         [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])

    def test_positive_positive(self):
        for max_order in range(1, 11):
            for min_order in range(1, max_order + 1):
                expect = list()
                for order in range(min_order, max_order + 1):
                    expect.extend(list_order_clusts(order))
                self.assertEqual(list_orders_clusts(max_order, min_order),
                                 expect)


class TestIndexOrderClusts(ut.TestCase):

    def test_valid(self):
        for order in range(11):
            index = index_order_clusts(order)
            self.assertEqual(index.to_list(),
                             list_order_clusts(order))
            self.assertEqual(index.names,
                             [ORDER_NAME, CLUST_NAME])


class TestIndexOrdersClusts(ut.TestCase):

    def test_valid(self):
        for max_order in range(1, 11):
            for min_order in range(1, max_order + 1):
                index = index_orders_clusts(max_order, min_order)
                self.assertEqual(index.to_list(),
                                 list_orders_clusts(max_order, min_order))
                self.assertEqual(index.names,
                                 [ORDER_NAME, CLUST_NAME])


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
        self.assertEqual(RelHeader.level_keys(), ["rel"])

    def test_level_names(self):
        self.assertEqual(RelHeader.level_names(), [REL_NAME])

    def test_rels_normal(self):
        rels = RelHeader(rels=list("qwerty")).rels
        self.assertIsInstance(rels, np.ndarray)
        self.assertTrue(np.array_equal(rels, np.array(list("qwerty"))))

    def test_rels_duplicated(self):
        rels = RelHeader(rels=list("banana")).rels
        self.assertIsInstance(rels, np.ndarray)
        self.assertTrue(np.array_equal(rels, np.array(list("ban"))))

    def test_rels_empty(self):
        self.assertRaisesRegex(ValueError,
                               "Got no relationships for header",
                               RelHeader,
                               rels=list())

    def test_max_order(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.max_order, 0)

    def test_min_order(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.min_order, 1)

    def test_orders(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.orders.to_list(), [0])

    def test_clusts(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.clusts.to_list(), [(0, 0)])

    def test_names(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertEqual(header.names, [AVERAGE_PREFIX])

    def test_signature(self):
        sig = RelHeader(rels=list("qwerty")).signature
        self.assertEqual(list(sig.keys()), ["max_order", "min_order", "rels"])
        sig_max = sig["max_order"]
        self.assertIsInstance(sig_max, int)
        self.assertEqual(sig_max, 0)
        sig_min = sig["min_order"]
        self.assertIsInstance(sig_min, int)
        self.assertEqual(sig_min, 1)
        sig_rels = sig["rels"]
        self.assertIsInstance(sig_rels, np.ndarray)
        self.assertEqual(sig_rels.tolist(), ["q", "w", "e", "r", "t", "y"])

    def test_index(self):
        index = RelHeader(rels=list("qwerty")).index
        self.assertIsInstance(index, pd.Index)
        self.assertNotIsInstance(index, pd.MultiIndex)
        self.assertEqual(index.name, REL_NAME)
        self.assertEqual(list(index.names), [REL_NAME])
        self.assertEqual(index.to_list(), list("qwerty"))

    def test_iter_clust_indexes(self):
        header = RelHeader(rels=list("qwerty"))
        clust_indexes = list(header.iter_clust_indexes())
        self.assertEqual(len(clust_indexes), 1)
        self.assertEqual(len(clust_indexes), header.clusts.size)
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
        selection = header.select(rel='w')
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(), ['w'])

    def test_select_one_rels(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select(rel=['w'])
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(), ['w'])

    def test_select_two_rels(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select(rel=['t', 'w'])
        self.assertIsInstance(selection, pd.Index)
        self.assertNotIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(), ['w', 't'])

    def test_select_invalid(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertRaisesRegex(ValueError,
                               "Expected rel to be one of",
                               header.select,
                               rel='v')

    def test_select_extra(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertRaisesRegex(TypeError,
                               "Unexpected keyword arguments",
                               header.select,
                               order=1)

    def test_select_extra_zero(self):
        header = RelHeader(rels=list("qwerty"))
        selection = header.select(order=0)
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

    def test_modified_empty_rels(self):
        header = RelHeader(rels=list("qwerty"))
        self.assertRaisesRegex(TypeError,
                               "No header for rels",
                               header.modified,
                               rels=[])

    def test_modified_max_order(self):
        header = RelHeader(rels=list("qwerty"))
        for max_order in range(4):
            modified = header.modified(max_order=max_order)
            if max_order == 0:
                self.assertIsInstance(modified, RelHeader)
                self.assertEqual(modified, header)
            else:
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified, make_header(rels=list("qwerty"),
                                                       max_order=max_order))

    def test_modified_min_order(self):
        header = RelHeader(rels=list("qwerty"))
        for min_order in range(4):
            modified = header.modified(max_order=min_order, min_order=min_order)
            if min_order == 0:
                self.assertIsInstance(modified, RelHeader)
                self.assertEqual(modified, header)
            else:
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified, make_header(rels=list("qwerty"),
                                                       max_order=min_order,
                                                       min_order=min_order))


class TestClustHeader(ut.TestCase):

    def test_clustered(self):
        self.assertTrue(ClustHeader.clustered())

    def test_levels(self):
        self.assertEqual(ClustHeader.levels(),
                         dict(order=ORDER_NAME, clust=CLUST_NAME))

    def test_num_levels(self):
        self.assertEqual(ClustHeader.num_levels(), 2)

    def test_level_keys(self):
        self.assertEqual(ClustHeader.level_keys(), ["order", "clust"])

    def test_level_names(self):
        self.assertEqual(ClustHeader.level_names(), [ORDER_NAME, CLUST_NAME])

    def test_max_order(self):
        for max_order in range(1, 11):
            header = ClustHeader(max_order=max_order)
            self.assertEqual(header.max_order, max_order)
        self.assertRaisesRegex(ValueError,
                               "max_order must be ≥ 1, but got 0",
                               ClustHeader,
                               max_order=0)

    def test_min_order(self):
        for min_order in range(1, 11):
            header = ClustHeader(max_order=min_order, min_order=min_order)
            self.assertEqual(header.min_order, min_order)

    def test_orders(self):
        for max_order in range(1, 11):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                self.assertTrue(header.orders.equals(index_orders(max_order,
                                                                  min_order)))

    def test_clusts(self):
        for max_order in range(1, 11):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                self.assertTrue(header.clusts.equals(
                    index_orders_clusts(max_order, min_order)
                ))

    def test_names(self):
        header = ClustHeader(max_order=4, min_order=3)
        self.assertEqual(header.names, ["cluster 3-1",
                                        "cluster 3-2",
                                        "cluster 3-3",
                                        "cluster 4-1",
                                        "cluster 4-2",
                                        "cluster 4-3",
                                        "cluster 4-4"])

    def test_signature(self):
        for max_order in range(1, 11):
            for min_order in range(1, max_order + 1):
                sig = ClustHeader(max_order=max_order,
                                  min_order=min_order).signature
                self.assertEqual(list(sig.keys()), ["max_order", "min_order"])
                sig_max = sig["max_order"]
                self.assertIsInstance(sig_max, int)
                self.assertEqual(sig_max, max_order)
                sig_min = sig["min_order"]
                self.assertIsInstance(sig_min, int)
                self.assertEqual(sig_min, min_order)

    def test_index(self):
        for max_order in range(1, 11):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                index = header.index
                self.assertIsInstance(index, pd.MultiIndex)
                self.assertEqual(list(index.names), [ORDER_NAME, CLUST_NAME])
                self.assertTrue(index.equals(header.clusts))

    def test_iter_clust_indexes(self):
        for max_order in range(1, 6):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                clust_indexes = list(header.iter_clust_indexes())
                self.assertEqual(len(clust_indexes), header.clusts.size)
                for index, clust in zip(clust_indexes,
                                        header.clusts,
                                        strict=True):
                    self.assertIsInstance(index, pd.MultiIndex)
                    self.assertEqual(index.to_list(), [clust])

    def test_select_none(self):
        header = ClustHeader(max_order=4)
        selection = header.select()
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertTrue(selection.equals(header.index))

    def test_select_order(self):
        header = ClustHeader(max_order=4)
        selection = header.select(order=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [(3, 1), (3, 2), (3, 3)])

    def test_select_orders(self):
        header = ClustHeader(max_order=4)
        selection = header.select(order=[3, 2])
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])

    def test_select_clust(self):
        header = ClustHeader(max_order=4)
        selection = header.select(clust=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [(3, 3), (4, 3)])

    def test_select_clusts(self):
        header = ClustHeader(max_order=4)
        selection = header.select(clust=[3, 2])
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [(2, 2), (3, 2), (3, 3), (4, 2), (4, 3)])

    def test_select_order_clust_exist(self):
        header = ClustHeader(max_order=4)
        selection = header.select(order=3, clust=1)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [(3, 1)])

    def test_select_orders_clusts_exist(self):
        header = ClustHeader(max_order=4)
        selection = header.select(order=[3, 2], clust=[1, 3])
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [(2, 1), (3, 1), (3, 3)])

    def test_select_order_clust_empty(self):
        header = ClustHeader(max_order=4)
        selection = header.select(order=1, clust=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [])

    def test_select_invalid_order(self):
        header = ClustHeader(max_order=4)
        self.assertRaisesRegex(ValueError,
                               "Expected order to be one of",
                               header.select,
                               order=5)

    def test_select_invalid_clust(self):
        header = ClustHeader(max_order=4)
        self.assertRaisesRegex(ValueError,
                               "Expected clust to be one of",
                               header.select,
                               clust=5)

    def test_select_extra(self):
        header = ClustHeader(max_order=4)
        self.assertRaisesRegex(TypeError,
                               "Unexpected keyword arguments",
                               header.select,
                               rel='w')

    def test_modified_none(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                self.assertEqual(header.modified(), header)

    def test_modified_rels(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                modified = header.modified(rels=list("qwerty"))
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified,
                                 make_header(rels=list("qwerty"),
                                             max_order=max_order,
                                             min_order=min_order))

    def test_modified_max_order(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                for new_max_order in range(1, 4):
                    modified = header.modified(max_order=new_max_order)
                    self.assertIsInstance(modified, ClustHeader)
                    self.assertEqual(modified,
                                     make_header(max_order=new_max_order,
                                                 min_order=min_order))
                self.assertRaisesRegex(TypeError,
                                       "No header for rels",
                                       header.modified,
                                       max_order=0)

    def test_modified_min_order(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = ClustHeader(max_order=max_order, min_order=min_order)
                for new_min_order in range(1, 4):
                    modified = header.modified(max_order=new_min_order,
                                               min_order=new_min_order)
                    self.assertIsInstance(modified, ClustHeader)
                    self.assertEqual(modified,
                                     make_header(max_order=new_min_order,
                                                 min_order=new_min_order))


class TestRelClustHeader(ut.TestCase):

    def test_clustered(self):
        self.assertTrue(RelClustHeader.clustered())

    def test_levels(self):
        self.assertEqual(RelClustHeader.levels(),
                         dict(rel=REL_NAME, order=ORDER_NAME, clust=CLUST_NAME))

    def test_num_levels(self):
        self.assertEqual(RelClustHeader.num_levels(), 3)

    def test_level_keys(self):
        self.assertEqual(RelClustHeader.level_keys(), ["rel", "order", "clust"])

    def test_level_names(self):
        self.assertEqual(RelClustHeader.level_names(),
                         [REL_NAME, ORDER_NAME, CLUST_NAME])

    def test_signature(self):
        for max_order in range(1, 11):
            for min_order in range(1, max_order + 1):
                sig = RelClustHeader(rels=list("qwerty"),
                                     max_order=max_order,
                                     min_order=min_order).signature
                self.assertEqual(list(sig.keys()),
                                 ["max_order", "min_order", "rels"])
                sig_max = sig["max_order"]
                self.assertIsInstance(sig_max, int)
                self.assertEqual(sig_max, max_order)
                sig_min = sig["min_order"]
                self.assertIsInstance(sig_min, int)
                self.assertEqual(sig_min, min_order)
                sig_rels = sig["rels"]
                self.assertIsInstance(sig_rels, np.ndarray)
                self.assertEqual(sig_rels.tolist(),
                                 ["q", "w", "e", "r", "t", "y"])

    def test_index(self):
        index = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2).index
        self.assertIsInstance(index, pd.MultiIndex)
        self.assertEqual(list(index.names), [REL_NAME, ORDER_NAME, CLUST_NAME])
        self.assertTrue(np.array_equal(index.get_level_values(REL_NAME),
                                       list("aaaaabbbbb")))
        self.assertTrue(np.array_equal(index.get_level_values(ORDER_NAME),
                                       [2, 2, 3, 3, 3, 2, 2, 3, 3, 3]))
        self.assertTrue(np.array_equal(index.get_level_values(CLUST_NAME),
                                       [1, 2, 1, 2, 3, 1, 2, 1, 2, 3]))

    def test_iter_clust_indexes(self):
        for max_order in range(1, 6):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                clust_indexes = list(header.iter_clust_indexes())
                self.assertEqual(len(clust_indexes), header.clusts.size)
                for index, clust in zip(clust_indexes,
                                        header.clusts,
                                        strict=True):
                    self.assertIsInstance(index, pd.MultiIndex)
                    self.assertEqual(index.size, len("qwerty"))
                    self.assertEqual(index.to_list(),
                                     [(rel, *clust) for rel in "qwerty"])

    def test_select_none(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select()
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertTrue(selection.equals(header.index))

    def test_select_rel(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(rel="b")
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [("b", 2, 1),
                          ("b", 2, 2),
                          ("b", 3, 1),
                          ("b", 3, 2),
                          ("b", 3, 3)])

    def test_select_order(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(order=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [("a", 3, 1),
                          ("a", 3, 2),
                          ("a", 3, 3),
                          ("b", 3, 1),
                          ("b", 3, 2),
                          ("b", 3, 3)])

    def test_select_clust(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(clust=2)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [("a", 2, 2),
                          ("a", 3, 2),
                          ("b", 2, 2),
                          ("b", 3, 2)])

    def test_select_order_clust_exist(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(rel="a", order=2, clust=1)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [("a", 2, 1)])

    def test_select_order_clust_empty(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(rel="a", order=2, clust=3)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [])

    def test_select_invalid_rel(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        self.assertRaisesRegex(ValueError,
                               "Expected rel to be one of",
                               header.select,
                               rel='c')

    def test_select_invalid_order(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        self.assertRaisesRegex(ValueError,
                               "Expected order to be one of",
                               header.select,
                               order=1)

    def test_select_invalid_clust(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        self.assertRaisesRegex(ValueError,
                               "Expected clust to be one of",
                               header.select,
                               clust=4)

    def test_select_extra(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        self.assertRaisesRegex(TypeError,
                               "Unexpected keyword arguments",
                               header.select,
                               extra="x")

    def test_select_extra_none(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(rel="a", order=2, clust=1, extra=None)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [("a", 2, 1)])

    def test_select_extra_zero(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(rel="a", order=2, clust=1, extra=0)
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [("a", 2, 1)])

    def test_select_extra_emptystr(self):
        header = RelClustHeader(rels=["a", "b"], max_order=3, min_order=2)
        selection = header.select(rel="a", order=2, clust=1, extra="")
        self.assertIsInstance(selection, pd.MultiIndex)
        self.assertEqual(selection.to_list(),
                         [("a", 2, 1)])

    def test_modified_none(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                self.assertEqual(header.modified(), header)

    def test_modified_rels(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                modified = header.modified(rels=list("uiop"))
                self.assertIsInstance(modified, RelClustHeader)
                self.assertEqual(modified,
                                 make_header(rels=list("uiop"),
                                             max_order=max_order,
                                             min_order=min_order))

    def test_modified_empty_rels(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                modified = header.modified(rels=[])
                self.assertIsInstance(modified, ClustHeader)
                self.assertEqual(modified,
                                 make_header(max_order=max_order,
                                             min_order=min_order))

    def test_modified_max_order(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                for new_max_order in range(1, 4):
                    modified = header.modified(max_order=new_max_order)
                    self.assertIsInstance(modified, RelClustHeader)
                    self.assertEqual(modified,
                                     make_header(rels=list("qwerty"),
                                                 max_order=new_max_order,
                                                 min_order=min_order))

    def test_modified_max_order_0(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                modified = header.modified(max_order=0)
                self.assertIsInstance(modified, RelHeader)
                self.assertEqual(modified,
                                 make_header(rels=list("qwerty"),
                                             max_order=0,
                                             min_order=min_order))

    def test_modified_all(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                for new_min_order in range(1, 4):
                    modified = header.modified(rels=list("uiop"),
                                               max_order=new_min_order,
                                               min_order=new_min_order)
                    self.assertIsInstance(modified, ClustHeader)
                    self.assertEqual(modified,
                                     make_header(rels=list("uiop"),
                                                 max_order=new_min_order,
                                                 min_order=new_min_order))

    def test_modified_nullified(self):
        for max_order in range(1, 4):
            for min_order in range(1, max_order + 1):
                header = RelClustHeader(rels=list("qwerty"),
                                        max_order=max_order,
                                        min_order=min_order)
                self.assertRaisesRegex(TypeError,
                                       "No header for rels",
                                       header.modified,
                                       rels=[],
                                       max_order=0)


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
        for omax1, omax2 in product(range(1, 4), repeat=2):
            for omin1, omin2 in product(range(1, 4), repeat=2):
                header1 = ClustHeader(max_order=omax1, min_order=omin1)
                header2 = ClustHeader(max_order=omax2, min_order=omin2)
                if (omax1, omin1) == (omax2, omin2):
                    self.assertEqual(header1, header2)
                else:
                    self.assertNotEqual(header1, header2)

    def test_relclustheaders(self):
        for rels1, rels2 in product(["qwerty", "uiop", "asdf", "ghjkl"],
                                    repeat=2):
            for omax1, omax2 in product(range(1, 4), repeat=2):
                for omin1, omin2 in product(range(1, 4), repeat=2):
                    header1 = RelClustHeader(rels=list(rels1),
                                             max_order=omax1,
                                             min_order=omin1)
                    header2 = RelClustHeader(rels=list(rels2),
                                             max_order=omax2,
                                             min_order=omin2)
                    if (rels1, omax1, omin1) == (rels2, omax2, omin2):
                        self.assertEqual(header1, header2)
                    else:
                        self.assertNotEqual(header1, header2)

    def test_different_types(self):
        for rels in ["qwerty", "uiop", "asdf", "ghjkl"]:
            for omax in range(1, 4):
                for omin in range(1, 4):
                    rh = RelHeader(rels=list(rels))
                    ch = ClustHeader(max_order=omax,
                                     min_order=omin)
                    rch = RelClustHeader(rels=list(rels),
                                         max_order=omax,
                                         min_order=omin)
                    for header1, header2 in combinations([rh, ch, rch], 2):
                        self.assertNotEqual(header1, header2)


class TestMakeHeader(ut.TestCase):

    def test_none(self):
        self.assertRaisesRegex(TypeError,
                               r"No header for rels=\[\] and max_order=0",
                               make_header)

    def test_rels(self):
        header = make_header(rels=["a", "b"])
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.tolist(), ["a", "b"])

    def test_max_order(self):
        header = make_header(max_order=2)
        self.assertIsInstance(header, ClustHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.tolist(),
                         [(1, 1), (2, 1), (2, 2)])

    def test_min_order(self):
        header = make_header(max_order=3, min_order=2)
        self.assertIsInstance(header, ClustHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.tolist(),
                         [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3)])

    def test_all(self):
        header = make_header(rels=["a", "b"], max_order=3, min_order=2)
        self.assertIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.tolist(),
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

    def test_none(self):
        self.assertRaisesRegex(TypeError,
                               r"No header for rels=\[\] and max_order=0",
                               parse_header,
                               pd.Index([]))

    def test_rel_index(self):
        header = parse_header(pd.Index(["a", "b"]))
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.to_list(), ["a", "b"])

    def test_rel_index_valid_name(self):
        header = parse_header(pd.Index(["a", "b"], name=REL_NAME))
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.to_list(), ["a", "b"])

    def test_rel_index_invalid_name(self):
        self.assertRaisesRegex(ValueError,
                               f"Expected index named {repr(REL_NAME)}",
                               parse_header,
                               pd.Index(["a", "b"], name=ORDER_NAME))

    def test_rel_multiindex(self):
        header = parse_header(pd.MultiIndex.from_arrays([["a", "b"]],
                                                        names=[REL_NAME]))
        self.assertIsInstance(header, RelHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.to_list(), ["a", "b"])

    def test_clust(self):
        header = parse_header(pd.MultiIndex.from_tuples([("1", "1"),
                                                         ("2", "1"),
                                                         ("2", "2")],
                                                        names=[ORDER_NAME,
                                                               CLUST_NAME]))
        self.assertIsInstance(header, ClustHeader)
        self.assertNotIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.to_list(),
                         [(1, 1), (2, 1), (2, 2)])

    def test_relclust(self):
        header = parse_header(pd.MultiIndex.from_tuples([("a", "1", "1"),
                                                         ("a", "2", "1"),
                                                         ("a", "2", "2"),
                                                         ("b", "1", "1"),
                                                         ("b", "2", "1"),
                                                         ("b", "2", "2")],
                                                        names=[REL_NAME,
                                                               ORDER_NAME,
                                                               CLUST_NAME]))
        self.assertIsInstance(header, RelClustHeader)
        self.assertEqual(header.index.to_list(),
                         [("a", 1, 1),
                          ("a", 2, 1),
                          ("a", 2, 2),
                          ("b", 1, 1),
                          ("b", 2, 1),
                          ("b", 2, 2)])

    def test_rel_index_repeated(self):
        self.assertRaisesRegex(ValueError,
                               "Indexes do not match",
                               parse_header,
                               pd.Index(["a", "b", "a"]))

    def test_missing_index_names(self):
        self.assertRaisesRegex(TypeError,
                               r"No header for rels=\[\] and max_order=0",
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
                                                         names=[ORDER_NAME,
                                                                CLUST_NAME]))

    def test_extra_values(self):
        self.assertRaisesRegex(ValueError,
                               "Invalid index level",
                               parse_header,
                               pd.MultiIndex.from_tuples([("1", "1"),
                                                          ("2", "1"),
                                                          ("2", "1"),
                                                          ("2", "2")],
                                                         names=[ORDER_NAME,
                                                                CLUST_NAME]))

    def test_nonnumeric(self):
        self.assertRaisesRegex(ValueError,
                               "Expected index names",
                               parse_header,
                               pd.MultiIndex.from_tuples([("a", "x", "x"),
                                                          ("a", "y", "x"),
                                                          ("a", "y", "y"),
                                                          ("b", "x", "x"),
                                                          ("b", "y", "x"),
                                                          ("b", "y", "y")],
                                                         names=[REL_NAME,
                                                                ORDER_NAME,
                                                                CLUST_NAME]))

    def test_invalid_numeric(self):
        self.assertRaisesRegex(ValueError,
                               "Expected index names",
                               parse_header,
                               pd.MultiIndex.from_tuples([("a", '0', '0'),
                                                          ("b", '0', '0')],
                                                         names=[REL_NAME,
                                                                ORDER_NAME,
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
