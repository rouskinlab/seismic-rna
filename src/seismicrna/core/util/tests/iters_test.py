import unittest as ut

from seismicrna.core.util.iters import filterboth, partition


class TestFilterboth(ut.TestCase):

    def test_none(self):
        items = [False, True, True, False, False]
        trues, falses = map(list, filterboth(None, items))
        self.assertEqual(trues, [True, True])
        self.assertEqual(falses, [False, False, False])

    def test_callable(self):
        items = [8, 6, 7, 5, 3, 0, 9]
        even, odd = map(list, filterboth(lambda i: i % 2 == 0, items))
        self.assertEqual(even, [8, 6, 0])
        self.assertEqual(odd, [7, 5, 3, 9])


class TestPartition(ut.TestCase):

    def test_none(self):
        items = [False, True, True, False, False]
        trues, falses = partition(None, items)
        self.assertEqual(trues, [True, True])
        self.assertEqual(falses, [False, False, False])

    def test_callable(self):
        items = [8, 6, 7, 5, 3, 0, 9]
        even, odd = partition(lambda i: i % 2 == 0, items)
        self.assertEqual(even, [8, 6, 0])
        self.assertEqual(odd, [7, 5, 3, 9])


if __name__ == "__main__":
    ut.main()
