import unittest as ut
from types import SimpleNamespace

from seismicrna.cluster.write import choose_best_k
from seismicrna.core.error import NoDataError
from seismicrna.core.logs import Level, restore_config, set_config


def make_runs_k(k: int, passing: bool, bic: float):
    """Fake the EM runs for one number of clusters, with only the
    attributes that find_best_k reads."""
    return SimpleNamespace(
        k=k, passing=lambda **kwargs: passing, best=SimpleNamespace(bic=bic)
    )


class TestChooseBestK(ut.TestCase):
    """The best number of clusters is the one with the lowest BIC among
    those that passed the filters, of any K. If none passed, the ensemble
    average (K = 1) is used, but only if it was run: it alone can fail
    only as underclustered, so it can never contribute a cluster that the
    filters rejected."""

    @restore_config
    def setUp(self):
        # The fallback logs a warning, which is not the subject here.
        set_config(verbosity=Level.ERROR)

    def test_best_of_those_passing(self):
        """A K greater than 1 is chosen whenever it passed and has the
        lowest BIC."""
        runs_ks = {
            1: make_runs_k(1, True, 100.0),
            2: make_runs_k(2, True, 50.0),
            3: make_runs_k(3, True, 75.0),
        }
        self.assertEqual(choose_best_k(runs_ks), 2)

    def test_skips_those_not_passing(self):
        """A K with a better BIC is skipped if it failed the filters."""
        runs_ks = {1: make_runs_k(1, True, 100.0), 2: make_runs_k(2, False, 10.0)}
        self.assertEqual(choose_best_k(runs_ks), 1)

    def test_only_a_high_k_passing(self):
        """If only a K greater than 1 passed, it is chosen even though
        the ensemble average did not pass."""
        runs_ks = {1: make_runs_k(1, False, 10.0), 2: make_runs_k(2, True, 90.0)}
        self.assertEqual(choose_best_k(runs_ks), 2)

    def test_none_passing_falls_back_to_one(self):
        """With no K passing, the ensemble average is used, never a K
        greater than 1 that failed and so could be overclustered."""
        runs_ks = {
            1: make_runs_k(1, False, 100.0),
            2: make_runs_k(2, False, 50.0),
            3: make_runs_k(3, False, 10.0),
        }
        self.assertEqual(choose_best_k(runs_ks), 1)

    def test_none_passing_without_one_writes_nothing(self):
        """With no K passing and K = 1 never run, there is nothing safe to
        fall back to, so nothing is written."""
        runs_ks = {2: make_runs_k(2, False, 100.0), 3: make_runs_k(3, False, 50.0)}
        self.assertRaises(NoDataError, choose_best_k, runs_ks, 2)

    def test_none_passing_names_min_clusters(self):
        """The error names the option that ruled out the fallback."""
        runs_ks = {3: make_runs_k(3, False, 100.0)}
        self.assertRaisesRegex(
            NoDataError, r"--min-clusters 3", choose_best_k, runs_ks, 3
        )

    def test_no_runs_is_an_error(self):
        self.assertRaisesRegex(ValueError, "from no runs", choose_best_k, dict())


if __name__ == "__main__":
    ut.main(verbosity=2)
