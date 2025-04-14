import os
from abc import ABC, abstractmethod
from concurrent import futures
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import beta as beta_dist, norm, gamma as gamma_dist
from scipy.optimize import minimize

from seismicrna.core.rna import RNAProfile, RNAStructure, RNAState
from seismicrna.mask.table import MaskPositionTableLoader

PARAMS_ORDER = [f"{pu}_{ab}"
                for pu in ["paired", "unpaired"]
                for ab in ["alpha", "beta"]]
MAIN_DIR = Path(__file__).resolve().parent.parent
MASK_DIR = MAIN_DIR.joinpath("dms", "out-astrophage", "pooled", "mask")
STRUCTS_FILE = Path("structs.csv")
MUS_FILE_095 = Path("mus_auc-0.95.csv")
MUS_FILE_090 = Path("mus_auc-0.90.csv")
MUTATION_CDF_PLOT_FILE = Path("mutation-cdfs.pdf")
MUTATION_COUNTS_PLOT_FILE = Path("mutation-counts.pdf")
CROSS_VALIDATE_FILE = Path("cross-validation.csv")
CROSS_VALIDATE_PLOT_FILE = Path("cross-validation.pdf")
CROSS_VALIDATE_PLOT_DIR = Path("cross-validation-plots")
MUTATION_HISTOGRAM_FILE = Path("mutation-histograms.pdf")
BOOTSTRAP_HISTOGRAMS_N2_N2_K2_FILE = Path("bootstrap-histograms_n2-n2-k2.pdf")
BOOTSTRAP_HISTOGRAMS_K1_N1_K2_FILE = Path("bootstrap-histograms_k1-n1-k2.pdf")
BOOTSTRAP_HISTOGRAMS_G1_G1_G2_FILE = Path("bootstrap-histograms_g1-g1-g2.pdf")
BOOTSTRAP_PARAMS_N2_N2_K2_FILE = Path("bootstrap-params_n2-n2-k2.csv")
BOOTSTRAP_PARAMS_K1_N1_K2_FILE = Path("bootstrap-params_k1-n1-k2.csv")
BOOTSTRAP_PARAMS_G1_G1_G2_FILE = Path("bootstrap-params_g1-g1-g2.csv")
BOOTSTRAP_PSEUDOENERGY_K1_N1_K2_FILE = Path("pseudoenergies_k1-n1-k2.pdf")

MAX_MU = 0.5
EPS = 1.e-6

rng = np.random.default_rng()


def rsf(std: float = 1.0):
    """ Random scale factor. """
    return np.exp(rng.normal(0, std))


def get_archiveii():
    return pd.read_csv(MAIN_DIR.joinpath("structures", "archiveii.csv"),
                       index_col="id")


def read_mask_position_table_file(ref: str):
    table_file = MASK_DIR.joinpath(ref, "full", "mask-position-table.csv")
    return MaskPositionTableLoader(table_file)


def get_sequence(table: MaskPositionTableLoader):
    return str(table.refseq.tr())


def map_refs_to_structs():
    """ Map the reference names to their structures. """
    refs_to_structs = dict()
    archiveii = get_archiveii()
    for ref_dir in MASK_DIR.iterdir():
        table_ref = ref_dir.name
        if table_ref in refs_to_structs:
            raise ValueError(table_ref)
        # Find the secondary structures of the entries in ArchiveII with
        # the same reference sequence.
        refseq = get_sequence(read_mask_position_table_file(table_ref))
        structs = archiveii.loc[archiveii["sequence"] == refseq,
                                "secondary_structure"]
        n_structs = len(set(structs))
        if n_structs != 1:
            print(f"Skipped {repr(table_ref)} with {n_structs} structures")
            continue
        struct = structs.iloc[0]
        refs_to_structs[table_ref] = struct
    refs_to_structs = pd.Series(refs_to_structs)
    refs_to_structs.name = "secondary_structure"
    refs_to_structs.index.name = "reference"
    return refs_to_structs.to_frame()


def get_mutation_rates(table: MaskPositionTableLoader):
    return table.fetch_ratio(rel="Mutated", squeeze=True).dropna()


def gather_mus(refs_to_structs: pd.DataFrame, min_auc=0.95):
    mus = list()
    for ref, db_string in refs_to_structs["secondary_structure"].items():
        table = read_mask_position_table_file(ref)
        refseq = get_sequence(table)
        struct = RNAStructure.from_db_string(db_string,
                                             refseq,
                                             ref=ref,
                                             reg="full",
                                             title=ref,
                                             branch="")
        is_internal = struct.is_paired_internally
        n_internal = np.count_nonzero(is_internal)
        is_terminal = struct.is_paired_terminally
        n_terminal = np.count_nonzero(is_terminal)
        is_unpaired = struct.is_unpaired
        n_unpaired = np.count_nonzero(is_unpaired)
        n_total = n_internal + n_terminal + n_unpaired
        if n_total > 0:
            f_paired = (n_internal + n_terminal) / n_total
        else:
            f_paired = np.nan
        # Use not >= (instead of <) to remove structures for which
        # fraction_paired is NaN.
        if not f_paired >= 1/3:
            print(f"Skipped {repr(ref)} with {f_paired * 100} % "
                  f"paired: {db_string}")
            continue
        mus_table = get_mutation_rates(table)
        state = RNAState.from_struct_profile(struct,
                                             RNAProfile(sample="",
                                                        branches={},
                                                        data_reg="full",
                                                        data_name="average",
                                                        data=mus_table,
                                                        region=struct.region))
        auc = state.calc_auc(terminal_pairs=False)
        # Use not >= (instead of <) to remove structures for which
        # AUC is NaN.
        if not auc >= min_auc:
            print(f"Skipped {repr(ref)} with AUC-ROC = {auc}")
            continue
        for base in ["A", "C"]:
            mus_base = mus_table.loc[mus_table.index.get_level_values("Base") == base]
            mus_internal = mus_base.loc[is_internal.loc[mus_base.index]]
            mus_terminal = mus_base.loc[is_terminal.loc[mus_base.index]]
            mus_unpaired = mus_base.loc[is_unpaired.loc[mus_base.index]]
            for pairing, mus_status in [("Internal", mus_internal),
                                        ("Terminal", mus_terminal),
                                        ("Unpaired", mus_unpaired)]:
                for (pos, base_), mu in mus_status.items():
                    assert base_ == base
                    mus.append({"Ref": ref,
                                "Position": pos,
                                "Base": base,
                                "Pairing": pairing,
                                "Mutation": mu})
    return pd.DataFrame.from_records(mus)


class FittedDistribution(ABC):

    @classmethod
    def get_name(cls):
        assert cls.__name__.endswith("Distribution")
        return cls.__name__[:-len("Distribution")]

    @classmethod
    @abstractmethod
    def get_bounds(cls) -> list[tuple[float | None, float | None]]:
        pass

    @classmethod
    def get_num_params(cls) -> int:
        return len(cls.get_bounds())

    @classmethod
    @abstractmethod
    def _calc_logpdf(cls, x: np.ndarray, *params: float) -> np.ndarray:
        """ Calculate the log probability density function. """

    @classmethod
    @abstractmethod
    def _initial_guess(cls, x: np.ndarray) -> tuple[float, ...]:
        """ Get an initial guess for the distribution parameters. """

    @classmethod
    def _fit_params(cls,
                    x: np.ndarray,
                    n_stderr: int = 9,
                    max_stderr: float = 0.1,
                    max_consecutive_failures: int = 300,
                    max_iter_per_try: int = 1000) -> tuple[float, ...]:
        """ Fit the distribution parameters to the data. """
        print(f"Fitting {cls.__name__} to {len(x)} data points")
        best_result = None
        loglikes = list()
        stderr = np.inf
        tries = 0
        consecutive_failures = 0
        while consecutive_failures < max_consecutive_failures and (len(loglikes) < n_stderr or stderr > max_stderr):
            tries += 1
            result = minimize(
                lambda params: -np.sum(cls._calc_logpdf(x, *params)),
                np.asarray(cls._initial_guess(x)),
                method='L-BFGS-B',
                bounds=cls.get_bounds(),
                options={"maxiter": max_iter_per_try}
            )
            if result.success:
                if best_result is None or result.fun < best_result.fun:
                    best_result = result
                    consecutive_failures = 0
                else:
                    consecutive_failures += 1
                loglikes.append(-result.fun)
                if len(loglikes) >= n_stderr:
                    loglikes.sort(reverse=True)
                    stderr = np.std(loglikes[:n_stderr]) / np.sqrt(n_stderr)
            else:
                consecutive_failures += 1
            print(f"\t{tries}\t{consecutive_failures}\t{-result.fun}\t{-best_result.fun if best_result is not None else np.nan}\t{stderr}")
        if best_result is None:
            raise ValueError(f"Failed to fit parameters for {cls.get_name()}")
        return tuple(map(float, best_result.x))

    def __init__(self, x: Any, eps: float = EPS, **kwargs):
        if not 0.0 <= eps <= 0.5:
            raise ValueError(f"Invalid eps: {eps}")
        self.x = np.clip(np.asarray_chkfinite(x), eps, 1. - eps)
        self.params = tuple(map(float, self._fit_params(self.x, **kwargs)))

    def calc_logpdf(self, x: np.ndarray | None = None) -> np.ndarray:
        """ Log probability density function. """
        if x is not None:
            x = np.clip(np.asarray_chkfinite(x), EPS, 1. - EPS)
        else:
            x = self.x
        return self._calc_logpdf(x, *self.params)

    def calc_loglike(self, x: np.ndarray | None = None) -> float:
        """ Log likelihood function. """
        return np.sum(self.calc_logpdf(x))

    def calc_pdf(self, x: np.ndarray | None = None) -> np.ndarray:
        """ Probability density function. """
        return np.exp(self.calc_logpdf(x))

    def calc_aic(self, x: np.ndarray | None = None) -> float:
        """ Akaike Information Criterion. """
        return 2 * self.get_num_params() - 2 * self.calc_loglike(x)


class LogitNormalDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(-1. / EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, mean: float, std: float):
        log_x = np.log(x)
        log_1_minus_x = np.log(1. - x)
        return norm.logpdf(log_x - log_1_minus_x, mean, std) - (log_x + log_1_minus_x)

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        raise NotImplementedError

    @classmethod
    def _fit_params(cls, x: np.ndarray):
        return norm.fit(np.log(x / (1. - x)))


class BetaDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, a: float, b: float):
        return beta_dist.logpdf(x, a, b)

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        raise NotImplementedError

    @classmethod
    def _fit_params(cls, x: np.ndarray):
        return beta_dist.fit(x, floc=0, fscale=1)[:2]


class GammaDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, a: float, scale: float):
        return gamma_dist.logpdf(x, a, scale=scale)

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        raise NotImplementedError

    @classmethod
    def _fit_params(cls, x: np.ndarray):
        a, _, scale = gamma_dist.fit(x, floc=0)
        return a, scale


class KumaraswamyDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, a: float, b: float):
        return (np.log(a) + np.log(b)) + (a-1)*np.log(x) + (b-1)*np.log(1-x**a)

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        mean = np.mean(x)
        if mean <= 1. - EPS:
            a_init = 1. / (1. - mean)
        else:
            a_init = cls.get_bounds()[0][1]
        if mean >= EPS:
            b_init = 1. / mean
        else:
            b_init = cls.get_bounds()[1][1]
        a_init = np.clip(rsf() * a_init, cls.get_bounds()[0][0], cls.get_bounds()[0][1])
        b_init = np.clip(rsf() * b_init, cls.get_bounds()[1][0], cls.get_bounds()[1][1])
        return a_init, b_init


class DoubleKumaraswamyDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(EPS, 1. - EPS), (EPS, 1. / EPS), (EPS, 1. / EPS), (EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, w: float, a1: float, b1: float, a2: float, b2: float):
        logpdf1 = KumaraswamyDistribution._calc_logpdf(x, a1, b1)
        logpdf2 = KumaraswamyDistribution._calc_logpdf(x, a2, b2)
        return np.log(np.maximum(w * np.exp(logpdf1) + (1-w) * np.exp(logpdf2), EPS))

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        mean = np.mean(x)
        if mean <= 1. - EPS:
            a_init = 1. / (1. - mean)
        else:
            a_init = cls.get_bounds()[0][1]
        if mean >= EPS:
            b_init = 1. / mean
        else:
            b_init = cls.get_bounds()[1][1]
        w_init = rng.uniform(*cls.get_bounds()[0])
        a1_init = np.clip(rsf() * a_init, *cls.get_bounds()[1])
        b1_init = np.clip(rsf() * b_init, *cls.get_bounds()[2])
        a2_init = np.clip(rsf() * a_init, *cls.get_bounds()[3])
        b2_init = np.clip(rsf() * b_init, *cls.get_bounds()[4])
        return w_init, a1_init, b1_init, a2_init, b2_init


class DoubleLogitNormalDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(EPS, 1. - EPS), (-1. / EPS, 1. / EPS), (EPS, 1. / EPS), (-1. / EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, w: float, a1: float, b1: float, a2: float, b2: float):
        logpdf1 = LogitNormalDistribution._calc_logpdf(x, a1, b1)
        logpdf2 = LogitNormalDistribution._calc_logpdf(x, a2, b2)
        return np.log(np.maximum(w * np.exp(logpdf1) + (1-w) * np.exp(logpdf2), EPS))

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        log_x = np.log(x)
        mean = np.mean(log_x)
        std = np.std(log_x)
        w_init = rng.uniform(*cls.get_bounds()[0])
        mean1_init = np.clip(rsf() * mean, *cls.get_bounds()[1])
        std1_init = np.clip(rsf() * std, *cls.get_bounds()[2])
        mean2_init = np.clip(rsf() * mean, *cls.get_bounds()[3])
        std2_init = np.clip(rsf() * std, *cls.get_bounds()[4])
        return w_init, mean1_init, std1_init, mean2_init, std2_init


class DoubleBetaDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(EPS, 1. - EPS), (EPS, 1. / EPS), (EPS, 1. / EPS), (EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, w: float, a1: float, b1: float, a2: float, b2: float):
        pdf1 = beta_dist.pdf(x, a1, b1)
        pdf2 = beta_dist.pdf(x, a2, b2)
        return np.log(np.maximum(w * pdf1 + (1-w) * pdf2, EPS))

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        a, b = beta_dist.fit(x, floc=0, fscale=1)[:2]
        w_init = rng.uniform(*cls.get_bounds()[0])
        a1_init = np.clip(rsf() * a, *cls.get_bounds()[1])
        b1_init = np.clip(rsf() * b, *cls.get_bounds()[2])
        a2_init = np.clip(rsf() * a, *cls.get_bounds()[3])
        b2_init = np.clip(rsf() * b, *cls.get_bounds()[4])
        return w_init, a1_init, b1_init, a2_init, b2_init


class DoubleGammaDistribution(FittedDistribution):

    @classmethod
    def get_bounds(cls):
        return [(EPS, 1. - EPS), (EPS, 1. / EPS), (EPS, 1. / EPS), (EPS, 1. / EPS), (EPS, 1. / EPS)]

    @classmethod
    def _calc_logpdf(cls, x: np.ndarray, w: float, a1: float, scale1: float, a2: float, scale2: float):
        pdf1 = gamma_dist.pdf(x, a1, scale=scale1)
        pdf2 = gamma_dist.pdf(x, a2, scale=scale2)
        return np.log(np.maximum(w * pdf1 + (1-w) * pdf2, EPS))

    @classmethod
    def _initial_guess(cls, x: np.ndarray):
        a, _, scale = gamma_dist.fit(x, floc=0)
        w_init = rng.uniform(*cls.get_bounds()[0])
        a1_init = np.clip(rsf() * a, *cls.get_bounds()[1])
        scale1_init = np.clip(rsf() * scale, *cls.get_bounds()[2])
        a2_init = np.clip(rsf() * a, *cls.get_bounds()[3])
        scale2_init = np.clip(rsf() * scale, *cls.get_bounds()[4])
        return w_init, a1_init, scale1_init, a2_init, scale2_init


DISTS = {KumaraswamyDistribution: "b:",
         BetaDistribution: "r:",
         LogitNormalDistribution: "g:",
         DoubleKumaraswamyDistribution: "b--",
         DoubleBetaDistribution: "r--",
         DoubleLogitNormalDistribution: "g--",
         DoubleGammaDistribution: "k-"}


# Graphing functions


def plot_counts(mus: pd.DataFrame):
    """ Plot the counts of the mutation rates. """
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    fig, ax = plt.subplots(figsize=(width_in, height_in))
    # Get RNA types from reference names
    rna_types = mus["Ref"].str.split("_").str[0]
    unique_rna_types = rna_types.unique()

    # Set up the plot
    bases = ["A", "C"]
    pairings = ["Internal", "Terminal", "Unpaired"]
    x_pos = np.arange(len(bases) * len(pairings))
    width = 0.8

    # Set up colors for RNA types
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_rna_types)))

    # Initialize bottom of stacked bars
    bottom = np.zeros(len(x_pos))

    # Create stacked bars for each RNA type
    for i, rna_type in enumerate(unique_rna_types):
        counts = []
        for base in bases:
            for pairing in pairings:
                mask = ((mus["Base"] == base) &
                       (mus["Pairing"] == pairing) &
                       (rna_types == rna_type))
                counts.append(np.sum(mask))

        ax.bar(x_pos, counts, width, bottom=bottom,
               label=rna_type, color=colors[i])
        bottom += counts

    # Customize plot
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"{b} - {p}" for b in bases for p in pairings],
                       rotation=0)
    ax.set_ylabel("Number of Mutation Rate Measurements")
    ax.set_xlabel("Base - Pairing")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title="RNA Type")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(MUTATION_COUNTS_PLOT_FILE)
    plt.close()


def plot_mutation_cdfs(mus: pd.DataFrame):
    """ Plot the cumulative distribution functions of the mutation rates. """
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    rna_types = mus["Ref"].str.split("_").str[0]
    unique_rna_types = rna_types.unique()
    fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(width_in, height_in))
    for i, base in enumerate(["A", "C"]):
        for j, pairing in enumerate(["Internal", "Terminal", "Unpaired"]):
            print(f"CDF: {base} - {pairing}")
            mask_ij = (mus["Base"] == base) & (mus["Pairing"] == pairing)
            x_ij = np.sort(np.clip(mus.loc[mask_ij, "Mutation"].values, EPS, 1. - EPS))
            x_ij = np.concatenate([[0.], x_ij, [MAX_MU]])
            ax = axes[i, j]
            for rna_type in unique_rna_types:
                mask_rna = (rna_type == rna_types) & mask_ij
                mutation_rates = np.sort(np.clip(mus.loc[mask_rna, "Mutation"].values, EPS, 1. - EPS))
                # Calculate empirical CDF using all mutation rates as x coordinates
                cdf = np.searchsorted(mutation_rates, x_ij, side='right') / len(mutation_rates)
                # Plot the distribution as step function
                ax.step(x_ij, cdf, where='post', lw=1, label=rna_type)
                ax.set_title(f"{base} - {pairing}")
                ax.legend(loc='lower right')
                ax.set_xlabel("Mutation Rate")
                ax.set_ylabel("Cumulative Probability")
    plt.tight_layout()
    plt.savefig(MUTATION_CDF_PLOT_FILE)
    plt.close()


def cross_validate(mus: pd.DataFrame):
    """ Cross-validate the mutation rate distributions. """
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    results = list()
    rna_types = mus["Ref"].str.split("_").str[0]
    unique_rna_types = rna_types.unique()
    for base in ["A", "C"]:
        for pairing in ["Internal", "Terminal", "Unpaired"]:
            mask = (mus["Base"] == base) & (mus["Pairing"] == pairing)
            fig, axes = plt.subplots(nrows=len(unique_rna_types), ncols=1, figsize=(width_in, height_in))
            for i, rna_type in enumerate(unique_rna_types):
                ax = axes[i]
                mask_train = mask & (rna_types != rna_type)
                mask_test = mask & (rna_types == rna_type)
                mutation_rates_train = mus.loc[mask_train, "Mutation"].values
                mutation_rates_test = mus.loc[mask_test, "Mutation"].values
                # Plot the histogram of the raw data.
                bins = np.linspace(0, MAX_MU, 101)
                ax.hist(mutation_rates_train, bins=bins, density=True, alpha=1/3, label="Train Data")
                ax.hist(mutation_rates_test, bins=bins, density=True, alpha=1/3, label="Test Data")
                for dist_type, style in DISTS.items():
                    # Fit this type of distribution.
                    dist = dist_type(mutation_rates_train)
                    dist_name = dist_type.get_name()
                    loglike_train = dist.calc_loglike(mutation_rates_train) / len(mutation_rates_train)
                    loglike_test = dist.calc_loglike(mutation_rates_test) / len(mutation_rates_test)
                    results.append({"base": base,
                                    "pairing": pairing,
                                    "dist": dist_name,
                                    "rna_type": rna_type,
                                    "loglike_train": loglike_train,
                                    "loglike_test": loglike_test})
                    # Plot the fitted distribution.
                    x = np.linspace(0.005, MAX_MU - 0.005, 491)
                    ax.plot(x, dist.calc_pdf(x), style, label=f"{dist_name} (LL Train = {loglike_train:.3f}; LL Test = {loglike_test:.3f})")

                title = f"Test Group: {rna_type} (N Train = {len(mutation_rates_train)}; N Test = {len(mutation_rates_test)})"
                print(title)
                print()
                ax.set_title(title)
                ax.set_xlabel("Mutation Rate")
                ax.set_ylabel("Density")
                # Sort legend entries by loglike_test value
                handles, labels = ax.get_legend_handles_labels()
                # Extract loglike values from labels that contain them
                loglikes = []
                other_items = []
                other_handles = []
                for h, l in zip(handles, labels):
                    if "LL Test =" in l:
                        ll = float(l.split("LL Test = ")[-1].strip(")"))
                        loglikes.append((h, l, ll))
                    else:
                        other_items.append(l)
                        other_handles.append(h)
                # Sort distribution entries by loglike_test
                loglikes.sort(key=lambda x: x[2], reverse=True)
                sorted_handles = other_handles + [x[0] for x in loglikes]
                sorted_labels = other_items + [x[1] for x in loglikes]
                ax.legend(sorted_handles, sorted_labels, bbox_to_anchor=(1.05, 1), loc='upper left', ncol=1)
            plt.tight_layout()
            CROSS_VALIDATE_PLOT_DIR.mkdir(exist_ok=True)
            plt.savefig(CROSS_VALIDATE_PLOT_DIR.joinpath(f"{base}-{pairing}.pdf"))
            plt.close()
    return pd.DataFrame.from_records(results)


def plot_cross_validation(cv_results: pd.DataFrame):
    """Plot stacked bar charts of cross-validation results."""
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    fig, axes = plt.subplots(2, 3, figsize=(width_in, height_in))

    # Set up the plot grid
    bases = ["A", "C"]
    pairings = ["Internal", "Terminal", "Unpaired"]

    # Get unique distributions and RNA types
    rna_types = cv_results["rna_type"].unique()

    # Set up colors for RNA types
    colors = plt.cm.Set3(np.linspace(0, 1, len(rna_types)))

    # Create plots for each base and pairing combination
    for i, base in enumerate(bases):
        for j, pairing in enumerate(pairings):
            ax = axes[i, j]

            # Filter data for this combination
            mask = (cv_results["base"] == base) & (cv_results["pairing"] == pairing)
            subset = cv_results[mask]

            # Calculate total loglike_test per distribution and sort
            dist_totals = subset.groupby("dist")["loglike_test"].sum()
            sorted_distributions = dist_totals.sort_values(ascending=False).index.tolist()

            # Create stacked bars using the sorted order
            bottom = np.zeros(len(sorted_distributions))
            bars = []
            for k, rna_type in enumerate(rna_types):
                values = []
                for dist in sorted_distributions:  # Use sorted order
                    mask = (subset["dist"] == dist) & (subset["rna_type"] == rna_type)
                    val = subset.loc[mask, "loglike_test"].values
                    values.append(val[0] if len(val) > 0 else 0)

                bar = ax.bar(sorted_distributions, values, bottom=bottom,  # Use sorted order for x-axis
                           color=colors[k], label=rna_type)
                bottom += values
                bars.append(bar)

            # Customize subplot
            ax.set_title(f"{base} - {pairing}")
            ax.set_xticks(range(len(sorted_distributions))) # Set tick positions
            ax.set_xticklabels(sorted_distributions, ha='right', rotation_mode='anchor')
            ax.tick_params(axis='x', rotation=60)
            # Get handles and labels, reverse them, and set the legend
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.set_xlabel("Distribution")
            ax.set_ylabel("Mean Test Log-Likelihood per Data Point")

    # Adjust layout
    plt.tight_layout()
    plt.savefig(CROSS_VALIDATE_PLOT_FILE)
    plt.close(fig)


def plot_mutation_histograms(mus: pd.DataFrame):
    """Create histograms comparing mutation rates with various distribution fits."""
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    fig, axes = plt.subplots(3, 2, figsize=(width_in, height_in))

    # Create a grid of subplots for each combination
    for i, base in enumerate(["A", "C"]):
        for j, pairing in enumerate(["Internal", "Terminal", "Unpaired"]):
            ax = axes[j, i]

            # Get mutation rates for this combination
            mask = (mus["Base"] == base) & (mus["Pairing"] == pairing)
            mutation_rates = mus.loc[mask, "Mutation"].values
            n_mutations = len(mutation_rates)

            # Create histogram
            ax.hist(mutation_rates,
                    bins=np.linspace(0, MAX_MU, 101),
                    density=True,
                    alpha=0.6,
                    label="Observed",
                    color="skyblue")

            # Fit various distributions
            x = np.linspace(0.005, MAX_MU - 0.005, 491)
            distributions = []
            for dist_type, style in DISTS.items():
                dist = dist_type(mutation_rates)
                aic = dist.calc_aic()
                pdf = dist.calc_pdf(x)
                distributions.append((dist_type.get_name(), pdf, aic, style))

            # Sort distributions by AIC (smallest to largest)
            distributions.sort(key=lambda x: x[2])

            # Plot distributions in order of AIC
            for name, pdf, aic, style in distributions:
                ax.plot(x, pdf, style, lw=1, label=f"{name} (AIC={aic:.1f})")

            # Add count information to title
            ax.set_title(f"{base} - {pairing}\nN={n_mutations:,}")
            ax.set_xlabel("Mutation Rate")
            ax.set_ylabel("Density")
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(MUTATION_HISTOGRAM_FILE)
    plt.close()


def process_bootstrap(dist_type: type[FittedDistribution],
                      unique_rna_types: Iterable[str],
                      rna_type_indexes: dict[str, np.ndarray],
                      mutation_rates_ij: np.ndarray,
                      category_size: int,
                      seed: int):
    print(f"\tBootstrapping with seed {seed}")
    rng = np.random.default_rng(seed)
    rna_type_mutation_rates = []
    for rna_type in unique_rna_types:
        rna_type_mutation_rates.append(mutation_rates_ij[rng.choice(rna_type_indexes[rna_type],
                                                                    size=category_size,
                                                                    replace=True)])
    bootstrap_rates = np.concatenate(rna_type_mutation_rates)
    return dist_type(bootstrap_rates)


def bootstrap_mutation_histograms(mus: pd.DataFrame,
                                  plot_file: Path,
                                  dist_types: dict[str, type[FittedDistribution]],
                                  n_bootstraps: int = 120):
    """Bootstrap the mutation histograms."""
    if n_bootstraps < 2:
        raise ValueError("n_bootstraps must be at least 2")
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    fig, axes = plt.subplots(3, 2, figsize=(width_in, height_in))
    rna_types = mus["Ref"].str.split("_").str[0]
    unique_rna_types = rna_types.unique()
    x = np.linspace(0.005, MAX_MU - 0.005, 491)
    params = list()
    for i, base in enumerate(["A", "C"]):
        for j, pairing in enumerate(["Internal", "Terminal", "Unpaired"]):
            print(f"Bootstrapping {base} - {pairing}")
            # Get the mutation rates for this combination.
            mask_ij = (mus["Base"] == base) & (mus["Pairing"] == pairing)
            mutation_rates_ij = mus.loc[mask_ij, "Mutation"].values
            category_size = round(len(mutation_rates_ij) / len(unique_rna_types))
            rna_types_ij = rna_types[mask_ij]
            rna_type_indexes = {rna_type: np.flatnonzero(rna_types_ij == rna_type)
                                for rna_type in unique_rna_types}

            # Generate the bootstrap datasets and fit the distributions.
            with futures.ProcessPoolExecutor() as executor:
                futures_list = [executor.submit(process_bootstrap,
                                                dist_types[pairing],
                                                unique_rna_types,
                                                rna_type_indexes,
                                                mutation_rates_ij,
                                                category_size,
                                                rng.integers(1000000))
                                for _ in range(n_bootstraps)]
                dists = [future.result() for future in futures.as_completed(futures_list)]
            # Select the distribution with the best AIC on the other distributions' data.
            best_b = -1
            best_aic = np.inf
            for b, dist in enumerate(dists):
                other_data = np.concatenate([d.x for d in dists[:b] + dists[b+1:]])
                aic = dist.calc_aic(other_data)
                if aic < best_aic:
                    best_aic = aic
                    best_b = b
            if best_b == -1:
                raise ValueError("No distribution found")
            best_dist = dists[best_b]
            params.append({"base": base,
                           "pairing": pairing,
                           "dist_type": best_dist.get_name(),
                           "aic": best_aic}
                           | {f"param{i}": p for i, p in enumerate(best_dist.params)})

            # Plot the distributions.
            ax = axes[j, i]
            for dist in dists:
                ax.hist(dist.x, bins=np.linspace(0, MAX_MU, 101), density=True, color="skyblue", alpha=0.1)
                if dist != best_dist:
                    ax.plot(x, dist.calc_pdf(x), "k-", lw=1, alpha=0.1)
            ax.plot(x, best_dist.calc_pdf(x), "r-", lw=2, alpha=0.8)

            # Add count information to title
            ax.set_title(f"{base} - {pairing}\nN={category_size * len(unique_rna_types)}")
            ax.set_xlabel("Mutation Rate")
            ax.set_ylabel("Density")
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()
    return pd.DataFrame.from_records(params)


def calc_rnastructure_defaults_nt(x: np.ndarray):
    """Calculate the default nucleotide-specific pseudoenergy for RNAStructure."""
    R = 1.9872041e-3  # kcal/mol/K
    T = 310.15  # Kelvin
    dist_file = Path(os.environ["DATAPATH"]).joinpath("dists", "DMSdist_nt.txt")
    print(dist_file)
    params = dict()
    with open(dist_file) as f:
        f.readline()
        for base in ["A", "C", "G", "U"]:
            for pairing in ["Paired", "Unpaired"]:
                (shape1,
                 location1,
                 scale1,
                 shape2,
                 location2,
                 scale2,
                 weight1,
                 weight2) = map(float, f.readline().split())
                if abs(weight1 + weight2 - 1.0) > 1e-6:
                    raise ValueError(f"Weights for {base, pairing} do not sum to 1.0")
                if location1 != 0:
                    raise ValueError(f"Location1 for {base, pairing} is not 0")
                if location2 != 0:
                    raise ValueError(f"Location2 for {base, pairing} is not 0")
                params[base, pairing] = weight1, shape1, scale1, shape2, scale2
    logpdfs_paired = dict()
    logpdfs_unpaired = dict()
    pseudoenergies = dict()
    for base in ["A", "C"]:
        logpdfs_paired[base] = DoubleGammaDistribution._calc_logpdf(x, *params[base, "Paired"])
        logpdfs_unpaired[base] = DoubleGammaDistribution._calc_logpdf(x, *params[base, "Unpaired"])
        pseudoenergies[base] = (-R * T) * (logpdfs_paired[base] - logpdfs_unpaired[base])
    return logpdfs_paired, logpdfs_unpaired, pseudoenergies


def calc_rnastructure_defaults(x: np.ndarray):
    """Calculate the default pseudoenergy for RNAStructure."""
    R = 3.8078201e-3  # kcal/mol/K
    T = 310.15  # Kelvin
    dist_file = Path(os.environ["DATAPATH"]).joinpath("dists", "DMSdist.txt")
    params = dict()
    with open(dist_file) as f:
        f.readline()
        for pairing in ["Paired", "Unpaired"]:
            (shape1,
             location1,
             scale1,
             shape2,
             location2,
             scale2,
             weight1,
             weight2) = map(float, f.readline().split())
            if abs(weight1 + weight2 - 1.0) > 1e-6:
                raise ValueError(f"Weights for {pairing} do not sum to 1.0")
            if location1 != 0:
                raise ValueError(f"Location1 for {pairing} is not 0")
            if location2 != 0:
                raise ValueError(f"Location2 for {pairing} is not 0")
            params[pairing] = weight1, shape1, scale1, shape2, scale2
    logpdfs_paired = dict()
    logpdfs_unpaired = dict()
    pseudoenergies = dict()
    logpdfs_paired = DoubleGammaDistribution._calc_logpdf(x, *params["Paired"])
    logpdfs_unpaired = DoubleGammaDistribution._calc_logpdf(x, *params["Unpaired"])
    pseudoenergies = (-R * T) * (logpdfs_paired - logpdfs_unpaired)
    return logpdfs_paired, logpdfs_unpaired, pseudoenergies


def plot_pseudoenergies(params: pd.DataFrame,
                        dist_type_paired: type[FittedDistribution],
                        dist_type_unpaired: type[FittedDistribution],
                        plot_file: Path):
    """Plot the pseudoenergies."""
    R = 1.9872041e-3  # kcal/mol/K
    T = 310.15  # Kelvin

    # Set up the plot.
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex="col", sharey="row", figsize=(width_in, height_in))

    # Calculate the default log PDFs and pseudoenergies.
    x = np.linspace(0.003, 0.250, 248)
    x_default = np.linspace(0.003, 0.997, 995)
    default_logpdfs_paired, default_logpdfs_unpaired, default_pseudoenergies = calc_rnastructure_defaults(x_default)

    # Draw the zero energy lines.
    axes[0, 0].set_title("New Parameters")
    axes[1, 0].plot([x[0], x[-1]], [0, 0], linestyle="--", color="gray", lw=2)
    axes[0, 1].set_title("Default Parameters")
    axes[1, 1].plot([x_default[0], x_default[-1]], [0, 0], linestyle="--", color="gray", lw=2)

    # Label the axes.
    axes[1, 0].set_xlabel("Mutation Rate (raw)")
    axes[1, 1].set_xlabel("Mutation Rate (normalized)")
    axes[0, 0].set_ylabel("Probability Density Function")
    axes[1, 0].set_ylabel("Pseudoenergy per Paired Base (kcal/mol)")

    # Plot the new and default PDFs and pseudoenergies.
    for base, color in zip(["A", "C"], ["C1", "C0"]):
        internal = params.loc[(params["base"] == base) & (params["pairing"] == "Internal")]
        unpaired = params.loc[(params["base"] == base) & (params["pairing"] == "Unpaired")]
        params_internal = [float(internal[f"param{i}"].iloc[0]) for i in range(dist_type_paired.get_num_params())]
        params_unpaired = [float(unpaired[f"param{i}"].iloc[0]) for i in range(dist_type_unpaired.get_num_params())]
        logpdf_internal = dist_type_paired._calc_logpdf(x, *params_internal)
        logpdf_unpaired = dist_type_unpaired._calc_logpdf(x, *params_unpaired)
        pseudoenergy = (-R * T) * (logpdf_internal - logpdf_unpaired)
        # Plot log PDFs.
        axes[0, 0].plot(x, np.exp(logpdf_internal), linestyle="--", color=color, lw=2, label=f"{base} - Internal")
        axes[0, 0].plot(x, np.exp(logpdf_unpaired), linestyle=":", color=color, lw=2, label=f"{base} - Unpaired")
        axes[0, 0].legend(loc="upper right")
        # Plot pseudoenergy.
        axes[1, 0].plot(x, pseudoenergy, linestyle="-", color=color, lw=2, label=f"{base}")
        axes[1, 0].legend(loc="upper left")
        # Plot the default log PDFs.
        axes[0, 1].plot(x_default, np.exp(default_logpdfs_paired[base]), linestyle="--", color=color, lw=2, label=f"{base} - Paired")
        axes[0, 1].plot(x_default, np.exp(default_logpdfs_unpaired[base]), linestyle=":", color=color, lw=2, label=f"{base} - Unpaired")
        axes[0, 1].legend(loc="upper right")
        # Plot the default pseudoenergy.
        axes[1, 1].plot(x_default, default_pseudoenergies[base], linestyle="-", color=color, lw=2, label=f"{base}")
        axes[1, 1].legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()


def plot_rnastructure_parameters():
    """Plot the RNAStructure parameters."""
    width_cm = 48
    height_cm = 28
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(width_in, height_in))
    ax1, ax2 = axes


def run():
    # Filter and process structure data.
    try:
        refs_to_structs = pd.read_csv(STRUCTS_FILE, index_col="reference")
    except FileNotFoundError:
        refs_to_structs = map_refs_to_structs()
        refs_to_structs.to_csv(STRUCTS_FILE, index=True)

    # Get mutation rates.
    try:
        mus_095 = pd.read_csv(MUS_FILE_095)
    except FileNotFoundError:
        mus_095 = gather_mus(refs_to_structs, min_auc=0.95)
        mus_095.to_csv(MUS_FILE_095, index=False)
    # Remove excessive mutation rates.
    mus_095 = mus_095.loc[mus_095["Mutation"] <= MAX_MU]

    #plot_counts(mus_095)
    # Fit and plot mutation rate distributions.
    #plot_mutation_cdfs(mus_095)
    #plot_mutation_histograms(mus_095)

    # Cross-validate.
    try:
        loglikes = pd.read_csv(CROSS_VALIDATE_FILE)
    except FileNotFoundError:
        loglikes = cross_validate(mus_095)
        loglikes.to_csv(CROSS_VALIDATE_FILE, index=False)
    #plot_cross_validation(loglikes)

    # Bootstrap the mutation histograms to get the parameters.
    if not BOOTSTRAP_PARAMS_N2_N2_K2_FILE.exists():
        params_n2_n2_k2 = bootstrap_mutation_histograms(
            mus_095,
            BOOTSTRAP_HISTOGRAMS_N2_N2_K2_FILE,
            {"Internal": DoubleLogitNormalDistribution,
             "Terminal": DoubleLogitNormalDistribution,
             "Unpaired": DoubleKumaraswamyDistribution}
        )
        params_n2_n2_k2.to_csv(BOOTSTRAP_PARAMS_N2_N2_K2_FILE, index=False)
    if not BOOTSTRAP_PARAMS_K1_N1_K2_FILE.exists():
        params_k1_n1_k2 = bootstrap_mutation_histograms(
            mus_095,
            BOOTSTRAP_HISTOGRAMS_K1_N1_K2_FILE,
            {"Internal": KumaraswamyDistribution,
             "Terminal": LogitNormalDistribution,
             "Unpaired": DoubleKumaraswamyDistribution}
        )
        params_k1_n1_k2.to_csv(BOOTSTRAP_PARAMS_K1_N1_K2_FILE, index=False)
    else:
        params_k1_n1_k2 = pd.read_csv(BOOTSTRAP_PARAMS_K1_N1_K2_FILE)
    plot_pseudoenergies(params_k1_n1_k2,
                        dist_type_paired=KumaraswamyDistribution,
                        dist_type_unpaired=DoubleKumaraswamyDistribution,
                        plot_file=BOOTSTRAP_PSEUDOENERGY_K1_N1_K2_FILE)
    if not BOOTSTRAP_PARAMS_G1_G1_G2_FILE.exists():
        params_g1_g1_g2 = bootstrap_mutation_histograms(
            mus_095,
            BOOTSTRAP_HISTOGRAMS_G1_G1_G2_FILE,
            {"Internal": GammaDistribution,
             "Terminal": GammaDistribution,
             "Unpaired": DoubleGammaDistribution}
        )
        params_g1_g1_g2.to_csv(BOOTSTRAP_PARAMS_G1_G1_G2_FILE, index=False)

    #

if __name__ == "__main__":
    run()

