import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

from .compare import EMRunsK
from .names import LOG_EXP_NAME, LOG_OBS_NAME
from .uniq import UniqReads
from ..core import path
from ..core.header import NUM_CLUSTS_NAME

EXP_COUNT_PRECISION = 3  # Number of digits to round expected log counts


def format_exp_count_col(k: int):
    return f"{LOG_EXP_NAME}: {NUM_CLUSTS_NAME}={k}"


def parse_exp_count_col(col: str):
    match = re.match(f"^{LOG_EXP_NAME}: {NUM_CLUSTS_NAME}=([0-9]+)$", col)
    if not match:
        raise ValueError(f"Invalid expected count column: {repr(col)}")
    k, = match.groups()
    return int(k)


def calc_log_obs_exp(uniq_reads: UniqReads, ks: list[EMRunsK]):
    """ Calculate the expected and observed log counts of each read. """
    # Compute the observed log counts of the bit vectors.
    log_obs = pd.Series(np.log(uniq_reads.counts_per_uniq),
                        index=uniq_reads.get_uniq_names())
    # For each number of clusters, compute the expected log counts.
    log_exp = [(format_exp_count_col(runs.k),
                pd.Series(runs.best.logn_exp, index=log_obs.index))
               for runs in ks]
    # Assemble all log counts into one DataFrame.
    log_obs_exp = pd.DataFrame.from_dict({LOG_OBS_NAME: log_obs}
                                         | dict(log_exp))
    # Sort by the read identifier (index) and round the log counts.
    return log_obs_exp.sort_index().round(EXP_COUNT_PRECISION)


def write_log_obs_exp(log_obs_exp: pd.DataFrame, to_dir: Path):
    """ Write the expected and observed log counts of unique reads. """
    file = to_dir.joinpath(f"read-counts{path.CSVZIP_EXT}")
    log_obs_exp.to_csv(file)
    return file


def graph_log_obs_exp(log_obs_exp: pd.DataFrame, to_dir: Path):
    """ Graph the expected vs. observed log counts of unique reads. """
    for column in log_obs_exp.columns:
        if column != LOG_OBS_NAME:
            k = parse_exp_count_col(column)
            fig = px.scatter(log_obs_exp,
                             x=column,
                             y=LOG_OBS_NAME,
                             title=f"{LOG_OBS_NAME} vs. {column}")
            file = to_dir.joinpath(f"log-obs-exp_k{k}{path.PDF_EXT}")
            fig.write_image(file)


def write_jackpotting(uniq_reads: UniqReads, ks: list[EMRunsK], to_dir: Path):
    log_obs_exp = calc_log_obs_exp(uniq_reads, ks)
    write_log_obs_exp(log_obs_exp, to_dir)
    graph_log_obs_exp(log_obs_exp, to_dir)
