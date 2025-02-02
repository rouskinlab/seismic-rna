from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

from .emk import EMRunsK
from .jackpot import calc_semi_g_anomaly
from .names import LOG_EXP_NAME, LOG_OBS_NAME, SEMI_G_ANOMALY, JACKPOT_QUOTIENT
from .uniq import UniqReads
from ..core import path
from ..core.header import NUM_CLUSTS_NAME
from ..core.logs import logger

PRECISION = 3


def format_exp_count_col(k: int):
    return f"{LOG_EXP_NAME}: {NUM_CLUSTS_NAME}={k}"


def assemble_log_obs_exp(uniq_reads: UniqReads, ks: list[EMRunsK]):
    """ Assemble the expected and observed log counts of each read. """
    # For each number of clusters, compute the expected log counts.
    log_exp = [(format_exp_count_col(runs.k),
                pd.Series(runs.best.log_exp, uniq_reads.uniq_names))
               for runs in ks]
    # Assemble all log counts into one DataFrame.
    log_obs_exp = pd.DataFrame.from_dict({LOG_OBS_NAME: uniq_reads.log_obs}
                                         | dict(log_exp))
    # Sort by the read identifier (index) and round the log counts.
    return log_obs_exp.sort_index().round(PRECISION)


def table_log_obs_exp(log_obs_exp: pd.DataFrame, to_dir: Path):
    """ Write the expected and observed log counts of unique reads. """
    file = to_dir.joinpath(f"log-obs-exp{path.CSVZIP_EXT}")
    log_obs_exp.to_csv(file)
    return file


def graph_log_obs_exp(log_obs_exp: pd.DataFrame,
                      ks: list[EMRunsK],
                      to_dir: Path):
    """ Graph the expected vs. observed log counts of unique reads. """
    log_obs = log_obs_exp[LOG_OBS_NAME]
    num_obs = np.exp(log_obs)
    max_log_obs = np.nanmax(log_obs) if log_obs.size > 0 else np.nan
    for runs in ks:
        try:
            k = runs.k
            column = format_exp_count_col(k)
            log_exp = log_obs_exp[column]
            semi_g_anomaly = calc_semi_g_anomaly(num_obs, log_exp)
            fig = px.scatter(log_obs_exp,
                             x=column,
                             y=LOG_OBS_NAME,
                             color=semi_g_anomaly,
                             color_continuous_scale="rdbu_r",
                             color_continuous_midpoint=0.,
                             labels={"color": SEMI_G_ANOMALY},
                             title=f"{LOG_OBS_NAME} vs. {column}")
            jackpot_quotient = runs.best.jackpot_quotient
            if not np.isnan(jackpot_quotient) and not np.isnan(max_log_obs):
                min_log_exp = np.nanmin(log_exp)
                fig.add_annotation(x=min_log_exp,
                                   y=max_log_obs,
                                   text=" = ".join(map(str,
                                                       [JACKPOT_QUOTIENT,
                                                        round(jackpot_quotient,
                                                              PRECISION)])),
                                   xanchor="left",
                                   yanchor="bottom",
                                   showarrow=False)
            file = to_dir.joinpath(f"log-obs-exp_k{k}{path.PDF_EXT}")
            fig.write_image(file)
        except Exception as error:
            logger.error(error)


def write_obs_exp_counts(uniq_reads: UniqReads, ks: list[EMRunsK], to_dir: Path):
    log_obs_exp = assemble_log_obs_exp(uniq_reads, ks)
    table_log_obs_exp(log_obs_exp, to_dir)
    graph_log_obs_exp(log_obs_exp, ks, to_dir)
