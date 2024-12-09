from pathlib import Path

import pandas as pd
import plotly.express as px

from .emk import EMRunsK, NOCONV
from .names import JACKPOT_QUOTIENT
from ..core.header import NUM_CLUSTS_NAME

EM_RUN_NAME = "Run"
K_RUN_NAMES = [NUM_CLUSTS_NAME, EM_RUN_NAME]
RUN_PASSING = "run_passing"
ATTRS = {
    RUN_PASSING: "Whether the run passed filters",
    "log_likes": "Final log likelihood",
    "bics": "Bayesian information criterion",
    "jackpot_quotients": JACKPOT_QUOTIENT,
    "min_nrmsds": "Minimum normalized RMSD between any two clusters",
    "max_pearsons": "Maximum Pearson correlation between any two clusters",
    "nrmsds_vs_best": "Normalized RMSD versus the best run",
    "pearsons_vs_best": "Pearson correlation versus the best run",
    "converged": f"Iterations ({NOCONV} if the run did not converge)",
}


def tabulate_attr(ks: list[EMRunsK], attr: str):
    """ Tabulate the values for one attribute. """
    # Cast runs.k and run from int to str to make them both categorial.
    # If runs.k is numeric, then the bars will be stacked, not grouped.
    # If run is numeric, then the run number will be indicated with a
    # color bar rather than a label.
    runs_values = dict()
    for runs in ks:
        runs_value = getattr(runs, attr)
        if callable(runs_value):
            runs_value = runs_value()
        for run, run_value in enumerate(runs_value):
            runs_values[str(runs.k), str(run)] = run_value
    if runs_values:
        runs_values = pd.Series(runs_values)
        runs_values.index.set_names(K_RUN_NAMES, inplace=True)
    else:
        runs_values = pd.Series([],
                                pd.MultiIndex.from_arrays([[], []],
                                                          names=K_RUN_NAMES))
    return runs_values


def tabulate(ks: list[EMRunsK]):
    """ Tabulate all attributes. """
    return pd.DataFrame.from_dict({title: tabulate_attr(ks, key)
                                   for key, title in ATTRS.items()})


def write_table(table: pd.DataFrame, cluster_dir: Path):
    table_file = cluster_dir.joinpath("summary.csv")
    table.to_csv(table_file)
    return table_file


def graph_attr(attr: pd.Series, passing_text: list[str] | None = None):
    """ Graph one attribute. """
    return px.bar(attr.to_frame().reset_index(names=K_RUN_NAMES),
                  x=NUM_CLUSTS_NAME,
                  y=attr.name,
                  title=attr.name,
                  color=EM_RUN_NAME,
                  text=passing_text,
                  barmode="group")


def graph_attrs(table: pd.DataFrame, to_dir: Path):
    """ Graph every attribute. """
    passing = table[ATTRS[RUN_PASSING]]
    if passing.all():
        passing_text = None
    else:
        passing_text = ["" if p else "FAIL" for p in passing]
    for key, attr in ATTRS.items():
        if key == RUN_PASSING:
            continue
        fig = graph_attr(table[attr], passing_text)
        fig.write_image(to_dir.joinpath(f"{key}.pdf"))


def write_summaries(ks: list[EMRunsK], to_dir: Path):
    table = tabulate(ks)
    write_table(table, to_dir)
    graph_attrs(table, to_dir)

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
