from pathlib import Path

import pandas as pd
import plotly.express as px

from .compare import EMRunsK, NOCONV
from ..core.header import NUM_CLUSTS_NAME

EM_RUN_NAME = "Run"
K_RUN_NAMES = [NUM_CLUSTS_NAME, EM_RUN_NAME]
RUN_PASSING = "run_passing"
ATTRS = {
    RUN_PASSING: "Whether the run passed filters",
    "log_likes": "Final log likelihood",
    "bics": "Bayesian information criterion",
    "min_nrmsds": "Minimum NRMSD between any two clusters",
    "max_pearsons": "Maximum Pearson correlation between any two clusters",
    "nrmsds_vs_best": "NRMSD versus the best run",
    "pearsons_vs_best": "Pearson correlation versus the best run",
    "converged": f"Iterations ({NOCONV} if the run did not converge)",
}


def tabulate_attr(ks: list[EMRunsK], attr: str):
    """ Tabulate the values for one attribute. """
    # Cast runs.k and run from int to str to make them both categorial.
    # If runs.k is numeric, then the bars will be stacked, not grouped.
    # If run is numeric, then the run number will be indicated with a
    # color bar rather than a label.
    values = pd.Series({(str(runs.k), str(run)): value
                        for runs in ks
                        for run, value in enumerate(getattr(runs, attr))})
    values.index.set_names(K_RUN_NAMES)
    return values


def tabulate(ks: list[EMRunsK]):
    """ Tabulate all attributes. """
    return pd.DataFrame.from_dict({title: tabulate_attr(ks, key)
                                   for key, title in ATTRS.items()})


def write_table(table: pd.DataFrame, cluster_dir: Path):
    table_file = cluster_dir.joinpath("graph.csv")
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


def graph_attrs(table: pd.DataFrame, cluster_dir: Path):
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
        fig.write_html(cluster_dir.joinpath(f"graph_{key}.html"))


def write_results(ks: list[EMRunsK], cluster_dir: Path):
    table = tabulate(ks)
    write_table(table, cluster_dir)
    graph_attrs(table, cluster_dir)
