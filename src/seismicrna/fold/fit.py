from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from click import command

from .main import load_foldable_tables
from ..core.arg import (CMD_FOLD,
                        arg_input_path,
                        opt_verify_times,
                        opt_num_cpus,
                        opt_force,
                        extra_defaults)
from ..core.rna import DEFAULT_MEAN, DEFAULT_COV, fit_beta_mixture_model
from ..core.run import run_func
from ..core.seq import BASE_NAME, DNA
from ..core.table import MUTAT_REL


@run_func(CMD_FOLD, extra_defaults=extra_defaults)
def run(input_path: Iterable[str | Path], *,
        verify_times: bool,
        num_cpus: int,
        force: bool):
    """ Predict RNA secondary structures using mutation rates. """
    mus = {base: np.array([], dtype=float) for base in DNA.four()}
    for table in load_foldable_tables(input_path, verify_times=verify_times):
        table_mus = table.fetch_ratio(rel=MUTAT_REL)
        for base in list(mus):
            base_mus = table_mus.loc[
                table_mus.index.get_level_values(BASE_NAME) == base
                ]
            for cluster in base_mus.columns:
                mus[base] = np.concatenate([mus[base],
                                            base_mus[cluster].dropna()])
    beta_params = pd.DataFrame.from_dict(
        {base: fit_beta_mixture_model(base_mus,
                                      DEFAULT_MEAN[base],
                                      DEFAULT_COV[base])
         for base, base_mus in mus.items()},
        orient="index"
    )
    print(beta_params)
    beta_params.to_csv("beta-params.csv")


params = [
    arg_input_path,
    opt_verify_times,
    opt_num_cpus,
    opt_force,
]


@command("fit", params=params)
def cli(*args, **kwargs):
    """ Fit parameters to mutational profiles. """
    return run(*args, **kwargs)
