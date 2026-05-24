import os
from pathlib import Path
from typing import Iterable

import numpy as np
from click import command

from .clusts import load_pclust
from .ends import load_pends
from .muts import load_pmut
from ..core import path
from ..core.arg import (
    DEFAULT_INJECTED_MUT_PROBS,
    DEFAULT_MIN_MUT_GAP_WEIGHTS,
    opt_param_dir,
    opt_profile_name,
    opt_sample_sim,
    opt_branch,
    opt_paired_end,
    opt_read_length,
    opt_reverse_fraction,
    opt_probe,
    opt_min_mut_gap,
    opt_min_mut_gap_weights,
    opt_mut_collisions,
    opt_injected_mut_probs,
    opt_num_reads,
    opt_batch_size,
    opt_write_read_names,
    opt_brotli_level,
    opt_force,
    opt_num_cpus,
    opt_seed,
)
from ..core.logs import logger
from ..core.rna import find_ct_region
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..filter.main import set_mut_gap_params
from ..idmut.sim import simulate_idmut

COMMAND = __name__.split(os.path.extsep)[-1]


def set_sim_mut_params(
    probe: str,
    min_mut_gap_weights: str | None = None,
    injected_mut_probs: str | None = None,
):
    """Resolve simulation mutation parameters based on probe type.

    Parameters
    ----------
    probe: str
        Probe type (one of the values in ``PROBES``), used to set
        defaults when a parameter is ``None``.
    min_mut_gap_weights: str or None, optional
        Comma-separated gap:weight pairs; if None, a probe-specific
        default is used. Pass an empty string to disable.
    injected_mut_probs: str or None, optional
        Comma-separated probabilities for injecting mutations 5' of
        existing mutations; if None, a probe-specific default is used.
        Pass an empty string to disable.

    Returns
    -------
    tuple[str, str]
        Resolved ``(min_mut_gap_weights, injected_mut_probs)`` strings.
    """
    if min_mut_gap_weights is None:
        min_mut_gap_weights = DEFAULT_MIN_MUT_GAP_WEIGHTS[probe]
        logger.detail(
            f"Auto-selected min_mut_gap_weights={repr(min_mut_gap_weights)} "
            f"for probe {repr(probe)}"
        )
    if injected_mut_probs is None:
        injected_mut_probs = DEFAULT_INJECTED_MUT_PROBS[probe]
        logger.detail(
            f"Auto-selected injected_mut_probs={repr(injected_mut_probs)} "
            f"for probe {repr(probe)}"
        )
    return min_mut_gap_weights, injected_mut_probs


def parse_min_mut_gap_weights(min_mut_gap_weights: str) -> dict[int, float]:
    """Parse a comma-separated 'gap:weight' string into a dict."""
    weights = dict()
    for pair in min_mut_gap_weights.split(","):
        if not pair.strip():
            continue
        items = tuple(pair.split(":"))
        if len(items) != 2:
            raise ValueError(
                f"Each pair must have exactly 2 items, but got {len(items)} ({items})"
            )
        try:
            gap = int(items[0])
            weight = float(items[1])
        except ValueError:
            raise ValueError(
                f"Each pair must be an integer and a float, but got {items}"
            ) from None
        if gap < 0:
            raise ValueError(f"gap must be ≥ 0, but got {gap}")
        if not 0 <= weight <= 1:
            raise ValueError(f"weight must be in [0, 1], but got {weight}")
        if gap in weights:
            raise ValueError(f"gap {gap} is repeated")
        weights[gap] = weight
    if weights:
        weights_sum = sum(weights.values())
        if not np.isclose(weights_sum, 1.0):
            raise ValueError(
                f"weights must sum to 1, but got {weights_sum} ({weights})"
            )
    return {gap: weight for gap in sorted(weights) if (weight := weights[gap]) > 0.0}


def parse_injected_mut_probs(injected_mut_probs: str) -> dict[int, float]:
    """Parse a comma-separated 'offset:prob' string into a dict."""
    probs = dict()
    for pair in injected_mut_probs.split(","):
        if not pair.strip():
            continue
        items = tuple(pair.split(":"))
        if len(items) != 2:
            raise ValueError(
                f"Each pair must have exactly 2 items, but got {len(items)} ({items})"
            )
        try:
            offset = int(items[0])
            prob = float(items[1])
        except ValueError:
            raise ValueError(
                f"Each pair must be an integer and a float, but got {items}"
            ) from None
        if offset < 1:
            raise ValueError(f"offset must be ≥ 1, but got {offset}")
        if not 0 <= prob <= 1:
            raise ValueError(f"prob must be in [0, 1], but got {prob}")
        if offset in probs:
            raise ValueError(f"offset {offset} is repeated")
        probs[offset] = prob
    return {offset: prob for offset in sorted(probs) if (prob := probs[offset]) > 0.0}


def _get_param_dir_fields(param_dir: Path):
    fields = path.parse(param_dir, [path.RefSeg, path.RegSeg])
    params_dir = fields[path.TOP]
    if params_dir.name != path.SIM_PARAM_DIR:
        raise ValueError(
            f"Expected parameter directory named {repr(path.SIM_PARAM_DIR)}, "
            f"but got {repr(params_dir.name)}"
        )
    return params_dir.parent, fields[path.REF], fields[path.REG]


def _load_param_dir(param_dir: Path, profile: str):
    """Load all parameters for a profile in a directory."""
    prefix = param_dir.joinpath(profile)
    region = find_ct_region(prefix.with_suffix(path.CT_EXT))
    pmut = load_pmut(prefix.with_suffix(path.PARAM_MUTS_EXT))
    u5s, u3s, pends = load_pends(prefix.with_suffix(path.PARAM_ENDS_EXT))
    pclust = load_pclust(prefix.with_suffix(path.PARAM_CLUSTS_EXT))
    return region, pmut, u5s, u3s, pends, pclust


def _from_param_dir(param_dir: Path, profile: str, **kwargs):
    """Simulate an IDmut dataset given parameter files."""
    sim_dir, _, _ = _get_param_dir_fields(param_dir)
    region, pmut, u5s, u3s, pends, pclust = _load_param_dir(param_dir, profile)
    return simulate_idmut(
        out_dir=sim_dir.joinpath(path.SIM_SAMPLES_DIR),
        ref=region.ref,
        refseq=region.seq,
        pmut=pmut,
        uniq_end5s=u5s,
        uniq_end3s=u3s,
        pends=pends,
        pclust=pclust,
        **kwargs,
    )


@run_func(COMMAND, with_tmp=True)
def run(
    *,
    param_dir: Iterable[str | Path],
    profile_name: str,
    sample: str,
    branch: str,
    paired_end: bool,
    read_length: int,
    reverse_fraction: float,
    probe: str,
    min_mut_gap: int | None,
    min_mut_gap_weights: str | None,
    mut_collisions: str,
    injected_mut_probs: str | None,
    num_reads: int,
    batch_size: int,
    write_read_names: bool,
    brotli_level: int,
    tmp_dir: Path,
    force: bool,
    num_cpus: int,
    seed: int | None,
):
    """Simulate an IDmut dataset."""
    min_mut_gap, mut_collisions = set_mut_gap_params(probe, min_mut_gap, mut_collisions)
    min_mut_gap_weights, injected_mut_probs = set_sim_mut_params(
        probe, min_mut_gap_weights, injected_mut_probs
    )
    injected_mut_probs_dict = parse_injected_mut_probs(injected_mut_probs)
    min_mut_gap_weights_dict = parse_min_mut_gap_weights(min_mut_gap_weights)
    return dispatch(
        _from_param_dir,
        num_cpus=num_cpus,
        pass_num_cpus=False,
        as_list=True,
        ordered=False,
        raise_on_error=False,
        args=as_list_of_tuples(map(Path, param_dir)),
        kwargs=dict(
            sample=sample,
            branch=branch,
            profile=profile_name,
            paired=paired_end,
            read_length=read_length,
            p_rev=reverse_fraction,
            min_mut_gap=min_mut_gap,
            min_mut_gap_weights=min_mut_gap_weights_dict,
            injected_mut_probs=injected_mut_probs_dict,
            mut_collisions=mut_collisions,
            num_reads=num_reads,
            batch_size=batch_size,
            write_read_names=write_read_names,
            brotli_level=brotli_level,
            tmp_dir=tmp_dir,
            force=force,
            seed=seed,
        ),
    )


params = [
    opt_param_dir,
    opt_profile_name,
    opt_sample_sim,
    opt_branch,
    opt_paired_end,
    opt_read_length,
    opt_reverse_fraction,
    opt_probe,
    opt_min_mut_gap,
    opt_min_mut_gap_weights,
    opt_mut_collisions,
    opt_injected_mut_probs,
    opt_num_reads,
    opt_batch_size,
    opt_write_read_names,
    opt_brotli_level,
    opt_force,
    opt_num_cpus,
    opt_seed,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """Simulate an IDmut dataset."""
    run(*args, **kwargs)
