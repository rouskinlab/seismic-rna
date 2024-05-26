from pathlib import Path

from .params import (sim_pmut,
                     sim_paired,
                     make_pmut_means_paired,
                     make_pmut_means_unpaired,
                     sim_pends,
                     sim_pclust)
from ..core.seq import DNA
from ..relate.sim import simulate_relate


def test_simulate():
    sample = "mysample"
    num_reads = 2**20
    batch_size = 2**16
    ref = "myref"
    seq = DNA.random(500)
    nclust = 3
    paired = sim_paired(seq, nclust)
    pclust = sim_pclust(nclust)
    print(pclust)

    u5s, u3s, pends = sim_pends(1, len(seq), len(seq) * 0.8, 250, 0.1)
    pm = make_pmut_means_paired()
    um = make_pmut_means_unpaired()

    pmut = [sim_pmut(paired[cluster], pm, um, 0.05, 0.05) for cluster in paired.columns]

    out_dir = Path.cwd().joinpath("out")
    report_file = simulate_relate(out_dir, sample, ref, seq, batch_size, num_reads, pmut, u5s, u3s, pends, pclust.values,
                                  brotli_level=10, force=True)
