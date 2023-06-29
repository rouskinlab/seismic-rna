from matplotlib import pyplot as plt
import numpy as np

from dreem.core.mu import calc_mu_adj, calc_mu_obs


def obs_uniform_adj(min_mu: float, max_mu: float, step_mu: float,
                    n_pos: int, min_mut_gap: int):
    """ Calculate the observed mutation rate given that the adjusted
    rate of mutations is uniform. """
    mus = np.arange(min_mu, max_mu, step_mu)
    mus_adj = np.broadcast_to(mus, (n_pos, mus.size))
    return calc_mu_obs(mus_adj, min_mut_gap)


def scale_mus(mus_adj: np.ndarray,
              dms_levels: np.ndarray,
              min_mut_gap: int):
    if mus_adj.ndim != 1:
        raise ValueError(mus_adj.ndim)
    # Compute the first-order rate constant of each base reacting with
    # DMS, assuming the DMS is sufficiently concentrated that the RNA
    # methylation can be approximated as first-order: m = 1 - exp(-kt)
    # where m is fraction methylated, k is rate constant, and t is time.
    # Solving for the rate constant yields k = -log(1 - m) / t:
    rates = -np.log(1. - mus_adj)
    # Compute the mutation rate at each level of DMS.
    scaled_mus_adj = 1. - np.exp(-rates.reshape((-1, 1))
                                 * dms_levels.reshape((1, -1)))
    scaled_mus_obs = calc_mu_obs(scaled_mus_adj, min_mut_gap)
    return scaled_mus_obs, dms_levels


if __name__ == "__main__":
    mus = np.array(
        [0.000618, 0.000635, 0.000501, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000297, 0.000258, 9.90E-05, 0.000214, 0.0, 0.0, 0.0,
         0.0, 0.000653, 0.032032, 0.003446, 0.01077, 0.0, 0.016023, 0.0, 0.0, 0.007609, 0.007695, 0.004478, 0.0,
         0.03114, 0.0, 0.012211, 0.0, 0.0, 0.0201, 0.0, 0.0, 0.0, 0.025912, 0.044596, 0.050585, 0.0, 0.0, 0.01611, 0.0,
         0.0, 0.00478, 0.0, 0.0, 0.03147, 0.0, 0.0, 0.0, 0.043469, 0.0, 0.03958, 0.096238, 0.069482, 0.01157, 0.0,
         0.003118, 0.005384, 0.0, 0.004087, 0.0, 0.004159, 0.005352, 0.004021, 0.005685, 0.012131, 0.087332, 0.0,
         0.004355, 0.0, 0.003033, 0.007747, 0.0, 0.070469, 0.083698, 0.0, 0.005183, 0.0, 0.006795, 0.0, 0.003571,
         0.043786, 0.117954, 0.060346, 0.055846, 0.0, 0.064388, 0.0, 0.0, 0.0, 0.019454, 0.042051, 0.038631, 0.042683,
         0.0, 0.0, 0.0, 0.0, 0.002815, 0.0, 0.0, 0.0, 0.064011, 0.039948, 0.0, 0.0, 0.012705, 0.014947, 0.0, 0.017714,
         0.036835, 0.031001, 0.0, 0.024237, 0.0, 0.027261, 0.05271, 0.028318, 0.017275, 0.012951, 0.0, 0.0, 0.02213,
         0.0, 0.002602, 0.008778, 0.008386, 0.014188, 0.0, 0.017224, 0.046624, 0.053452, 0.0, 0.038046, 0.0, 0.038197,
         0.034338, 0.0, 0.042533, 0.0, 0.0, 0.011521, 0.0, 0.033208, 0.0, 0.020505, 0.015081, 0.025978, 0.0, 0.0,
         0.050031, 0.0, 0.0, 0.0, 0.026061, 0.010114, 0.02219, 0.0, 0.037478, 0.0, 0.044335, 0.06296, 0.049044,
         0.055637, 0.0, 0.031897, 0.0, 0.022962, 0.022047, 0.021793, 0.0, 0.013352, 0.0, 0.010578, 0.0, 0.0, 0.0, 0.0,
         0.060872, 0.0, 0.023116, 0.029301, 0.022334, 0.0, 0.007726, 0.018374, 0.0, 0.0, 0.0, 0.023918, 0.0, 0.016198,
         0.023367, 0.0, 0.0, 0.0, 0.0, 0.000562, 0.0, 0.000366, 0.0, 0.000303, 0.000423, 0.000219, 0.0003, 0.0,
         0.000202, 0.000199, 0.000448, 0.000691, 0.0, 0.0003, 0.000637, 0.0, 0.000187, 0.002153])
    scaled_obs, dms = scale_mus(mus, np.arange(0., 100., 1.), 3)
    modal = np.argmax(scaled_obs, axis=1) < scaled_obs.shape[1] - 1
    plt.plot(scaled_obs[modal].T)
    plt.xlabel("Relative DMS concentration")
    plt.ylabel("Modeled DMS reactivity (observed)")
    plt.show()
