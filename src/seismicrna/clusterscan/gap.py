"""Validate the gap between two adjacent filterscan domains by testing
whether the domains' clusters are *independent* on the reads that span the
gap.

filterscan cuts where anti-correlated "bridge" pairs are depleted, but that
depletion happens both at a genuine independent gap and across a
positively-correlated (co-folding) region -- so a spurious cut can fragment
one module into two. This test distinguishes the two cases: it reuses the
already-computed per-domain cluster mutation-rate profiles (the expensive
part) as *fixed* components and fits only the joint cluster proportions on the
spanning reads -- a concave, stable, fast reduced EM (no re-clustering, none
of the merged-EM instability). It then compares, by BIC, an independent
(product) model against a free-association (saturated) model. If association
genuinely improves the fit, the gap is spurious and the domains should be
merged.

Because a spanning read's likelihood factorizes into an A-side part (its
A-domain positions under the A-profiles) and a B-side part, each read's
component log-likelihood is just ``A_i(a) + B_i(b)`` and both models are fit
over the same fixed components.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

from ..core.logs import logger

if TYPE_CHECKING:
    import pandas as pd
    from ..cluster.uniq import UniqReads

# Mutation rates are clipped away from 0 and 1 to keep the log-likelihood
# finite.
EPS = 1e-6
# The reduced EM over proportions is concave, so it converges quickly; these
# limits keep it snappy even when a near-empty cell drifts slowly toward zero.
_MAX_ITER = 1000
_TOL = 1e-4


def _side_loglik(mut, cov, mus: "pd.DataFrame", positions: list[int], region_end5: int):
    """Bernoulli log-likelihood of each read's covered positions on one side
    under each of that side's fixed cluster profiles.

    Returns an ``(n_reads x n_clusters)`` array.
    """
    import numpy as np

    cols = [p - region_end5 for p in positions]
    M = mut[:, cols].astype(float)  # mutated indicator
    C = cov[:, cols].astype(float)  # covered indicator
    NM = C - M  # covered & not mutated (mut implies cov)
    mu = np.clip(mus.loc[positions].to_numpy(dtype=float).T, EPS, 1.0 - EPS)  # (K x P)
    return M @ np.log(mu).T + NM @ np.log1p(-mu).T


def _fit_saturated(comp_ll, weights):
    """Fit free mixing proportions over fixed components (E/M on weights only).
    ``comp_ll`` is ``(n x K)``; returns (logL, proportions)."""
    import numpy as np
    from scipy.special import logsumexp

    n, k = comp_ll.shape
    logw = np.log(np.full(k, 1.0 / k))
    total = weights.sum()
    prev = -np.inf
    ll = prev
    for _ in range(_MAX_ITER):
        logr = logw[None, :] + comp_ll
        logz = logsumexp(logr, axis=1)
        ll = float((weights * logz).sum())
        resp = np.exp(logr - logz[:, None])
        wk = (weights[:, None] * resp).sum(axis=0)
        logw = np.log(np.clip(wk / total, 1e-300, None))
        if ll - prev < _TOL:
            break
        prev = ll
    return ll, np.exp(logw)


def _fit_product(a_ll, b_ll, weights):
    """Fit the independent product proportions ``pi_ab = u_a v_b`` over the
    fixed components (E/M on the marginals only). Returns (logL, u, v)."""
    import numpy as np
    from scipy.special import logsumexp

    n, k_a = a_ll.shape
    k_b = b_ll.shape[1]
    u = np.full(k_a, 1.0 / k_a)
    v = np.full(k_b, 1.0 / k_b)
    total = weights.sum()
    comp = a_ll[:, :, None] + b_ll[:, None, :]  # (n x k_a x k_b)
    prev = -np.inf
    ll = prev
    for _ in range(_MAX_ITER):
        logpi = np.log(u)[:, None] + np.log(v)[None, :]
        logr = logpi[None] + comp
        logz = logsumexp(logr.reshape(n, -1), axis=1)
        ll = float((weights * logz).sum())
        resp = np.exp(logr - logz[:, None, None])
        u = np.clip(
            (weights[:, None] * resp.sum(axis=2)).sum(axis=0) / total, 1e-300, None
        )
        u /= u.sum()
        v = np.clip(
            (weights[:, None] * resp.sum(axis=1)).sum(axis=0) / total, 1e-300, None
        )
        v /= v.sum()
        if ll - prev < _TOL:
            break
        prev = ll
    return ll, u, v


def evaluate_gap(
    merged_uniq_reads: "UniqReads",
    mus_a: "pd.DataFrame",
    mus_b: "pd.DataFrame",
    left_end3: int,
    right_end5: int,
    gap_min_assoc: float,
) -> tuple[bool, dict]:
    """Test whether the gap between domain A and domain B is genuine.

    Parameters
    ----------
    merged_uniq_reads: UniqReads
        Unique reads over the merged (A ∪ B) region.
    mus_a, mus_b: pd.DataFrame
        Cluster mutation rates for domain A and domain B, indexed by position
        (1-indexed) with one column per cluster.
    left_end3: int
        Last position of domain A.
    right_end5: int
        First position of domain B. May be greater than ``left_end3 + 1``: a
        gap-mode of "omit" (the default) leaves an uncalled background region
        between adjacent domains, rather than domain B starting immediately
        after domain A.
    gap_min_assoc: float
        Minimum deviation of the fitted joint proportions from independence
        (max |pi_sat - outer(u, v)|) required to call a gap spurious.

    Returns
    -------
    (passed, info): tuple[bool, dict]
        ``passed`` is True to keep the gap, False to merge the domains.
        ``info`` holds the fitted statistics for logging.
    """
    import numpy as np

    k_a = mus_a.shape[1]
    k_b = mus_b.shape[1]
    info: dict = dict(k_a=k_a, k_b=k_b)
    # With a single cluster on either side there is no association to estimate.
    if k_a < 2 or k_b < 2:
        info["reason"] = "single-cluster side"
        return True, info
    ur = merged_uniq_reads
    # A read spans the gap only if it actually reaches into both domains
    # (not merely into the background gap between them, if any).
    mask = (ur.read_end5s <= left_end3) & (ur.read_end3s >= right_end5)
    n_span = int(mask.sum())
    info["n_spanning_uniq"] = n_span
    # Too few spanning reads to populate the K_A x K_B joint cells: no evidence
    # of coupling, so keep the gap.
    if n_span < k_a * k_b:
        info["reason"] = "too few spanning reads"
        return True, info
    mut = ur.get_mut_matrix()[mask]
    cov = ur.get_cov_matrix()[mask]
    region_end5 = ur.region.end5
    unmasked = list(ur.region.unmasked_int)
    a_pos = [
        p
        for p in unmasked
        if p <= left_end3 and p in mus_a.index and bool(np.isfinite(mus_a.loc[p]).all())
    ]
    b_pos = [
        p
        for p in unmasked
        if p >= right_end5
        and p in mus_b.index
        and bool(np.isfinite(mus_b.loc[p]).all())
    ]
    if not a_pos or not b_pos:
        info["reason"] = "no usable positions on a side"
        return True, info
    a_ll = _side_loglik(mut, cov, mus_a, a_pos, region_end5)
    b_ll = _side_loglik(mut, cov, mus_b, b_pos, region_end5)
    # Unique-read weighting (each distinct read pattern counts once) so PCR
    # duplication does not inflate the evidence.
    weights = np.ones(n_span)
    ln_n = np.log(float(n_span))
    ll_product, u, v = _fit_product(a_ll, b_ll, weights)
    comp = (a_ll[:, :, None] + b_ll[:, None, :]).reshape(n_span, k_a * k_b)
    ll_saturated, w_saturated = _fit_saturated(comp, weights)
    pi_saturated = w_saturated.reshape(k_a, k_b)
    params_product = (k_a - 1) + (k_b - 1)
    params_saturated = k_a * k_b - 1
    bic_product = -2.0 * ll_product + params_product * ln_n
    bic_saturated = -2.0 * ll_saturated + params_saturated * ln_n
    assoc = float(np.max(np.abs(pi_saturated - np.outer(u, v))))
    # A gap is spurious (fail -> merge) when the free-association model beats
    # the independent model by BIC AND the fitted association is large enough
    # to matter (the effect-size guard against BIC's sensitivity at large N).
    passed = not (bic_saturated < bic_product and assoc > gap_min_assoc)
    info.update(
        n_spanning_uniq=n_span,
        bic_product=bic_product,
        bic_saturated=bic_saturated,
        coupling=bic_product - bic_saturated,
        assoc=assoc,
        passed=passed,
    )
    logger.debug(
        "Gap test at {}|{}: K_A={} K_B={} n_span={} "
        "BIC(product)={:.1f} BIC(saturated)={:.1f} assoc={:.3f} -> {}",
        left_end3,
        right_end5,
        k_a,
        k_b,
        n_span,
        bic_product,
        bic_saturated,
        assoc,
        "keep" if passed else "merge",
    )
    return passed, info
