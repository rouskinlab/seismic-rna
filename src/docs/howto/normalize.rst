
Normalize Mutation Rates
================================================================================

Some applications (e.g. predicting structures, computing differences between
mutation rates) require that the mutation rates be normalized first.

.. _norm_quantile:

How to normalize mutation rates to a quantile
--------------------------------------------------------------------------------

You can normalize a mutational profile by choosing a quantile of the mutation
rates to set to 1 using ``--quantile`` (``-q``).
All mutation rates less than that quantile will then be scaled up linearly,
while mutation rates greater than that quantile will be set to 1.

.. image::
    normalize.png

For example, if you chose quantile 1, then the value at quantile 1 (i.e. the
maximum mutation rate) would be set to 1, and all other mutation rates would be
scaled proportionally.
If the maximum mutation rate were 0.20, for instance, then it would be scaled to
1.00; proportionally, other mutations rates of 0.01, 0.05, and 0.10 would be
scaled to 0.05, 0.25, and 0.50, respectively.

If you chose quantile 0.5 instead, then the value at quantile 0.5 (i.e. the
median mutation rate) would be set to 1; all other mutation rates less than the
median would be scaled proportionally, and all greater than the median would be
capped at 1.
If the median mutation rate were 0.10, for instance, then it would be scaled to
1.00; proportionally, other mutations rates of 0.01 and 0.05 would be scaled to
0.10 and 0.50, respectively, while a mutation rate of 0.20 would be set to 1.00.
