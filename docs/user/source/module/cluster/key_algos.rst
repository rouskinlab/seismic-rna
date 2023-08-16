
Key algorithms
++++++++++++++++++++++++

Preprocessing the bitvector
---------------------------

The bitvector is pre-processed by the following steps:
    - Remove the G and U bases (if the option ``--include-gu`` is set to False)
    - Remove the low Mutation fraction bases (threshold defined by ``--signal-thresh``)
    - Remove the reads with too many mutations (more than 3 standard deviations from the mean)
    - Remove the reads with deletions (if the option ``--include-del`` is set to False)
    - Remove the reads with mutations closer than 4 bases apart
    - Convert the bitvector to binary format by counting substitution as 1
    - Keep only unique reads, and keep track of their occurence


Here is an example of non processed bitvector, as outputted by the vectoring module. Please refer to :ref:`bitvector` for more information about the bitvector format.
::

                        G1  C2  G3  T4  C5  T6  T7  A8  G9  T10  G11  C12  A13  T14  A15
    __index_level_0__                                                                   
    reference_1_read_1   1   1   1   1  64   1   1   1   1    1    1    1    1    1    1
    reference_1_read_2   1   1   1   1   1   1   1   1   1   16    1    1    1    1    1
    reference_1_read_3   1   1   1   1   1   1   1   1   1    1    1    1    1    1   32
    reference_1_read_4   1   1   1   1   1   1   1   1   1    1    1    1    1    1    1
    reference_1_read_5   1   1   1   1   1   1   1   1   1    1    1    1    1    1    1
    reference_1_read_6   1   1   1   1   1   1   128 1   1    1    1    1    1    1    1
    reference_1_read_7   1   1   1   1   1   1   1   1   1    1    1    1    1    1    1
    reference_1_read_8   1   1   1   1   1   1   1   1   1    1    1    256  1    1    1
    reference_1_read_9   1   4   1   1   1   1   1   1   1    1    1    1    1    1    1


After pre-processing:
::

                         C5   T7   C12  A15  occurence
    __index_level_0__                                                                   
    idx_0                 1   0    0    0    1
    idx_1                 0   1    0    0    1
    idx_3                 0   0    1    0    1
    idx_4                 0   0    0    1    1
    idx_5                 0   0    0    0    4


Clustering algorithm
--------------------

We use Expectation-Maximisation (EM) algorithm with a Bernoulli mixture model. 
This model is adapted to account for the fact that one read cannot have mutations closer than 4 bases apart (i.e. the bases mutations are not independent).

The expectation step computes the probability that each read belongs to each cluster:

.. math::

   P(x_n|\mu_k) = \frac{ \prod_{i=1}^{D} \mu_{ki}^{x_{ni}} * (1-\mu_{ki})^{1-x_{ni}} } {d_k}

with :
    - :math:`x_n` the bitvector of the read n
    - :math:`\mu_k` the mutation profile of cluster k
    - :math:`D` the number of bases in the read
    - :math:`d_k` the correction factor to account for the fact that one read cannot have consecutive mutations

The maximisation step computes the new mutation profile of each cluster :math:`\mu_k` and their proportion :math:`\pi_k`. 
We use a numerical solver to find the maximum likelihood solution.

See the following paper for more details: `Determination of RNA structural diversity and its role in HIV-1 RNA splicing  <https://www.nature.com/articles/s41586-020-2253-5>`_.