# SEISMIC-RNA Release Schedule

## Planned updates

### Functionality
Addition of new functions into the software (use sparingly).

#### v0.7.0
- Introduce variation of information to quantify the reproducibility of clusters.
- Introduce option to filter out reads with less than an integer number of informative positions: `min_ninfo_read`.
- Introduce functions to simulate entire datasets for each step of the pipeline (primarily for testing purposes).


### Performance
Optimization of existing functions to improve speed, memory, disk usage, etc.

#### v0.7.1
- Decrease the sizes of the relation vector batch files by storing relation vectors as sparse data structures, i.e. `pandas.arrays.SparseArray`.

#### 

### Interface
Additions or updates to the command line and/or application programming interface for existing functions.

#### v0.7.0
- Introduce "amplicon mode" and "random mode" (default) that each bundles a set of default options.
  - Defaults for random mode (suitable for randomly fragmented and/or primed libraries):
    - `bt2_local=True`: Use local mode in `align` so that any extra flanking sequences are not misinterpreted as mutations (causing false positives).
    - `min_finfo_read=0`: No minimum fraction of informative bases is required to keep a read during `mask`.
    - `min_mut_gap=0`: Do not remove any reads for having mutations too close together during `mask` -- not because they should not be removed, but rather because the mutation unbiasing algorithm assumes that the reads come from amplicons. I aim to enable the unbiasing algorithm to handle non-amplicon data properly in a future release, but it's a non-trivial challenge.
  - Defaults for amplicon mode (suitable for PCR amplicons and other homogenous samples):
    - `bt2_local=False`: Use end-to-end mode in `align` so that mutations near the ends of reads are not dropped (causing false negatives).
    - `min_finfo_read=0.95`: At least 95% of the read must cover the amplicon to keep it during `mask`.
    - `min_mut_gap=3`: Remove any reads that have a pair of mutations separated by fewer than 3 non-mutated positions during `mask`.
    - Additionally, `mask` will not default to the full reference if no coordinates or primers are given because this behavior is generally unsuitable for amplicons, whose ends are obscured by primers.

### Debugging
Fixes for buggy parts of the software.

#### v0.7.0
- Fix the version attribute so that `seismicrna.__version__` returns the version instead of `unknown`.


### Testing
Additions or upgrades to the test suite.

#### v0.7.0
- Introduce end-to-end tests that send a small, simulated dataset through every step of the pipeline.


### Documentation
Additions or upgrades to the documentation.

#### v0.7.0
- Introduce the LaTeX skeleton of the (eventually) highly technical developer manual.
- Reorganize the files in `docs/` into `docs/dev/` (for the developer manual) and `docs/user` (for the user manual).
- Introduce the SEISMIC-RNA Release Schedule.
