# Changelog

## 0.25.0 (2026-05-25)

### Breaking changes

#### Renamed subcommands
- `seismic relate` is now `seismic idmut` (Identify Mutations).
- `seismic mask` is now `seismic filter`.
- Output directories and report fields from older versions can be migrated automatically with `seismic migrate out -o out_new`.

#### Renamed options
- `--fold-table` → `--fold-table-region`; `--fold-full` → `--fold-full-region` (in `seismic fold`, `seismic graph`, `seismic draw`, `seismic wf`).
- `--quantile` split into `--fold-quantile` (folding) and `--graph-quantile` (graphs).
- Demultiplexing options corrected: `--mismatch-tolerence` → `--mismatch-tolerance`; `--index-tolerence` → `--index-tolerance`.
- Mutation-counting options renamed to `--no-mut` and `--only-mut` for clarity; deletions and insertions are now counted by default (`--count-del/--no-del`, `--count-ins/--no-ins`).
- In `seismic filter`, masking options renamed `--mask-*` → `--drop-*` for reads (`--drop-read`, `--drop-read-file`, `--drop-discontig`); the iteration cap is now `--max-filter-iter`.
- The bulk-base mask option `--mask-gu` is replaced by per-base options `--mask-a/--keep-a`, `--mask-c/--keep-c`, `--mask-g/--keep-g`, `--mask-u/--keep-u`.

#### Other behavior changes
- `seismic pool` now takes the name of the new pooled sample as a positional argument (previously `--pooled`).
- `seismic join` now takes the name of the new joined region as a positional argument (previously `--joined`).
- `seismic ensembles` now creates a new step directory called `ensembles` and generates an ensemble report.
- `seismic sim abstract` no longer accepts JSON files from `seismic export` -- only table files -- because it now also accepts options `--fold-coords`, `--fold-primers`, `--fold-regions-file`, and `--fold-table-region/--fold-full-region`, which can only be used with tables.

### New subcommands
- **`seismic collate`** — bundle a set of outputs (with optional graphs/SVGs) into a portable directory tree for sharing.
- **`seismic importmm`** — import Mutation Map files from RNA Framework and DRACO.

### Rewritten subcommand
- **`seismic demult`** — complete rewrite for correctness and speed; no longer requires fixed-length barcodes. New options: `--barcode`, `--read-pos`, `--allow-n/--no-allow-n`.

### New options

#### Reproducibility
- **`--seed`** — Optional seed for the random number generator, to ensure reproduciblity. Accepted by every command that uses random numbers, such as `cluster` and most `sim` subcommands.

#### Filter (formerly Mask)
- **`--probe {DMS|SHAPE|ETC|none}`** — Declare the chemical probe to be DMS, SHAPE, ETC (a carbodiimide reagent), or none (i.e. an untreated control sample). Default is `DMS`. Enables probe-specific defaults across the pipeline (see below for details).
- **Per-base masking** (`--mask-a/c/g/u`) for fine-grained control over which bases are masked. By default, setting `--probe` to `DMS` masks out G and T/U, `ETC` masks out A and C, and `SHAPE` and `none` do not mask out any bases. The defaults can be overridden by passing `--mask-a/c/g/u` or `--keep-a/c/g/u` explicitly.
- **`--mut-collisions {drop|merge|auto}`** — choose whether reads with overly close mutations are dropped (default and recommended when `--probe` is `DMS`) or whether close mutations are merged without dropping the reads (default and recommended when `--probe` is `SHAPE` or `ETC`). The `merge` strategy is also implemented in the clustering and jackpotting algorithms.
- **`--min-fcov-read`** — Minimum fractional coverage threshold per read. Recommended use case is setting to `--min-fcov-read=1` when reads come from PCR amplicons to discard reads that do not cover the entire amplicon.
- **`--self-contained`** — produce self-contained batches. Compared to the default (`--no-self-contained`), these batches are larger files but may load faster because they can be loaded without also loading the corresponding IDmut batches. Advantageous when only a small fraction of reads or positions remain after filtering, e.g. slicing a small region from a large dataset or a long transcript.

#### Cluster
- **`--jackpot-max-data`** — limit the amount of data used for jackpot simulations to prevent large time and memory consumption.
- **`--max-jackpot-sims`** — limit the number of null simulations during the jackpotting calculation to prevent infinite loops.
- **`--self-contained`** — produce self-contained batches, like the new option in the Filter step.

#### Fold (RNA structure prediction)
- **ViennaRNA backend.** Added ViennaRNA's `RNAfold` and `RNAsubopt` as folding backends. Select with `--fold-backend viennarna`. Selected by default with `--fold-backend auto` if `--probe` is `SHAPE`, `ETC`, or `none`. (RNAstructure is the default backend when `--probe` is `DMS`.)
- **ShapeKnots backend.** Added RNAstructure's `ShapeKnots` as a folding backend for pseudoknot prediction. Select with `--fold-backend rnastructure --pseudoknots`. Note that ViennaRNA does not currently support pseudoknots, so the `--pseudoknots` flag can only be used with RNAstructure.
- **Energy method.** Control how mutation rates are used during structure prediction with the `--fold-energy-method` option. Currently, RNAstructure supports `cordero` (default if `--probe` is `DMS`) and `deigan`, and ViennaRNA supports `eddy` (default if probe is `SHAPE`, `ETC`, or `none`) and `deigan`.
- **Eddy priors.** New options `--eddy-prior-paired-file` and `--eddy-prior-unpaired-file` control the prior distributions of paired and unpaired bases when `--fold-energy-method` is `eddy`.
- **Deigan slope and intercept.** New options `--deigan-slope` and `--deigan-intercept` control the energy function when `--fold-energy-method` is `deigan`. Their defaults match those used by both RNAstructure and ViennaRNA.
- `--fold-isolated` — allow isolated pairs (default is not allowed).
- `--fold-dry-run` — write fold commands and input files without executing.

#### Pool
- **`--min-pearson`** / **`--max-marcd`** — Require sample mutation rates to be sufficiently similar in order for the samples to be pooled.

#### Simulation (`seismic sim`)
- **`--probe`** support, mirroring `seismic filter`.
- **`--min-mut-gap-weights`** — Produce a sample with a mixture of minimum gaps between mutations, more similar to real DMS-MaPseq data. By default, enabled when `--probe` is `DMS`.
- **`--injected-mut-probs`** — Produce a sample with extra spurious mutations within several positions 5' of a real mutation, more similar to SHAPE data. Only used when `--mut-collisions` is set to `merge`. By default, enabled when `--probe` is `SHAPE` or `ETC`.
- **`--fold-backend`** added to `sim fold` (ViennaRNA support throughout the simulation pipeline).
- **`--fold-coords`, `--fold-primers`, `--fold-regions-file`** — explicit region selection for `sim abstract` (replaces JSON-file workflow).
- **`--max-tries`, `--max-fraction-ident`, `--max-pearson-sim`, `--min-marcd-sim`** — controls for ensuring simulated cluster diversity.
- **`--min-aucroc`** — filter profiles in `sim abstract` by minimum AUC-ROC.

#### Ensembles
- **`--threshold-divisor`** — controls ensemble-detection sensitivity. Default is 1 (consistent with previous behavior). Increasing it above 1 lowers the threshold and increases sensitivity at the risk of more false positives. Vice versa for decreasing it below 1. Must be strictly > 0.

#### Graph
- **`--terminal-pairs/--no-terminal-pairs`** — include or exclude base pairs at the ends of contiguous stems in structure-based graphs (`roc` and `aucroll`).
- **`--graph-quantile`** — quantile for normalizing values in graphs.

### Bug fixes
- Fixed `xam_to_fastq_cmd` failing when the output FASTQ file path (`fq_out`) was specified.
- Fixed a redundant calculation in the EM clustering algorithm's joint probability function (Issue #33).
