#!/usr/bin/env bash
# Comprehensive CLI test for SEISMIC-RNA.
# Usage: bash test_seismic_cli.sh [WORKDIR]
# Default WORKDIR: /tmp/seismic-rna-test
set -eu -o pipefail

WORKDIR="${1:-tmp_seismic_cli_test}"

# -- Directories ---------------------------------------------------------------
SIM_DIR="$WORKDIR/sim"
OUT_DIR="$WORKDIR/out"
WF_OUT="$WORKDIR/wf"        # seismic wf outputs (separate subdirs per variant)

# -- Simulation parameters -----------------------------------------------------
REFS="test-refs"
REF="test-ref"
SAMPLE1="sample1"
SAMPLE2="sample2"
REFLEN=100
NUM_READS=5000
MUT_RATES="-u am 0.1 -u cm 0.1 -u gm 0.1 -u tm 0.1"

# -- Derived paths (populated after simulation) --------------------------------
FASTA="$SIM_DIR/refs/$REFS.fa"
CT_FILE="$SIM_DIR/params/$REF/full/simulated.ct"
PARAM_DIR="$SIM_DIR/params/$REF/full"
SAMPLES_DIR="$SIM_DIR/samples"
BAM1="$OUT_DIR/$SAMPLE1/align/$REF.bam"
BAM2="$OUT_DIR/$SAMPLE2/align/$REF.bam"
IDMUT1="$OUT_DIR/$SAMPLE1/idmut/$REF"
IDMUT2="$OUT_DIR/$SAMPLE2/idmut/$REF"
FILTER1="$OUT_DIR/$SAMPLE1/filter/$REF/full"
FILTER1_70="$OUT_DIR/$SAMPLE1/filter/$REF/1-70"
CLUSTER1="$OUT_DIR/$SAMPLE1/cluster/$REF"

# -- Helpers -------------------------------------------------------------------
step()    { echo; echo "=== $* ==="; }
substep() { echo "  -- $*"; }
ok()      { echo "  [OK] $*"; }

echo "SEISMIC-RNA comprehensive CLI test"
echo "WORKDIR: $WORKDIR"
mkdir "$WORKDIR"

# ===============================================================================
# PHASE 1 - seismic sim total  (guaranteed 2-structure data via retry logic)
# ===============================================================================
step "Phase 1: seismic sim total"

substep "1a. sim total - sample1 (primary data generation)"
seismic --log "" --exit-on-error sim total \
    --sim-dir "$SIM_DIR" \
    --refs "$REFS" \
    --ref  "$REF"  \
    --reflen "$REFLEN" \
    --sample "$SAMPLE1" \
    --num-reads $NUM_READS \
    --fold-min 2 --fold-max 2 \
    $MUT_RATES \
    --force --seed 0
ok "sim total (sample1)"

# ===============================================================================
# PHASE 2 - Individual sim subcommands  (re-run on existing data with --force)
# ===============================================================================
step "Phase 2: Individual sim subcommands"

substep "2a. sim ref"
seismic --log "" --exit-on-error sim ref \
    --sim-dir "$SIM_DIR" \
    --refs "$REFS" \
    --ref  "$REF"  \
    --reflen "$REFLEN" \
    --force --seed 0
ok "sim ref"

substep "2b. sim fold"
seismic --log "" --exit-on-error sim fold \
    "$FASTA" \
    --sim-dir "$SIM_DIR" \
    --fold-min 1 --fold-max 2 \
    --force
ok "sim fold"

substep "2c. sim muts"
seismic --log "" --exit-on-error sim muts \
    --ct-file "$CT_FILE" \
    $MUT_RATES \
    --force --seed 0
ok "sim muts"

substep "2d. sim ends"
seismic --log "" --exit-on-error sim ends \
    --ct-file "$CT_FILE" \
    --force
ok "sim ends"

substep "2e. sim clusts"
seismic --log "" --exit-on-error sim clusts \
    --ct-file "$CT_FILE" \
    --force --seed 0
ok "sim clusts"

substep "2f. sim params  (runs muts + ends + clusts together)"
seismic --log "" --exit-on-error sim params \
    --ct-file "$CT_FILE" \
    $MUT_RATES \
    --force --seed 0
ok "sim params"

substep "2g. sim fastq - sample1 (overwrite)"
seismic --log "" --exit-on-error sim fastq \
    --param-dir "$PARAM_DIR" \
    --sample "$SAMPLE1" \
    --num-reads $NUM_READS \
    --force --seed 0
ok "sim fastq (sample1)"

substep "2h. sim fastq - sample2 (new)"
seismic --log "" --exit-on-error sim fastq \
    --param-dir "$PARAM_DIR" \
    --sample "$SAMPLE2" \
    --num-reads $NUM_READS \
    --seed 1
ok "sim fastq (sample2)"

substep "2i. sim idmut - sample1"
seismic --log "" --exit-on-error sim idmut \
    --param-dir "$PARAM_DIR" \
    --sample "$SAMPLE1" \
    --num-reads $NUM_READS \
    --force --seed 0
ok "sim idmut (sample1)"

# ===============================================================================
# PHASE 3 - seismic align
# ===============================================================================
step "Phase 3: seismic align"

substep "3a. align - sample1 (local + fastp, default)"
seismic --log "" --exit-on-error align \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$OUT_DIR"
ok "align (sample1, default)"

substep "3b. align - sample2"
seismic --log "" --exit-on-error align \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE2" \
    --out-dir  "$OUT_DIR"
ok "align (sample2)"

substep "3c. align - sample1, end-to-end alignment"
seismic --log "" --exit-on-error align \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WORKDIR/align-e2e" \
    --bt2-end-to-end --force
ok "align (end-to-end)"

substep "3d. align - sample1, no fastp"
seismic --log "" --exit-on-error align \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WORKDIR/align-nofastp" \
    --no-fastp --force
ok "align (no fastp)"

# ===============================================================================
# PHASE 4 - seismic idmut
# ===============================================================================
step "Phase 4: seismic idmut"

substep "4a. idmut - sample1"
seismic --log "" --exit-on-error idmut \
    "$FASTA" \
    "$BAM1" \
    --out-dir "$OUT_DIR"
ok "idmut (sample1)"

substep "4b. idmut - sample2"
seismic --log "" --exit-on-error idmut \
    "$FASTA" \
    "$BAM2" \
    --out-dir "$OUT_DIR"
ok "idmut (sample2)"

# ===============================================================================
# PHASE 5 - seismic filter
# filter writes output alongside its input; no --out-dir option.
# All full-region runs after 5a overwrite the same path, so use --force.
# ===============================================================================
step "Phase 5: seismic filter"

# -- 5A. Probe variants --------------------------------------------------------
substep "5a. filter - DMS probe, full region (canonical)"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS
ok "filter (DMS, full)"

substep "5b. filter - DMS probe, region 1-70"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --region-coords "$REF" 1 70 \
    --probe DMS
ok "filter (DMS, region 1-70)"

substep "5c. filter - SHAPE probe"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe SHAPE \
    --force
ok "filter (SHAPE)"

substep "5d. filter - ETC probe"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe ETC \
    --force
ok "filter (ETC)"

substep "5e. filter - no probe"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe none \
    --force
ok "filter (none)"

# -- 5B. Read filter variants --------------------------------------------------
substep "5f. filter - --min-ncov-read 10"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --min-ncov-read 10 \
    --force
ok "filter (min-ncov-read 10)"

substep "5g. filter - --min-fcov-read 0.5"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --min-fcov-read 0.5 \
    --force
ok "filter (min-fcov-read 0.5)"

substep "5h. filter - --min-finfo-read 0.8"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --min-finfo-read 0.8 \
    --force
ok "filter (min-finfo-read 0.8)"

substep "5i. filter - --max-fmut-read 0.2"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --max-fmut-read 0.2 \
    --force
ok "filter (max-fmut-read 0.2)"

substep "5j. filter - --keep-discontig (requires --min-mut-gap 0)"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --keep-discontig \
    --min-mut-gap 0 \
    --force
ok "filter (keep-discontig)"

# -- 5C. Position filter variants ---------------------------------------------
substep "5k. filter - --min-ninfo-pos 100"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --min-ninfo-pos 100 \
    --force
ok "filter (min-ninfo-pos 100)"

substep "5l. filter - --max-fmut-pos 0.3"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --max-fmut-pos 0.3 \
    --force
ok "filter (max-fmut-pos 0.3)"

# -- 5D. Mutation specification variants --------------------------------------
substep "5m. filter - --no-del"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --no-del \
    --force
ok "filter (no-del)"

substep "5n. filter - --no-ins"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --no-ins \
    --force
ok "filter (no-ins)"

substep "5o. filter - --no-mut ad (don't count A-base deletions)"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --no-mut ad \
    --force
ok "filter (no-mut ad)"

substep "5p. filter - --no-mut ai (don't count A-base insertions)"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --no-mut ai \
    --force
ok "filter (no-mut ai)"

substep "5q. filter - --no-mut ag --no-mut ac (exclude A->G and A->C)"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --no-mut ag --no-mut ac \
    --force
ok "filter (no-mut ag ac)"

substep "5r. filter - --only-mut ag (only count A->G substitutions)"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --only-mut ag \
    --force
ok "filter (only-mut ag)"

substep "5s. filter - --only-mut ac (only count A->C substitutions)"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --only-mut ac \
    --force
ok "filter (only-mut ac)"

substep "5t. filter - --only-mut ag --only-mut ac"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --only-mut ag --only-mut ac \
    --force
ok "filter (only-mut ag ac)"

# -- 5E. Position masking variants --------------------------------------------
substep "5u. filter - --mask-a"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --mask-a \
    --force
ok "filter (mask-a)"

substep "5v. filter - --mask-g"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --mask-g \
    --force
ok "filter (mask-g)"

substep "5w. filter - --mask-polya 0"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --mask-polya 0 \
    --force
ok "filter (mask-polya 0)"

substep "5x. filter - --mask-pos test-ref 10"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --mask-pos "$REF" 10 \
    --force
ok "filter (mask-pos 10)"

# -- 5F. Observer bias correction ---------------------------------------------
substep "5y. filter - --exact-unbias"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --exact-unbias \
    --force
ok "filter (exact-unbias)"

# -- 5G. Mutation gap variants -------------------------------------------------
substep "5z. filter - --min-mut-gap 2"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --min-mut-gap 2 \
    --force
ok "filter (min-mut-gap 2)"

substep "5aa. filter - --min-mut-gap 1 --mut-collisions merge"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --min-mut-gap 1 \
    --mut-collisions merge \
    --force
ok "filter (mut-collisions merge)"

# Restore the canonical DMS/full filter for downstream steps.
substep "5ab. filter - restore canonical DMS full region"
seismic --log "" --exit-on-error filter \
    "$IDMUT1" \
    --probe DMS \
    --force
ok "filter (canonical restored)"

# Filter sample2 (needed for pool later).
substep "5ac. filter - sample2, DMS, full region"
seismic --log "" --exit-on-error filter \
    "$IDMUT2" \
    --probe DMS
ok "filter (sample2)"

# ===============================================================================
# PHASE 6 - seismic cluster
# cluster writes output alongside its input; no --out-dir option.
# ===============================================================================
step "Phase 6: seismic cluster"

substep "6a. cluster - full region, max 2 clusters"
seismic --log "" --exit-on-error cluster \
    "$FILTER1" \
    --max-clusters 2
ok "cluster (full, k<=2)"

substep "6b. cluster - region 1-70, max 2 clusters"
seismic --log "" --exit-on-error cluster \
    "$FILTER1_70" \
    --max-clusters 2
ok "cluster (1-70, k<=2)"

# ===============================================================================
# PHASE 7 - seismic table
# ===============================================================================
step "Phase 7: seismic table"

substep "7a. table - idmut output"
seismic --log "" --exit-on-error table \
    "$IDMUT1" \
    --force
ok "table (idmut)"

substep "7b. table - filter output"
seismic --log "" --exit-on-error table \
    "$OUT_DIR/$SAMPLE1/filter/$REF" \
    --force
ok "table (filter)"

substep "7c. table - cluster output"
seismic --log "" --exit-on-error table \
    "$CLUSTER1" \
    --force
ok "table (cluster)"

# ===============================================================================
# PHASE 8 - seismic fold  (outputs alongside input data)
# ===============================================================================
step "Phase 8: seismic fold"

substep "8a. fold - from filter tables (default: full region)"
seismic --log "" --exit-on-error fold \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "fold (full region default)"

substep "8b. fold - with --fold-coords to restrict region"
seismic --log "" --exit-on-error fold \
    "$OUT_DIR/$SAMPLE1/filter/$REF" \
    --fold-coords "$REF" 1 70
ok "fold (fold-coords 1-70)"

substep "8c. fold - with --fold-table-region"
seismic --log "" --exit-on-error fold \
    "$OUT_DIR/$SAMPLE1/filter/$REF" \
    --fold-table-region --force
ok "fold (fold-table-region)"

# ===============================================================================
# PHASE 9 - seismic draw  (outputs alongside fold data; requires RNARTISTCORE)
# ===============================================================================
step "Phase 9: seismic draw"

if [ -n "${RNARTISTCORE:-}" ]; then
    substep "9a. draw - from fold output"
    seismic --log "" --exit-on-error draw \
        "$OUT_DIR/$SAMPLE1/fold/$REF"
    ok "draw"
else
    substep "9a. draw - SKIPPED (RNARTISTCORE not set)"
fi

# ===============================================================================
# PHASE 10 - seismic graph  (all 13 subcommands)
# ===============================================================================
step "Phase 10: seismic graph subcommands"

substep "10a. graph profile - from filter (mutation profiles)"
seismic --log "" --exit-on-error graph profile \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph profile (filter)"

substep "10b. graph profile - from cluster"
seismic --log "" --exit-on-error graph profile \
    "$CLUSTER1"
ok "graph profile (cluster)"

substep "10c. graph histpos - per-position histograms"
seismic --log "" --exit-on-error graph histpos \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph histpos"

substep "10d. graph histread - per-read histograms"
seismic --log "" --exit-on-error graph histread \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph histread"

substep "10e. graph giniroll - rolling Gini coefficients"
seismic --log "" --exit-on-error graph giniroll \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph giniroll"

substep "10f. graph snrroll - rolling signal-to-noise"
seismic --log "" --exit-on-error graph snrroll \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph snrroll"

substep "10g. graph roc - ROC curves (needs fold output)"
seismic --log "" --exit-on-error graph roc \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph roc"

substep "10h. graph aucroll - rolling AUC"
seismic --log "" --exit-on-error graph aucroll \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph aucroll"

substep "10i. graph poscorr - position correlations"
seismic --log "" --exit-on-error graph poscorr \
    "$OUT_DIR/$SAMPLE1/filter"
ok "graph poscorr"

substep "10j. graph mutdist - mutation distances"
seismic --log "" --exit-on-error graph mutdist \
    "$OUT_DIR/$SAMPLE1/filter/$REF"
ok "graph mutdist"

substep "10k. graph abundance - cluster abundance"
seismic --log "" --exit-on-error graph abundance \
    "$CLUSTER1"
ok "graph abundance"

substep "10l. graph corroll - rolling correlation (pairwise)"
seismic --log "" --exit-on-error graph corroll \
    "$OUT_DIR/$SAMPLE1/filter"
ok "graph corroll"

substep "10m. graph delprof - delta profiles (pairwise)"
seismic --log "" --exit-on-error graph delprof \
    "$OUT_DIR/$SAMPLE1/filter"
ok "graph delprof"

substep "10n. graph scatter - scatter plots (pairwise)"
seismic --log "" --exit-on-error graph scatter \
    "$OUT_DIR/$SAMPLE1/filter"
ok "graph scatter"

# ===============================================================================
# PHASE 11 - seismic collate
# ===============================================================================
step "Phase 11: seismic collate"

substep "11a. collate - entire out dir"
seismic --log "" --exit-on-error collate \
    "$OUT_DIR"
ok "collate"

# ===============================================================================
# PHASE 12 - seismic wf  (multiple variants)
# ===============================================================================
step "Phase 12: seismic wf"
mkdir -p "$WF_OUT"

substep "12a. wf - basic (align + idmut + filter), DMS"
seismic --log "" --exit-on-error wf \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WF_OUT/basic" \
    --probe DMS
ok "wf (basic)"

substep "12b. wf - with --region-coords, DMS"
seismic --log "" --exit-on-error wf \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WF_OUT/region-coords" \
    --probe DMS \
    --region-coords "$REF" 1 70
ok "wf (region-coords)"

substep "12c. wf - with clustering"
seismic --log "" --exit-on-error wf \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WF_OUT/cluster" \
    --probe DMS \
    --cluster \
    --max-clusters 2
ok "wf (cluster)"

substep "12d. wf - with fold"
seismic --log "" --exit-on-error wf \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WF_OUT/fold" \
    --probe DMS \
    --fold
ok "wf (fold)"

if [ -n "${RNARTISTCORE:-}" ]; then
    substep "12e. wf - with fold + draw"
    seismic --log "" --exit-on-error wf \
        "$FASTA" \
        --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
        --out-dir  "$WF_OUT/fold-draw" \
        --probe DMS \
        --fold --draw
    ok "wf (fold + draw)"
else
    substep "12e. wf - with fold + draw - SKIPPED (RNARTISTCORE not set)"
fi

substep "12f. wf - SHAPE probe"
seismic --log "" --exit-on-error wf \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WF_OUT/shape" \
    --probe SHAPE
ok "wf (SHAPE probe)"

substep "12g. wf - end-to-end alignment"
seismic --log "" --exit-on-error wf \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --out-dir  "$WF_OUT/e2e" \
    --bt2-end-to-end \
    --probe DMS
ok "wf (end-to-end)"

substep "12h. wf - two samples, cluster + all graphs"
seismic --log "" --exit-on-error wf \
    "$FASTA" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE1" \
    --dmfastqx "$SAMPLES_DIR/$SAMPLE2" \
    --out-dir  "$WF_OUT/multi" \
    --probe DMS \
    --cluster \
    --max-clusters 2 \
    --graph-mprof \
    --graph-ncov \
    --graph-mhist \
    --graph-giniroll \
    --graph-roc \
    --graph-aucroll \
    --graph-poscorr \
    --graph-mutdist \
    --graph-abundance
ok "wf (multi-sample + all graphs)"

# ===============================================================================
# PHASE 13 - Utility commands
# ===============================================================================
step "Phase 13: Utility commands"

substep "13a. list - filter output"
seismic --log "" --exit-on-error list \
    "$OUT_DIR/$SAMPLE1/filter"
ok "list"

substep "13b. pool - vertically merge sample1 and sample2 idmut"
seismic --log "" --exit-on-error pool \
    "pooled" \
    "$IDMUT1" \
    "$IDMUT2"
ok "pool"

substep "13c. join - merge full and 1-70 filter regions horizontally"
seismic --log "" --exit-on-error join \
    "joined-1-100" \
    "$FILTER1" \
    "$FILTER1_70"
ok "join"

substep "13d. ensembles - idmut output"
seismic --log "" --exit-on-error ensembles \
    "$IDMUT1"
ok "ensembles"

substep "13e. sim abstract - from filter output"
seismic --log "" --exit-on-error sim abstract \
    "$FILTER1"
ok "sim abstract"

substep "13f. cleanfa - clean FASTA"
seismic --log "" --exit-on-error cleanfa \
    "$FASTA" \
    --out-dir "$WORKDIR/cleanfa"
ok "cleanfa"

substep "13g. ct2db - convert CT to dot-bracket"
seismic --log "" --exit-on-error ct2db \
    "$SIM_DIR/params/$REF/full"
ok "ct2db"

substep "13h. db2ct - convert dot-bracket back to CT"
seismic --log "" --exit-on-error db2ct \
    "$SIM_DIR/params/$REF/full" \
    --force
ok "db2ct"

substep "13i. renumct - renumber CT file"
seismic --log "" --exit-on-error renumct \
    --ct-pos-5 "$CT_FILE" 1 \
    --out-dir "$WORKDIR/renumct"
ok "renumct"

# ===============================================================================
rm -r "$WORKDIR"
