#!/bin/zsh

set -eu -o pipefail


# Directories
TOP_DIR=$HOME/git/seismic-rna
MOD_DIR=$TOP_DIR/src/seismicrna
DOC_DIR=$TOP_DIR/docs/user
SRC_DIR=$DOC_DIR/src
API_DIR=$SRC_DIR/api
OUT_DIR=$DOC_DIR/build


# Build the Python API source files.
sphinx-apidoc -s rst --no-toc --no-headings --module-first -o $API_DIR $MOD_DIR

# Build the HTML documentation files from the source files.
sphinx-build -b html $SRC_DIR $OUT_DIR
