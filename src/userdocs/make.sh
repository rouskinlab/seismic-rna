#!/bin/zsh

set -eu -o pipefail


# Directories
TOP_DIR=$HOME/git/seismic-rna
GHP_DIR=$TOP_DIR/docs
SRC_DIR=$TOP_DIR/src
MOD_DIR=$SRC_DIR/seismicrna
DOC_DIR=$SRC_DIR/userdocs
API_DIR=$DOC_DIR/api

# Delete the old GitHub Pages and Python API source files.
rm -rv $GHP_DIR $API_DIR

# Build the Python API source files.
sphinx-apidoc -s rst --no-toc --no-headings --module-first -o $API_DIR $MOD_DIR

# Build the GitHub Pages files from the source files.
sphinx-build -b html $DOC_DIR $GHP_DIR
