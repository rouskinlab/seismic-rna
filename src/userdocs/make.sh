#!/bin/zsh

set -eu -o pipefail


# Directories
TOP_DIR=$HOME/git/seismic-rna
GHP_DIR=$TOP_DIR/docs
IMG_DIR=$TOP_DIR/logo
SRC_DIR=$TOP_DIR/src
MOD_DIR=$SRC_DIR/seismicrna
DOC_DIR=$SRC_DIR/userdocs
API_DIR=$DOC_DIR/api

# Delete the old Python API source files, if any.
if [ -d $API_DIR ]
then
    rm -rv $API_DIR/*
fi

# Build the Python API source files.
sphinx-apidoc -s rst --no-toc --no-headings --module-first -o $API_DIR $MOD_DIR

# Delete the old GitHub Pages files, if any.
if [ -d $GHP_DIR ]
then
    rm -rv $GHP_DIR/*
fi

# Make an empty file called .nojekyll to tell GitHub Pages to not use Jekyll.
# Otherwise, some files including the style sheets are not copied.
# See https://github.blog/2009-12-29-bypassing-jekyll-on-github-pages/
touch $GHP_DIR/.nojekyll

# Build the GitHub Pages files from the source files.
sphinx-build -b html $DOC_DIR $GHP_DIR
