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

# Delete the old GitHub Pages and Python API source files.
rm -rv $GHP_DIR $API_DIR

# Build the Python API source files.
sphinx-apidoc -s rst --no-toc --no-headings --module-first -o $API_DIR $MOD_DIR

# Build the GitHub Pages files from the source files.
sphinx-build -b html $DOC_DIR $GHP_DIR

# Copy the logo images to the GitHub Pages directory.
cp $IMG_DIR/logo.png $IMG_DIR/logo-48.png $GHP_DIR

# Make an empty file called .nojekyll to tell GitHub Pages to not use Jekyll.
# Otherwise, some files including the style sheets are not copied.
# See https://github.blog/2009-12-29-bypassing-jekyll-on-github-pages/
touch $GHP_DIR/.nojekyll

