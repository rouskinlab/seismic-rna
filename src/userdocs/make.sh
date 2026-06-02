#!/bin/bash

# This script MUST be run from the directory containing it, i.e.
# cd src/userdocs; bash make.sh

set -eu -o pipefail

# Directories
ROOT_DIR=$PWD/../..
BUILD_DIR=$ROOT_DIR/docs
LOGO_DIR=$ROOT_DIR/logo
SRC_DIR=$ROOT_DIR/src
SEISMIC_DIR=$SRC_DIR/seismicrna
DOCSRC_DIR=$SRC_DIR/userdocs
PYAPI_DIR=$DOCSRC_DIR/api
PYAPI_AUTODOC_DIR=$PYAPI_DIR/_autodoc

# Delete the old auto-generated Python API source files, if any
# (preserve $PYAPI_DIR/index.rst which is hand-written).
if [ -d $PYAPI_AUTODOC_DIR ]
then
    rm -rv $PYAPI_AUTODOC_DIR
fi
mkdir -p $PYAPI_AUTODOC_DIR

# Build the Python API source files into the _autodoc subdirectory.
sphinx-apidoc -s rst --no-toc --no-headings --module-first -o $PYAPI_AUTODOC_DIR $SEISMIC_DIR

# Delete the old GitHub Pages files, if any.
if [ -d $BUILD_DIR ]
then
    rm -rv $BUILD_DIR
fi
mkdir $BUILD_DIR

# Make an empty file called .nojekyll to tell GitHub Pages to not use Jekyll.
# Otherwise, some files including the style sheets are not copied.
# See https://github.blog/2009-12-29-bypassing-jekyll-on-github-pages/
touch $BUILD_DIR/.nojekyll

# Build the GitHub Pages files from the source files.
sphinx-build -b html $DOCSRC_DIR $BUILD_DIR

# Brand color (#a00821) is set via html_theme_options in conf.py.
