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

# Delete the old Python API source files, if any.
if [ -d $PYAPI_DIR ]
then
    rm -rv $PYAPI_DIR
fi
mkdir $PYAPI_DIR

# Build the Python API source files.
sphinx-apidoc -s rst --no-toc --no-headings --module-first -o $PYAPI_DIR $SEISMIC_DIR

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

# Change the color of the upper left sidebar from the default (#2980b9)
# to the shade of blue in the logo (#2186d9).
DEFAULT_NAV_BACKGROUND="#2980b9"
CUSTOM_NAV_BACKGROUND="#2186d9"
for STYLESHEET in theme badge_only
do
  sed -i "s/$DEFAULT_NAV_BACKGROUND/$CUSTOM_NAV_BACKGROUND/g" $BUILD_DIR/_static/css/$STYLESHEET.css
done
