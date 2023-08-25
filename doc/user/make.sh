#!/bin/zsh

set -eu -o pipefail


# Remove existing output files.
for dirname in build src/_autosummary
do
    if [ -d $dirname ]
    then
        rm -rv $dirname/
    fi
done

clear


# Build the HTML files.
make html


# Open the main index (on macOS).
open build/html/index.html

