#!/bin/zsh

set -eu -o pipefail


# Remove any pre-existing intermediate files.
./_clean.sh

# Compile the latest files into a PDF.
./_build.sh

# Remove the latest intermediate files.
./_clean.sh
