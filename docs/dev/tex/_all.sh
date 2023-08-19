#!/bin/zsh

set -eu -o pipefail


# Remove any pre-existing intermediate files.
./_clean.sh

# Compile the latest files into a PDF.
./_build.sh

# Move the final PDF to the top-level directory.
mv main.pdf ../../../DEV-MANUAL.pdf

# Remove the latest intermediate files.
./_clean.sh

