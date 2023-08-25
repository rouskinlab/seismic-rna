#!/bin/zsh

set -eu -o pipefail


# Remove any existing intermediate files, leaving only source and output files.
echo "Removing existing intermediate files ..."
for file in *.aux(N) *.bbl(N) *.blg(N) *.err(N) *.lof(N) *.log(N) *.lot(N) *.lua(N) *.nlo(N) *.out(N) *.pdf(N) *.py(N) *.synctex*(N) *.timestamp(N) *.toc(N) *.xmpdata(N) *.xmpi(N)
do
    rm -v $file
done
