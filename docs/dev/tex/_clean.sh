#!/bin/zsh

set -eu -o pipefail


# Remove any existing intermediate files, leaving only source and output files.
echo "Removing existing intermediate files ..."
for file in *.aux(N) *.bbl(N) *.blg(N) *.lof(N) *.log(N) *.lot(N) *.lua(N) *.nlo(N) *.out(N) *.synctex*(N) *.timestamp(N) *.toc(N) *.xmpdata(N) *.xmpi(N)
do
    rm -v $file
done

