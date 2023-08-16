#!/bin/zsh

set -eu -o pipefail


PREFIX=main


build_main () {
    clear
    pdflatex $PREFIX.tex
}


build_refs () {
    clear
    bibtex $PREFIX
}


# Compile the new files.
build_main
#build_refs
build_main
build_main
