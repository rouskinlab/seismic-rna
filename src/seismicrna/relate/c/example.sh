#/bin/bash

set -eux -o pipefail

# Compile the C extension module.
python setup.py build

# Move the SO file to the current directory to be imported.
mv build/lib.macosx-10.13-x86_64-cpython-312/relate.cpython-312-darwin.so .
rm -r build

# Test the code.
python test.py

