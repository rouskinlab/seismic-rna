#/bin/bash

set -eux -o pipefail

# Compile the C extension module with assertions and no optimization.
clang -O2 -DNDEBUG -fno-strict-overflow -Wsign-compare -Wunreachable-code -Wall -fPIC -isystem /Users/mfa/conda/envs/seismic/include -fPIC -isystem /Users/mfa/conda/envs/seismic/include -I/Users/mfa/conda/envs/seismic/include/python3.12 -c relate.c -o relate.o
clang -bundle -undefined dynamic_lookup -Wl,-rpath,/Users/mfa/conda/envs/seismic/lib -L/Users/mfa/conda/envs/seismic/lib -Wl,-rpath,/Users/mfa/conda/envs/seismic/lib -L/Users/mfa/conda/envs/seismic/lib relate.o -o relate.cpython-312-darwin.so

# Test the code.
python example.py

