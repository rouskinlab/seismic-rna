#!/bin/bash

# Run all unit tests for SEISMIC-RNA.
# If all tests succeed, return 0; otherwise, return 1.

set -euo pipefail

# Run the tests and capture the results.
RESULTS=test-results.txt
seismic +test 2> $RESULTS
cat $RESULTS

# Check whether all the tests succeeded.
LAST=$(tail -n 1 $RESULTS)
if [ "$LAST" == "OK" ];
then
	  EXIT=0
else
	  EXIT=1
fi

# Clean up the output files.
rm $RESULTS
if [ -d log ]
then
	  rm -r log
fi

# Exit 0 if all tests succeeded, otherwise exit 1.
exit $EXIT
