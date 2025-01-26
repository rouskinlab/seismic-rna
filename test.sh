#!/bin/bash

# Run all unit tests for SEISMIC-RNA.
# If all tests succeed, return 0; otherwise, return 1.

set -eu -o pipefail

# Run the tests; log the stderr to the console and the $RESULTS file.
RESULTS=test-results.txt
seismic -q test -vv 2>&1 | tee $RESULTS
cat $RESULTS
echo

# Check whether all the tests succeeded.
STATUS=$(tail -n 1 $RESULTS)
if [[ $STATUS == OK* ]];
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
echo "Exiting test suite with status $EXIT"
exit $EXIT

