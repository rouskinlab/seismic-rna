#!/bin/bash

# Run all unit tests for SEISMIC-RNA.
# If all tests succeed, return 0; otherwise, return 1.

set -eu -o pipefail

# Run the tests and capture the results.
RESULTS=test-results.txt
seismic test 2> $RESULTS
cat $RESULTS
echo

# Check whether all the tests succeeded.
STATUS=$(tail -n 3 $RESULTS | head -n 1)
echo "STATUS=$STATUS"
if [ "$STATUS" == "OK" ];
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
