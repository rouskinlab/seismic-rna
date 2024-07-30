#!/bin/bash

# Run all unit tests for SEISMIC-RNA.
# If all tests succeed, return 0; otherwise, return 1.

set -euxo pipefail

RESULTS=test-results.txt

seismic +test -vv 2> $RESULTS
#cat $RESULTS
tail -n 60 $RESULTS
LAST=$(tail -n 1 $RESULTS)
if [ $LAST == "OK" ];
then
	EXIT=0
else
	EXIT=1
fi

if [ -d log ]
then
	rm -r log
fi
rm $RESULTS
exit $EXIT

