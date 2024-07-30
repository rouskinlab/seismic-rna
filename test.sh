#!/bin/bash

# Run all unit tests for SEISMIC-RNA.
# If all tests succeed, return 0; otherwise, return 1.

RESULTS=test-results.txt

seismic +test -vv 2> $RESULTS
cat $RESULTS
LAST=$(tail -n 1 $RESULTS)
if [ $LAST == "OK" ];
then
	EXIT=0
else
	EXIT=1
fi

rm -r log
rm $RESULTS
exit $EXIT

