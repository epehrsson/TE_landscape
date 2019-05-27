#!/bin/bash

# Combines per-sample files listing the number of annotated bases and bases per state into a single file

# States
states=$1

# Samples
header=$2

# Print row names (Total, states)
echo 'Total' > test.txt
cat $states >> test.txt

# Add a column for each sample
for file in ${@:3}; do paste test.txt $file > test; mv test test.txt; done

# Add column names
(printf "\t"; tr "\n" "\t" < $header | sed "s/\t$//"; printf "\n"; cat test.txt) > test
mv test test.txt
echo 'Combined states'
