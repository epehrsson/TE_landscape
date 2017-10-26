#!/bin/bash

states=$1
header=$2

echo 'Total' > test.txt
cat $states >> test.txt
for file in ${@:3}; do paste test.txt $file > test; mv test test.txt; done
(printf "\t"; tr "\n" "\t" < $header | sed "s/\t$//"; printf "\n"; cat test.txt) > test
mv test test.txt
echo 'Combined states'
