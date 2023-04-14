#!/bin/bash

set -e

NDUMPS=2
NSAMPLES=100

for i in $(seq 1 $NDUMPS)
do 
	python3 dump_data.py $NSAMPLES
done	
