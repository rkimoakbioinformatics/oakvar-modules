#!/bin/bash

for chr in {1..22} X Y; do
    echo "Processing Chromosome $chr..."
    python get_data_all.py -n $chr
done